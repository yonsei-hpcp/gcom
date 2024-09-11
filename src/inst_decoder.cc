#include <iostream>
#include <fstream>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/serialization.hpp>
#include <filesystem>
#include "worker.h"
#include "inst_decoder.h"

using namespace GCoM;

#define WARP_SIZE 32 // Strictly speaking it is HW specification. But let's consider it HW independent for now since it is implicitly included in SASS trace.
#define TURING_BINARY_VERSION 75
#define VOLTA_BINARY_VERSION 70
const unsigned long long LOCAL_MEM_SIZE_MAX = 1 << 14; // Copied from Accel-sim gpu-simulator/gpgpu-sim/src/abstract_hardware_model.h: Volta max local mem is 16kB

GCoM::SASSDecoder::SASSDecoder(OptionParser opp)
{
    mSASSTracePath = opp.mSASSTracePath;
    mAccelsimTracer = new trace_parser((mSASSTracePath / "kernelslist.g").c_str());
    mCommandlist = mAccelsimTracer->parse_commandlist_file();
    mCurrentGPUCommand = mCommandlist.begin();
}

GCoM::SASSDecoder::~SASSDecoder()
{
    delete mAccelsimTracer;
}

int GCoM::SASSDecoder::ParseKernalTrace()
{
    while (mCurrentGPUCommand != mCommandlist.end() && mCurrentGPUCommand->m_type != command_type::kernel_launch)
    {
        mCurrentGPUCommand++;
    }
    if (mCurrentGPUCommand == mCommandlist.end())
            return -1;
    
    kernel_trace_t *kernel_trace_info = mAccelsimTracer->parse_kernel_info(mCurrentGPUCommand->command_string);
    assert(kernel_trace_info->binary_verion == TURING_BINARY_VERSION || kernel_trace_info->binary_verion == VOLTA_BINARY_VERSION);

    mKernelInfo = KernelInfo();
    unsigned ctaSize = kernel_trace_info->tb_dim_x * kernel_trace_info->tb_dim_y * kernel_trace_info->tb_dim_z;
    mKernelInfo.numCTAs = kernel_trace_info->grid_dim_x * kernel_trace_info->grid_dim_y * kernel_trace_info->grid_dim_z;
    unsigned kernelSize = ctaSize * mKernelInfo.numCTAs;
    unsigned endThread = ctaSize; //number of warps per CTA

    FreeWarps(); // To ensure mWarps has been freed
    for (; endThread <= kernelSize; endThread += ctaSize)
    {
        // for each CTA
        std::vector<std::vector<inst_trace_t> *> threadblockTraces; // insts before decoding
        for (unsigned i = 0; i < ctaSize / WARP_SIZE + ((ctaSize % WARP_SIZE) ? 1: 0); i++)
            threadblockTraces.push_back(new std::vector<inst_trace_t>);
            
        mAccelsimTracer->get_next_threadblock_traces(threadblockTraces,
                                                    kernel_trace_info->trace_verion,
                                                    kernel_trace_info->ifs);
        
        // initial decoding before cache model. TODO: execute in parallel
        unsigned prevNumWarps = mKernelInfo.numWarps;
        mKernelInfo.numWarps += threadblockTraces.size();
        mWarps.resize(mKernelInfo.numWarps);
        for (unsigned i = prevNumWarps; i < mKernelInfo.numWarps; i++)
        {
            // for each new warps in a CTA
            unsigned ctaWarpIdx = i - prevNumWarps;
            mWarps[i].mInsts.resize(threadblockTraces[ctaWarpIdx]->size());
            for (long unsigned int j = 0; j < threadblockTraces[ctaWarpIdx]->size(); j++)
            {
                DecocdeInst(threadblockTraces[ctaWarpIdx]->at(j), *kernel_trace_info, mWarps[i].mInsts[j]);
                mWarps[i].mInsts[j].mWarpPtr = &mWarps[i];
            }
        }

        for (auto vecPtr : threadblockTraces)
            delete vecPtr;
    }

    mAccelsimTracer->kernel_finalizer(kernel_trace_info);
    
    mKernelInfo.sharedMemSize = kernel_trace_info->shmem;
    mKernelInfo.numRegs = kernel_trace_info->nregs;
    if (mKernelInfo.numWarps % mKernelInfo.numCTAs != 0)
    {
        std::cout << "[Error] Number of warps is not multiple of CTA size" << std::endl;
        return -1;
    } else
    {
        mCurrentGPUCommand++;
        return 0;
    }
}

WarpInst &GCoM::SASSDecoder::GetWarpInst(int warpId, int instIdx)
{
    return mWarps[warpId].mInsts[instIdx];
}

Warp &GCoM::SASSDecoder::GetWarp(int warpId)
{
    return mWarps[warpId];
}

int GCoM::SASSDecoder::FreeWarps()
{
    mWarps.clear();
    mWarps.shrink_to_fit();
    return 0;
}

void GCoM::SASSDecoder::ExportWarps(std::filesystem::path exportPath)
{
    std::ofstream ofs(exportPath, std::ios::binary);
    boost::archive::binary_oarchive oa(ofs);
    oa & mWarps;
    ofs.close();
}

// Decode SASS instruction
// Reference: Accelsim bool trace_warp_inst_t::parse_from_trace_struct
int GCoM::SASSDecoder::DecocdeInst(inst_trace_t &instTrace,  kernel_trace_t kernel_trace_info, WarpInst &warpInst)
{
    // decode jobs for the the cache model: decode opcode, active mask, memory space, cache operation type
    warpInst.mDecoded.activeMask = instTrace.mask;
    warpInst.mDecoded.pc = (GCoM::Address) instTrace.m_pc;

    std::vector<std::string> opcode_tokens = instTrace.get_opcode_tokens();
    std::string opcode1 = opcode_tokens[0];

    if (kernel_trace_info.binary_verion == TURING_BINARY_VERSION)
        warpInst.mDecoded.op = TuringOpcodeMap.at(opcode1);
    else if (kernel_trace_info.binary_verion == VOLTA_BINARY_VERSION)
        warpInst.mDecoded.op = VoltaOpcodeMap.at(opcode1);
    // else need to be assert false on higher function
    
    if (warpInst.mDecoded.op == EUArchOp::LOAD_OP || warpInst.mDecoded.op == EUArchOp::STORE_OP)
    {
        if (opcode1 == "LDG" || opcode1 == "LDL") {
            if (opcode1 == "LDG")
                warpInst.mDecoded.space = EMemorySpace::GLOBAL_SPACE;
            else
                warpInst.mDecoded.space = EMemorySpace::LOCAL_SPACE;
            
            // check the cache scope, if its strong GPU, then bypass L1
            if (instTrace.check_opcode_contain(opcode_tokens, "STRONG") &&
                instTrace.check_opcode_contain(opcode_tokens, "GPU"))
            {
                warpInst.mDecoded.cacheOp = ECacheOpType::CACHE_GLOBAL;
            } else
            {
                warpInst.mDecoded.cacheOp = ECacheOpType::CACHE_ALL;
            }
        }
        else if (opcode1 == "STG" || 
                opcode1 == "ATOMG" || opcode1 == "RED" || opcode1 == "ATOM")
        {
            warpInst.mDecoded.space = EMemorySpace::GLOBAL_SPACE;

            if (opcode1 == "STG")
                warpInst.mDecoded.cacheOp = ECacheOpType::CACHE_ALL;
            else
            {
                warpInst.mDecoded.cacheOp = ECacheOpType::CACHE_GLOBAL; // all the atomics should be done at L2
                warpInst.mDecoded.isAtomic = true;
            }
        } 
        else if (opcode1 == "STL")
        {
            warpInst.mDecoded.space = EMemorySpace::LOCAL_SPACE;
            warpInst.mDecoded.cacheOp = ECacheOpType::CACHE_ALL;
        } 
        else if (opcode1 == "LDS" || opcode1 == "STS" || 
                opcode1 == "ATOMS" || opcode1 == "LDSM")
        {
            warpInst.mDecoded.space = EMemorySpace::SHARED_SPACE;
        } 
        else if (opcode1 == "LD" || opcode1 == "ST")
        {
            // resolve generic loads
            if (kernel_trace_info.shmem_base_addr == 0 || kernel_trace_info.local_base_addr == 0)
            {
                // shmem and local addresses are not set
                // assume all the mem reqs are shared by default
                warpInst.mDecoded.space = EMemorySpace::SHARED_SPACE;
            } 
            else
            {
                // check the first active address
                for (unsigned i = 0; i < warpInst.mDecoded.activeMask.size(); ++i)
                {
                    if (warpInst.mDecoded.activeMask.test(i)) 
                    {
                        if (instTrace.memadd_info->addrs[i] >= kernel_trace_info.shmem_base_addr &&
                                instTrace.memadd_info->addrs[i] < kernel_trace_info.local_base_addr)
                            warpInst.mDecoded.space = EMemorySpace::SHARED_SPACE;
                        else if (instTrace.memadd_info->addrs[i] >= kernel_trace_info.local_base_addr &&
                                instTrace.memadd_info->addrs[i] < kernel_trace_info.local_base_addr + LOCAL_MEM_SIZE_MAX) {
                            warpInst.mDecoded.space = EMemorySpace::LOCAL_SPACE;
                            warpInst.mDecoded.cacheOp = ECacheOpType::CACHE_ALL;
                        } 
                        else 
                        {
                            warpInst.mDecoded.space = EMemorySpace::GLOBAL_SPACE;
                            warpInst.mDecoded.cacheOp = ECacheOpType::CACHE_ALL;
                        }
                        break;
                    }
                }
            }
        }
    }

    // copy datas needed for cache model
    // memory coalescing (i.e., accessQ generation) will be handled in cache model
    if (instTrace.memadd_info != NULL)
    {
        warpInst.mDataSize = instTrace.memadd_info->width;
        warpInst.mPerThreadMemAddr.resize(WARP_SIZE);
        for (unsigned i = 0; i < WARP_SIZE; i++)
            warpInst.mPerThreadMemAddr[i] = instTrace.memadd_info->addrs[i];
    }

    // decode jobs to generate interval profile: register dependency, latency, initiation interval(# functional units)
    warpInst.mDecoded.out.resize(instTrace.reg_dsts_num);
    for (unsigned i = 0; i < instTrace.reg_dsts_num; i++)
        warpInst.mDecoded.out[i] = instTrace.reg_dest[i] + 1; // +1 Aligning with Accel-sim where starts from R1
    
    warpInst.mDecoded.in.resize(instTrace.reg_srcs_num);
    for (unsigned i = 0; i < instTrace.reg_srcs_num; i++)
        warpInst.mDecoded.in[i] = instTrace.reg_src[i] + 1; // +1 Aligning with Accel-sim where starts from R1


    ClassifyLatencyInitIntv(opcode1, warpInst);

    // Count the number of active threads for the performance model
    mKernelInfo.numThreadInsts += warpInst.mDecoded.activeMask.count();
    return 0;
}

int GCoM::SASSDecoder::ClassifyLatencyInitIntv(std::string instTraceOpcode, WarpInst &warpInst)
{
    warpInst.mOpLatencyInitIntvType = EOpLatencyInitIntvType::DEFALUT;

    switch (warpInst.mDecoded.op)
    {
    case EUArchOp::ALU_OP:
    case EUArchOp::INTP_OP:
    case EUArchOp::BRANCH_OP:
    case EUArchOp::CALL_OPS:
    case EUArchOp::RET_OPS:
        warpInst.mOpLatencyInitIntvType = EOpLatencyInitIntvType::INT;
        break;
    case EUArchOp::SP_OP:
        warpInst.mOpLatencyInitIntvType = EOpLatencyInitIntvType::FP;
        if (instTraceOpcode == "HADD2" || instTraceOpcode == "HADD2_32I" ||
            instTraceOpcode == "HFMA2" || instTraceOpcode == "HFMA2_32I" ||
            instTraceOpcode == "HMUL2" || instTraceOpcode == "HMUL2_32I" ||
            instTraceOpcode == "HSET2" || instTraceOpcode == "HSETP2")
        {
            warpInst.mOpLatencyInitIntvType = EOpLatencyInitIntvType::FP16;
        }
        break;
    case EUArchOp::DP_OP:
        warpInst.mOpLatencyInitIntvType = EOpLatencyInitIntvType::DP;
        break;
    case EUArchOp::SFU_OP:
        warpInst.mOpLatencyInitIntvType = EOpLatencyInitIntvType::SFU;
        break;
    case EUArchOp::TENSOR_CORE_OP:
        warpInst.mOpLatencyInitIntvType = EOpLatencyInitIntvType::TENSOR;
        break;
    case EUArchOp::SPECIALIZED_UNIT_1_OP:
        warpInst.mOpLatencyInitIntvType = EOpLatencyInitIntvType::SPECIALIZED_UNIT_1;
        break;
    case EUArchOp::SPECIALIZED_UNIT_2_OP:
        warpInst.mOpLatencyInitIntvType = EOpLatencyInitIntvType::SPECIALIZED_UNIT_2;
        break;
    case EUArchOp::SPECIALIZED_UNIT_3_OP:
        warpInst.mOpLatencyInitIntvType = EOpLatencyInitIntvType::SPECIALIZED_UNIT_3;
        break;
    case EUArchOp::SPECIALIZED_UNIT_4_OP:
        warpInst.mOpLatencyInitIntvType = EOpLatencyInitIntvType::SPECIALIZED_UNIT_4;
        break; 
    case EUArchOp::SPECIALIZED_UNIT_5_OP:
        warpInst.mOpLatencyInitIntvType = EOpLatencyInitIntvType::SPECIALIZED_UNIT_5;
        break;
    case EUArchOp::SPECIALIZED_UNIT_6_OP:
        warpInst.mOpLatencyInitIntvType = EOpLatencyInitIntvType::SPECIALIZED_UNIT_6;
        break;
    case EUArchOp::SPECIALIZED_UNIT_7_OP:
        warpInst.mOpLatencyInitIntvType = EOpLatencyInitIntvType::SPECIALIZED_UNIT_7;
        break;
    case EUArchOp::SPECIALIZED_UNIT_8_OP:
        warpInst.mOpLatencyInitIntvType = EOpLatencyInitIntvType::SPECIALIZED_UNIT_8;
        break;
    default:
        break;
    }

    return 0;
}