#include <filesystem>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <bitset>
#include <cmath>
#include <chrono>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/vector.hpp>
#include "worker.h"
#include "interval_model.h"
#include "cache_model.h"
#include "inst_decoder.h"

using namespace GCoM;
using namespace std;

int GCoM::RunGCoMWithSimpleCacheModel(OptionParser opp)
{
    chrono::microseconds totTime = chrono::microseconds(0);
    chrono::microseconds traceReadTime = chrono::microseconds(0);
    chrono::microseconds cacheSimTime = chrono::microseconds(0);
    chrono::microseconds repwarpSelTime = chrono::microseconds(0);
    chrono::microseconds intervalProfilerTime = chrono::microseconds(0);
    chrono::microseconds perfModelTime = chrono::microseconds(0);
    chrono::system_clock::time_point startTime = chrono::system_clock::now();
    

    // run SASS parser identical with accel-sim
    SASSDecoder decoder(opp);
    
    HWConfig hwConfig;
    opp.ReadHWConfig(hwConfig);

    GCoMModel::Result appResult;
    vector<GCoMModel::Result> eachKernelResults;
    unsigned  kernelIdx = 1; // starting from 1
    while (decoder.ParseKernalTrace() != -1) // Parse and decode SASS instructions of a kernel
    {
        // for each kernel

        chrono::system_clock::time_point componentStartTime = chrono::system_clock::now();

        // initialize cache model and run
        // Emulate TB allocation and warp scheduling, access cache
        // input: cacheModel, decoder, hwConfig
        // output:
        //      Update decoder.GetWarpInst().mMemStat decoder.GetWarpInst().mDecoded.accessQ
        //      Update cacheModel.mWarpInstStatistics_PCInstIdx2CacheHitMiss
        CacheModel simpleCache(ECacheModelType::SIMPLE, hwConfig);
        if (simpleCache.SetSharedMemSizeMaxCta(decoder.mKernelInfo, hwConfig) != 0)
            return -1;
        if (RunSimpleCache(&simpleCache, &decoder, hwConfig) !=0)
            return -1;
        
        cacheSimTime += chrono::duration_cast<chrono::microseconds>(chrono::system_clock::now() - componentStartTime);
        componentStartTime = chrono::system_clock::now();

        // representative warp selection in parallel
        GCoMModel performanceModel(hwConfig);
        int representativeWarpID;
        if (performanceModel.SelectRepresentativeWarp(&decoder, opp.mExportWarpSelInputPath, opp.mRepWarpPath, kernelIdx, representativeWarpID) != 0)
            return -1;
        
        Warp representativeWarp = decoder.GetWarp(representativeWarpID);
        // decoder.FreeWarps();

        simpleCache.UpdateGlobalCacheStat(representativeWarp);
        // delete accelsimCache;

        repwarpSelTime += chrono::duration_cast<chrono::microseconds>(chrono::system_clock::now() - componentStartTime);
        componentStartTime = chrono::system_clock::now();

        if (performanceModel.ProfileInterval(representativeWarp))
            return -1;

        intervalProfilerTime += chrono::duration_cast<chrono::microseconds>(chrono::system_clock::now() - componentStartTime);
        componentStartTime = chrono::system_clock::now();

        // run performance model
        GCoMModel::Result kernelResult;
        if (performanceModel.RunPerformanceModel(decoder.mKernelInfo, kernelResult))
            return -1;

        eachKernelResults.push_back(kernelResult);
        appResult = appResult + kernelResult;
        cout << "Kernel " << kernelIdx << " done" << endl;
        kernelIdx ++;

        perfModelTime += chrono::duration_cast<chrono::microseconds>(chrono::system_clock::now() - componentStartTime);
    }

    for (long unsigned int i = 0; i < eachKernelResults.size(); i++)
    {
        cout << "K" << i + 1 << " ";
        cout << eachKernelResults[i] << endl;
    }
    
    cout << "Total ";
    cout << appResult << endl;

    if (DEBUG_LEVEL == 1)
    {
        chrono::system_clock::time_point endTime = chrono::system_clock::now();
        totTime = chrono::duration_cast<chrono::microseconds> (endTime - startTime);
        traceReadTime = totTime - cacheSimTime - repwarpSelTime - intervalProfilerTime - perfModelTime;
        cout << "Time\t,total\t,traceRead\t,cacheSim\t,repwarpSel\t,intervalProfiler\t,perfModel" << endl;
        cout << "uSeconds\t," << totTime.count() << "\t," << traceReadTime.count() << "\t," << cacheSimTime.count() << "\t," << repwarpSelTime.count() << "\t," << intervalProfilerTime.count() << "\t," << perfModelTime.count() << endl;
    }
    return 0;
}

int GCoMModel::SelectRepresentativeWarp(SASSDecoder *decoder, fs::path warpSelInputPath, fs::path repWarpPath, unsigned kernelIdx, int &representativeWarpID)
{
    // Export decoder.mWarps for external warp selector input
    if (warpSelInputPath != "")
        decoder->ExportWarps(warpSelInputPath);

    // read representative warp from the external warp selector
    ifstream ifs(repWarpPath, ios::binary);
    if(ifs.is_open() != true) {
        cout << "[Error] representative warp file not ready" << endl;
        return -1;
    }
    boost::archive::binary_iarchive ia(ifs);
    
    vector<KernelRepWarp> kernelRepWarpVector;
    ia & kernelRepWarpVector;

    KernelRepWarp kernelRepWarp;
    try 
    {
        kernelRepWarp = kernelRepWarpVector[kernelIdx - 1];
        assert(kernelRepWarp.kernelNumber == (int) kernelIdx);
    } catch (const char* msg) {
        cout << "[Error] representative warp not ready" << endl;
        cout << msg << endl;
        return -1;
    }

    representativeWarpID = kernelRepWarp.repWarpIdx;
    return 0;
}

int GCoMModel::ProfileInterval(Warp &representativeWarp)
{
    map<unsigned, DepTableEntry> dependencyTable; // {regID: doneCycle}

    unsigned prevIssueCycle = 0;
    IntervalProfile::IntervalStat intervalStat;
    double sumIntvWarpInsts = 0;
    double sumStallCycle = 0;

    for (unsigned instIdx = 0; instIdx < representativeWarp.mInsts.size(); instIdx++)
    {
        WarpInst &warpInst = representativeWarp.mInsts[instIdx];

        // check source operand dependency and decide issue cycle
        unsigned curIssueCycle;
        bool hasMemDep = false;
        CheckSrcDep(warpInst, dependencyTable, prevIssueCycle, curIssueCycle, hasMemDep);
        intervalStat.hasMemInst = intervalStat.hasMemInst || hasMemDep;

        // if issue cycle is not continuous, new interval starts. Store the previous interval
        if (curIssueCycle != prevIssueCycle + 1)
        {
            intervalStat.stallCycle = curIssueCycle - prevIssueCycle - 1;
            mIntervalProfile.push_back(intervalStat);
            
            sumIntvWarpInsts += intervalStat.nWarpInst;
            sumStallCycle += intervalStat.stallCycle;

            intervalStat = IntervalProfile::IntervalStat();
        }
        prevIssueCycle = curIssueCycle;

        // decide latency of warpInst and collect interval statistics
        if (DecideInstLatencyCollectIntervalStat(warpInst, instIdx, intervalStat) == -1)
            return -1;

        // update dependency table
        for (unsigned outRegID : warpInst.mDecoded.out)
        {
            DepTableEntry entry;
            entry.doneCycle = curIssueCycle + warpInst.mDecoded.latency;
            entry.warpInst = &warpInst;
            dependencyTable[outRegID] = entry;
        }
    }

    // store the last interval
    // modeling cycle waiting for last instruction to end
    // it does not consider latency hiding effect of another CTA, but this will not be dominant when CTA # is large
    unsigned doneCycle = 0;
    for (unsigned instIdx = representativeWarp.mInsts.size() - intervalStat.nWarpInst; instIdx < representativeWarp.mInsts.size(); instIdx++)
    {
        WarpInst &warpInst = representativeWarp.mInsts[instIdx];
        unsigned issueCycle = prevIssueCycle - representativeWarp.mInsts.size() + instIdx + 1;
        doneCycle = max(doneCycle, issueCycle + warpInst.mDecoded.latency);
    }
    intervalStat.stallCycle = doneCycle - prevIssueCycle - 1;
    mIntervalProfile.push_back(intervalStat);

    sumIntvWarpInsts += intervalStat.nWarpInst;
    sumStallCycle += intervalStat.stallCycle;
    mIntervalProfile.mIssueProb = sumIntvWarpInsts / (sumIntvWarpInsts + sumStallCycle);
    mIntervalProfile.mAvgIntvWarpInsts = sumIntvWarpInsts / mIntervalProfile.len();

    return 0;
}

GCoMModel::Result GCoMModel::Result::operator+(const GCoMModel::Result &other) const
{
    Result result;
    result.numThreadInst = numThreadInst + other.numThreadInst;
    result.numCycle = numCycle + other.numCycle;
    result.cpi = (double) result.numCycle / result.numThreadInst;
    result.cpiBase = (cpiBase * numThreadInst + other.cpiBase * other.numThreadInst) / result.numThreadInst;
    result.cpiComData = (cpiComData * numThreadInst + other.cpiComData * other.numThreadInst) / result.numThreadInst;
    result.cpiComStruct = (cpiComStruct * numThreadInst + other.cpiComStruct * other.numThreadInst) / result.numThreadInst;
    result.cpiMemData = (cpiMemData * numThreadInst + other.cpiMemData * other.numThreadInst) / result.numThreadInst;
    result.cpiMemStruct = (cpiMemStruct * numThreadInst + other.cpiMemStruct * other.numThreadInst) / result.numThreadInst;
    result.cpiIdle = (cpiIdle * numThreadInst + other.cpiIdle * other.numThreadInst) / result.numThreadInst;
    return result;
}

std::ostream& operator<<(std::ostream &os, const GCoMModel::Result &result)
{
    os << "Result\t,numThreadInst\t,numCycle\t,cpi\t,base\t,comData\t,comStruct\t,memData\t,memStruct\t,idle" << endl;
    os << "Value\t," << result.numThreadInst << "\t," << result.numCycle << "\t," << result.cpi << "\t," << result.cpiBase << "\t," << result.cpiComData << "\t," << result.cpiComStruct << "\t," << result.cpiMemData << "\t," << result.cpiMemStruct << "\t," << result.cpiIdle;

    return os;
}

int GCoMModel::RunPerformanceModel(KernelInfo kernelInfo, Result &result)
{
    CycleBreakdown GPUCycleBreakdown;
    unsigned numActiveSM = min(mHWConfig.numSMs, kernelInfo.numCTAs);
    unsigned numWarpsPerCTA = kernelInfo.numWarps / kernelInfo.numCTAs;
    unsigned maxConcurrentCTAPerSM;
    if (CalculateMaxConcurrentCTAPerSM(kernelInfo, mHWConfig, maxConcurrentCTAPerSM) != 0)
        return -1;

    int error = 0;

    // SMs running more CTAs when kernelInfo.numCTAs % numActiveSM != 0
    unsigned numMoreCTASMs = kernelInfo.numCTAs % numActiveSM;
    CycleBreakdown longSMCycleBreakdown;
    if (numMoreCTASMs > 0)
    {
        unsigned allocedCTAPerSM = kernelInfo.numCTAs / numActiveSM + 1;
        error += RunPerSMModel(allocedCTAPerSM, maxConcurrentCTAPerSM, numWarpsPerCTA, numActiveSM, longSMCycleBreakdown);
        GPUCycleBreakdown = GPUCycleBreakdown + (longSMCycleBreakdown / numActiveSM * numMoreCTASMs); // Average each SM
    }

    // SMs running less CTAs
    unsigned allocedCTAPerSM = kernelInfo.numCTAs / numActiveSM;
    unsigned numLessCTASMs = numActiveSM - numMoreCTASMs;
    CycleBreakdown smCycleBreakdown;
    error += RunPerSMModel(allocedCTAPerSM, maxConcurrentCTAPerSM, numWarpsPerCTA, numActiveSM, smCycleBreakdown);
    GPUCycleBreakdown = GPUCycleBreakdown + (smCycleBreakdown / numActiveSM * numLessCTASMs); // Average each SM

    if (error != 0)
        return -1;

    // SM load imbalance model
    // Unlike the ISCA'22 paper, GCoM code did not attribute inactive SM cycles as idle cycles
    double longestCycle = max(longSMCycleBreakdown.tot, smCycleBreakdown.tot);
    GPUCycleBreakdown.tot = longestCycle;
    GPUCycleBreakdown.idle += (longestCycle - longSMCycleBreakdown.tot) * numMoreCTASMs / numActiveSM;
    GPUCycleBreakdown.idle += (longestCycle - smCycleBreakdown.tot) * numLessCTASMs / numActiveSM;

    GPUCycleBreakdown.tot += 1;
    GPUCycleBreakdown.base += 1;

    assert(abs(GPUCycleBreakdown.tot - GPUCycleBreakdown.base - GPUCycleBreakdown.comData - GPUCycleBreakdown.comStruct
            - GPUCycleBreakdown.memData - GPUCycleBreakdown.memStruct - GPUCycleBreakdown.idle) < 10);
    assert(GPUCycleBreakdown.memData + 1e-6 >= GPUCycleBreakdown.memDataQueue + GPUCycleBreakdown.memDataScbLast);
    assert(abs(GPUCycleBreakdown.memDataQueue - GPUCycleBreakdown.memDataQueueNoC - GPUCycleBreakdown.memDataQueueDRAM) < 1e-6);
    assert(abs(GPUCycleBreakdown.memStruct - GPUCycleBreakdown.memStructL1Bank - GPUCycleBreakdown.memStructMSHR) < 1e-6);

    result.numThreadInst = kernelInfo.numThreadInsts;
    result.numCycle = GPUCycleBreakdown.tot;
    CycleBreakdown cpi = GPUCycleBreakdown / kernelInfo.numThreadInsts;
    result.cpi = cpi.tot;
    result.cpiBase = cpi.base;
    result.cpiComData = cpi.comData;
    result.cpiComStruct = cpi.comStruct;
    result.cpiMemData = cpi.memData;
    result.cpiMemStruct = cpi.memStruct;
    result.cpiIdle = cpi.idle;

    return 0;
}

GCoMModel::CycleBreakdown GCoMModel::CycleBreakdown::operator+(const GCoMModel::CycleBreakdown &other) const
{
    CycleBreakdown result;
    result.tot = tot + other.tot;
    result.base = base + other.base;
    result.comData = comData + other.comData;
    result.comStruct = comStruct + other.comStruct;
    result.memData = memData + other.memData;
    result.memStruct = memStruct + other.memStruct;
    result.idle = idle + other.idle;
    result.memDataScbLast = memDataScbLast + other.memDataScbLast;
    result.memDataQueue = memDataQueue + other.memDataQueue;
    result.memDataQueueNoC = memDataQueueNoC + other.memDataQueueNoC;
    result.memDataQueueDRAM = memDataQueueDRAM + other.memDataQueueDRAM;
    result.memStructMSHR = memStructMSHR + other.memStructMSHR;
    result.memStructL1Bank = memStructL1Bank + other.memStructL1Bank;
    return result;
}

GCoMModel::CycleBreakdown GCoMModel::CycleBreakdown::operator/(const unsigned long long &divisor) const
{
    CycleBreakdown result;
    result.tot = tot / divisor;
    result.base = base / divisor;
    result.comData = comData / divisor;
    result.comStruct = comStruct / divisor;
    result.memData = memData / divisor;
    result.memStruct = memStruct / divisor;
    result.idle = idle / divisor;
    result.memDataScbLast = memDataScbLast / divisor;
    result.memDataQueue = memDataQueue / divisor;
    result.memDataQueueNoC = memDataQueueNoC / divisor;
    result.memDataQueueDRAM = memDataQueueDRAM / divisor;
    result.memStructMSHR = memStructMSHR / divisor;
    result.memStructL1Bank = memStructL1Bank / divisor;
    return result;
}

GCoMModel::CycleBreakdown GCoMModel::CycleBreakdown::operator*(const unsigned &multiplier) const
{
    CycleBreakdown result;
    result.tot = tot * multiplier;
    result.base = base * multiplier;
    result.comData = comData * multiplier;
    result.comStruct = comStruct * multiplier;
    result.memData = memData * multiplier;
    result.memStruct = memStruct * multiplier;
    result.idle = idle * multiplier;
    result.memDataScbLast = memDataScbLast * multiplier;
    result.memDataQueue = memDataQueue * multiplier;
    result.memDataQueueNoC = memDataQueueNoC * multiplier;
    result.memDataQueueDRAM = memDataQueueDRAM * multiplier;
    result.memStructMSHR = memStructMSHR * multiplier;
    result.memStructL1Bank = memStructL1Bank * multiplier;
    return result;
}

int GCoMModel::CheckSrcDep(WarpInst &warpInst, std::map<unsigned, DepTableEntry> &dependencyTable, unsigned prevIssueCycle,
                unsigned &curIssueCycle, bool &hasMemDep)
{
    curIssueCycle = prevIssueCycle + 1;
    hasMemDep = false;
    for (unsigned srcRegID : warpInst.mDecoded.in)
    {
        auto it = dependencyTable.find(srcRegID);
        if (it != dependencyTable.end())
        {
            curIssueCycle = max(curIssueCycle, it->second.doneCycle);

            WarpInst *srcWarpInstPtr = it->second.warpInst;
            EUArchOp op = srcWarpInstPtr->mDecoded.op;
            EMemorySpace memSpace = srcWarpInstPtr->mDecoded.space;
            if ((curIssueCycle != prevIssueCycle + 1) &&
                (op == EUArchOp::LOAD_OP || op == EUArchOp::STORE_OP) 
                && (memSpace == EMemorySpace::GLOBAL_SPACE || memSpace == EMemorySpace::LOCAL_SPACE))
                hasMemDep = true;
        }
    }

    return 0;
}

int GCoMModel::DecideInstLatencyCollectIntervalStat(WarpInst &warpInst, unsigned instIdx, IntervalProfile::IntervalStat &intervalStat)
{
    intervalStat.nWarpInst++;
    unsigned numActiveThreadInst = warpInst.mDecoded.activeMask.count();
    intervalStat.nActiveThreadInst += numActiveThreadInst;

    if (ReadLatencyInitIntv(warpInst) == -1)
        return -1;
    

    EUArchOp op = warpInst.mDecoded.op;
    EMemorySpace memSpace = warpInst.mDecoded.space;
    if ((op == EUArchOp::LOAD_OP || op == EUArchOp::STORE_OP) 
                && (memSpace == EMemorySpace::GLOBAL_SPACE || memSpace == EMemorySpace::LOCAL_SPACE)
                && numActiveThreadInst != 0)
    {
        Warp::WarpInstCacheHitMissCount &globalCacheStat = warpInst.mWarpPtr->mGlobalCacheStat[instIdx];
        if ((globalCacheStat.l1MissWarp + globalCacheStat.l1HitWarp == 0)
                || (globalCacheStat.l2MissWarp + globalCacheStat.l2HitWarp == 0))
        {
            cout << "[Error] invalid global cache stat" << endl;
            return -1;
        }

        double avgLatency = 0;
        if (op == EUArchOp::LOAD_OP)
        {
            double l1MissRate = (double) globalCacheStat.l1MissWarp / (globalCacheStat.l1MissWarp + globalCacheStat.l1HitWarp);
            double l2MissRate = (double) globalCacheStat.l2MissWarp / (globalCacheStat.l2MissWarp + globalCacheStat.l2HitWarp);

            avgLatency += mHWConfig.l1DLatency;
            avgLatency += mHWConfig.L2Latency * l1MissRate
                    + mHWConfig.dramLatency * l1MissRate * l2MissRate;

            intervalStat.l1Rmiss += warpInst.mMemStat.coalescedL1Miss;
            intervalStat.l2Rhit += globalCacheStat.l2HitWarp;
            intervalStat.l2Rmiss += globalCacheStat.l2MissWarp;

            warpInst.mDecoded.latency = avgLatency + mHWConfig.loadPipeline;
        }
        else // op == EUArchOp::STORE_OP
        {
            if ((mHWConfig.l1DWritePolicy == EWritePolicy::LOCAL_WB_GLOBAL_WT && mHWConfig.l1DWriteAllocatePolicy == EWriteAllocatePolicy::NO_WRITE_ALLOCATE)
                    || (mHWConfig.l1DWritePolicy == EWritePolicy::WRITE_THROUGH && mHWConfig.l1DWriteAllocatePolicy == EWriteAllocatePolicy::NO_WRITE_ALLOCATE)
                    || (mHWConfig.l1DWritePolicy == EWritePolicy::WRITE_THROUGH && mHWConfig.l1DWriteAllocatePolicy == EWriteAllocatePolicy::LAZY_FETCH_ON_READ))
            {
                avgLatency += mHWConfig.l1DLatency + mHWConfig.L2Latency;
            }
            else
            {
                cout << "[Error] Unsupported L1 cache write policy" << endl;
                return -1;
            }

            if ((mHWConfig.L2WritePolicy == EWritePolicy::WRITE_BACK && mHWConfig.L2WriteAllocatePolicy == EWriteAllocatePolicy::LAZY_FETCH_ON_READ)
                    || (mHWConfig.L2WritePolicy == EWritePolicy::WRITE_BACK && mHWConfig.L2WriteAllocatePolicy == EWriteAllocatePolicy::FETCH_ON_WRITE))
            {
                avgLatency += 0;
            }
            else if (mHWConfig.L2WritePolicy == EWritePolicy::LOCAL_WB_GLOBAL_WT && mHWConfig.L2WriteAllocatePolicy == EWriteAllocatePolicy::NO_WRITE_ALLOCATE)
            {
                avgLatency += mHWConfig.dramLatency;
            }
            else
            {
                cout << "[Error] Unsupported L2 cache write policy" << endl;
                return -1;
            }

            intervalStat.l1Wmiss += warpInst.mMemStat.l1Miss;
            intervalStat.l2Whit += globalCacheStat.l2HitWarp;
            intervalStat.l2Wmiss += globalCacheStat.l2MissWarp;

            warpInst.mDecoded.latency = avgLatency + mHWConfig.storePipeline;
        }
        
        // L1 bank access latency
        unsigned bankAccessLatency = 1;
        bitset<64> usedBank{0x0};
        for (MemAccess memAccess : warpInst.mDecoded.accessQ)
        {
            unsigned bankId = HashAddress(memAccess.mAddr, mHWConfig.l1DnBank, LogB2(mHWConfig.l1DBankInterleaveByte), mHWConfig.l1DBankHashFunction);
            if (!usedBank.test(bankId))
            {
                usedBank.set(bankId);
            }
            else
            {
                usedBank.reset();
                bankAccessLatency++;
                usedBank.set(bankId);
            }
        }
        warpInst.mDecoded.initiationInterval = bankAccessLatency;
        warpInst.mDecoded.latency += bankAccessLatency;

        unsigned regReadLatency;
        ModelRegReadLatency(warpInst, regReadLatency);
        warpInst.mDecoded.latency += regReadLatency;
    } 
    else if ((op == EUArchOp::LOAD_OP || op == EUArchOp::STORE_OP) 
                && memSpace == EMemorySpace::SHARED_SPACE
                && numActiveThreadInst != 0)
    {
        // Simple shared memory model
        // ignore bank conflict, multiple access for simplicity
        warpInst.mDecoded.latency = mHWConfig.sharedMemLatency;
    }
    else if (numActiveThreadInst != 0)
    {
        warpInst.mDecoded.latency += mHWConfig.computePipeline;

        unsigned regReadLatency;
        ModelRegReadLatency(warpInst, regReadLatency);
        warpInst.mDecoded.latency += regReadLatency;
    }
    else // numActiveThreadInst == 0
    {
        warpInst.mDecoded.latency = mHWConfig.computePipeline;
    }

    auto it = TuringOpFUMap.find(op);
    if (it == TuringOpFUMap.end())
    {
        cout << "[Error] Unsupported opcode in TuringOpFUMap: " <<  static_cast<int>(op) << endl;
        return -1;
    }
    intervalStat.dispatch[it->second] += warpInst.mDecoded.initiationInterval;

    return 0;
}

int GCoMModel::ReadLatencyInitIntv(WarpInst &warpInst)
{
    switch (warpInst.mOpLatencyInitIntvType)
    {
    case EOpLatencyInitIntvType::DEFALUT:
        warpInst.mDecoded.latency = 1;
        warpInst.mDecoded.initiationInterval = 1;
        break;
    case EOpLatencyInitIntvType::INT:
        warpInst.mDecoded.latency = mHWConfig.intLatency;
        warpInst.mDecoded.initiationInterval = mHWConfig.intInitIntv;
        break;
    case EOpLatencyInitIntvType::FP:
        warpInst.mDecoded.latency = mHWConfig.fpLatency;
        warpInst.mDecoded.initiationInterval = mHWConfig.fpInitIntv;
        break;
    case EOpLatencyInitIntvType::FP16:
        warpInst.mDecoded.latency = mHWConfig.fpLatency;
        warpInst.mDecoded.initiationInterval = max(mHWConfig.fpInitIntv / 2, (unsigned) 1);
        break;
    case EOpLatencyInitIntvType::DP:
        warpInst.mDecoded.latency = mHWConfig.dpLatency;
        warpInst.mDecoded.initiationInterval = mHWConfig.dpInitIntv;
        break;
    case EOpLatencyInitIntvType::SFU:
        warpInst.mDecoded.latency = mHWConfig.sfuLatency;
        warpInst.mDecoded.initiationInterval = mHWConfig.sfuInitIntv;
        break;
    case EOpLatencyInitIntvType::SPECIALIZED_UNIT_1:
        warpInst.mDecoded.latency = mHWConfig.specializedUnit1Latency;
        warpInst.mDecoded.initiationInterval = mHWConfig.specializedUnit1IntiIntv;
        break;
    case EOpLatencyInitIntvType::SPECIALIZED_UNIT_2:
        warpInst.mDecoded.latency = mHWConfig.specializedUnit2Latency;
        warpInst.mDecoded.initiationInterval = mHWConfig.specializedUnit2IntiIntv;
        break;
    case EOpLatencyInitIntvType::SPECIALIZED_UNIT_3:
        warpInst.mDecoded.latency = mHWConfig.specializedUnit3Latency;
        warpInst.mDecoded.initiationInterval = mHWConfig.specializedUnit3IntiIntv;
        break;
    case EOpLatencyInitIntvType::SPECIALIZED_UNIT_4:
        warpInst.mDecoded.latency = mHWConfig.specializedUnit4Latency;
        warpInst.mDecoded.initiationInterval = mHWConfig.specializedUnit4IntiIntv;
        break;
    case EOpLatencyInitIntvType::SPECIALIZED_UNIT_5:
        warpInst.mDecoded.latency = mHWConfig.specializedUnit5Latency;
        warpInst.mDecoded.initiationInterval = mHWConfig.specializedUnit5IntiIntv;
        break;
    case EOpLatencyInitIntvType::SPECIALIZED_UNIT_6:
        warpInst.mDecoded.latency = mHWConfig.specializedUnit6Latency;
        warpInst.mDecoded.initiationInterval = mHWConfig.specializedUnit6IntiIntv;
        break;
    case EOpLatencyInitIntvType::SPECIALIZED_UNIT_7:
        warpInst.mDecoded.latency = mHWConfig.specializedUnit7Latency;
        warpInst.mDecoded.initiationInterval = mHWConfig.specializedUnit7IntiIntv;
        break;
    case EOpLatencyInitIntvType::SPECIALIZED_UNIT_8:
        warpInst.mDecoded.latency = mHWConfig.specializedUnit8Latency;
        warpInst.mDecoded.initiationInterval = mHWConfig.specializedUnit8IntiIntv;
        break;
    default:
        cout << "[Error] Unsupported operation latency and initiation interval type" << endl;
        return -1;
    }

    return 0;
}

void GCoMModel::ModelRegReadLatency(WarpInst &warpInst, unsigned &regReadLatency)
{
    unsigned minOCStep = (warpInst.mDecoded.in.size() ? 1 : 0); // Assuming minimal step for operand collection
    minOCStep += 2; // Static latency of operand collector from Accel-sim
    regReadLatency = ceil((double) minOCStep / mHWConfig.regFileReadThroughput);
}

int GCoMModel::RunPerSMModel(unsigned allocedCTAPerSM, unsigned maxConcurrentCTAPerSM, unsigned numWarpsPerCTA, unsigned numActiveSM, CycleBreakdown &smCycleBreakdown)
{
    int error = 0;
    unsigned concurrentCTAs = min(maxConcurrentCTAPerSM, allocedCTAPerSM);
    unsigned concurrentWarps = concurrentCTAs * numWarpsPerCTA;
    
    for (unsigned intervalIdx = 0; intervalIdx < mIntervalProfile.len(); intervalIdx++)
    {
        error += RunPerIntervalModel(intervalIdx, concurrentWarps, numActiveSM, smCycleBreakdown);
        if (DEBUG_LEVEL == 10)
            cout << "[" << intervalIdx << "] " << smCycleBreakdown.tot << " " << smCycleBreakdown.base << " " << smCycleBreakdown.comData << " " << smCycleBreakdown.comStruct << " " << smCycleBreakdown.memData << " " << smCycleBreakdown.memStruct << " " << smCycleBreakdown.idle << endl;
    }
    unsigned numRepeat = allocedCTAPerSM / concurrentCTAs;
    smCycleBreakdown = smCycleBreakdown * numRepeat;

    allocedCTAPerSM -= concurrentCTAs * numRepeat;

    if (allocedCTAPerSM > 0)
    {
        // Run the remaining CTAs
        concurrentCTAs = allocedCTAPerSM;
        concurrentWarps = concurrentCTAs * numWarpsPerCTA;
        for (unsigned intervalIdx = 0; intervalIdx < mIntervalProfile.len(); intervalIdx++)
             error += RunPerIntervalModel(intervalIdx, concurrentWarps, numActiveSM, smCycleBreakdown);
    }

    if (error != 0)
        return -1;
    else
        return 0;
}

int GCoMModel::RunPerIntervalModel(unsigned intervalIdx, unsigned smConcurrentWarps, unsigned numActiveSM, CycleBreakdown &smCycleBreakdown)
{
    CycleBreakdown intvCycleBreakdown;
    int error = 0;

    bool isLastInterval = (intervalIdx == mIntervalProfile.len() - 1);
    error += RunMultiSubCoreModel(mIntervalProfile[intervalIdx], isLastInterval, mHWConfig.numSubcorePerSM,
            mHWConfig.subcoreIssueRate, smConcurrentWarps, mHWConfig.wrapSchedulePolicy, intvCycleBreakdown);

    error += RunIntraSMContentionModel(mIntervalProfile[intervalIdx], mHWConfig.numSubcorePerSM, 
            smConcurrentWarps, mHWConfig.operandCollectorQue, mHWConfig.fuQue, intvCycleBreakdown);

    error += RunMemoryContentionModel(mIntervalProfile[intervalIdx], smConcurrentWarps, numActiveSM,
            intvCycleBreakdown);

    smCycleBreakdown = smCycleBreakdown + intvCycleBreakdown;

    if (error != 0)
        return -1;
    else
        return 0;
}

int GCoMModel::RunMultiSubCoreModel(IntervalProfile::IntervalStat intervalStat, bool isLast, unsigned numSubcore, 
        double issueRate, unsigned smConcurrentWarps, string warpSchedulePolicy,
        CycleBreakdown &intvCycleBreakdown)
{
    int scW = (int) ceil((double) smConcurrentWarps / numSubcore); // number of concurrent warps of a subcore
    int lessScW = max(scW - 1, 0); // number of concurrent warps of a subcore that runs less warps than others
    int nSc = smConcurrentWarps  % numSubcore; // number of subcores that runs more warps than others
    if (nSc == 0) 
        nSc = numSubcore;
    int lessNSc = numSubcore - nSc; // number of subcores that runs less warps than others

    double othBase = 0, lessOthBase = 0, issuedInstnInStall = 0, lessIssuedInstnInStall = 0, dataStallSub = 0, lessDataStallSub = 0;

    // latency hiding in subcores (each subcore has a warp scheduler)
    if (RunSubcoreDataStallModel(intervalStat, warpSchedulePolicy, issueRate, (double) max(scW - 1, 0), 
            othBase, dataStallSub, issuedInstnInStall) != 0)
        return -1;

    // latency hiding in subcores that runs less warps than others
    if (lessNSc > 0 && lessScW > 0)
    {
        if (RunSubcoreDataStallModel(intervalStat, warpSchedulePolicy, issueRate, (double) max(lessScW - 1, 0), 
                lessOthBase, lessDataStallSub, lessIssuedInstnInStall) != 0)
            return -1;
    }

    // average all subcores
    double base = (intervalStat.nWarpInst / issueRate + othBase) * nSc;
    if (lessScW > 0)
        base += (intervalStat.nWarpInst / issueRate + lessOthBase) * lessNSc;
    base /= numSubcore;
    double dataStall = (dataStallSub * nSc + lessDataStallSub * lessNSc) / numSubcore;

    // intra-SM load imbalance model
    double idle = 0;
    if (lessNSc > 0 && lessScW > 0)
    {
        idle = (othBase - issuedInstnInStall / issueRate); // non-overlapped instn of subcore
        idle -= (lessOthBase - lessIssuedInstnInStall / issueRate); // non-overlapped instn of subcore that runs less warps
        idle = idle * lessNSc / numSubcore; // average of subcores
    }
    else if (lessNSc > 0)
    {
        // When some subcores are inactive, they are idle for base + stall cycle of the other subcores
        idle = intervalStat.nWarpInst / issueRate + othBase + dataStallSub;
        idle = idle * lessNSc / numSubcore; // average of subcores
    }

    // last instruction latency
    double memDataScbLast = 0;
    if (isLast)
    {
        memDataScbLast = dataStall + base; // may not add base in the future
        dataStall += base; // may deprecated in the future
        if (lessNSc > 0 && lessScW > 0)
            idle += (lessDataStallSub - dataStallSub) * lessNSc / numSubcore;
        else if (lessNSc > 0)
            idle += (intervalStat.nWarpInst / issueRate + othBase) * lessNSc / numSubcore; // may deprecated in the future
    }

    // update cycle breakdown
    // classify dataStall
    if (intervalStat.hasMemInst || memDataScbLast > 0)
        intvCycleBreakdown.memData += dataStall;
    else
        intvCycleBreakdown.comData += dataStall;
    intvCycleBreakdown.memDataScbLast += memDataScbLast;

    intvCycleBreakdown.base += base;
    intvCycleBreakdown.idle += idle;
    intvCycleBreakdown.tot += base + dataStall + idle;

    return 0;
}

int GCoMModel::RunSubcoreDataStallModel(IntervalProfile::IntervalStat intervalStat, std::string warpSchedulePolicy,
                double issueRate, double numOthW, double &othBase, double &dataStallSub, double &issuedInstnInStall)
{
    if (warpSchedulePolicy == "gto") 
    {
        othBase = numOthW * mIntervalProfile.mAvgIntvWarpInsts / issueRate;
        issuedInstnInStall = min(intervalStat.stallCycle * mIntervalProfile.mIssueProb, 1.0) * numOthW * mIntervalProfile.mAvgIntvWarpInsts;
    }
    else if (warpSchedulePolicy == "lrr")
    {
        othBase = numOthW * intervalStat.nWarpInst / issueRate;
        issuedInstnInStall = numOthW * mIntervalProfile.mIssueProb; // Assume numOthW * issue_prob warps are ready
    }
    else
    {
        cout << "[Error] Unsupported warp schedule policy" << endl;
        return -1;
    }
    dataStallSub = max(intervalStat.stallCycle - issuedInstnInStall / issueRate, 0.0);

    return 0;
}

int GCoMModel::RunIntraSMContentionModel(IntervalProfile::IntervalStat intervalStat, unsigned numSubcore, 
                unsigned smConcurrentWarps, unsigned operandCollectorQue, unsigned fuQue, CycleBreakdown &intvCycleBreakdown)
{
    // Assumes LD/ST unit shared by all subcores and compute functional units are not shared
    // Assumes subcore issue rate is 1

    int scW = (int) ceil((double) smConcurrentWarps / numSubcore); // number of concurrent warps of a subcore
    int lessScW = max(scW - 1, 0); // number of concurrent warps of a subcore that runs less warps than others
    int nSc = smConcurrentWarps  % numSubcore; // number of subcores that runs more warps than others
    if (nSc == 0) 
        nSc = numSubcore;

    double nActSc; // number of active subcores
    if (lessScW > 0)
        nActSc = numSubcore;
    else
        nActSc = nSc;

    double memDispatch = 0, comDispatch = 0;
    double minDispatch = intervalStat.nWarpInst * smConcurrentWarps / nActSc;
    double maxDispatch = minDispatch;
    int nFU = 0; // number of functional units
    for (const pair<EFUType, unsigned> &it : intervalStat.dispatch)
    {
        EFUType fuType = it.first;
        unsigned sumIntvDptchCyc = it.second; // sum of dispatch cycles of instructions in the interval considering # of functional unit lanes
        double dispatch;
        if (fuType == EFUType::MEM)
        {
            dispatch = sumIntvDptchCyc * smConcurrentWarps;
            memDispatch += dispatch;
        }
        else
        {
            dispatch = sumIntvDptchCyc * smConcurrentWarps / nActSc;
            comDispatch += dispatch;
        }
        if (dispatch > 0)
            nFU++;
        if (dispatch > maxDispatch)
            maxDispatch = dispatch;
    }
    double intraSMstrcutStall = maxDispatch - minDispatch;

    // 2*FUs + #OC(2) = # of pipeline queue between warp scheduler and execution unit
    // # interval inst - # que inst's dispatch delay cannot be hiden.
    int nQue = fuQue*nFU + operandCollectorQue / numSubcore;
    double unoverlap = max( (double) intervalStat.nWarpInst * smConcurrentWarps - nQue, 0.0);
    double overlap = (double) intervalStat.nWarpInst * smConcurrentWarps - unoverlap;
    intraSMstrcutStall = intraSMstrcutStall * unoverlap / (unoverlap + overlap)
            + max(intraSMstrcutStall * overlap / (unoverlap + overlap) - intervalStat.stallCycle, 0.0);
    
    // classify intra-SM structural stall to memory and compute considering co-existing memory and compute warp instructions
    double denom = max(memDispatch - minDispatch, 0.0) + max(comDispatch - minDispatch, 0.0);
    denom = denom != 0.0 ? denom : 1.0;
    intvCycleBreakdown.memStruct += intraSMstrcutStall * max(memDispatch - minDispatch, 0.0) / denom;
    intvCycleBreakdown.memStructL1Bank += intraSMstrcutStall * max(memDispatch - minDispatch, 0.0) / denom;
    intvCycleBreakdown.comStruct += intraSMstrcutStall * max(comDispatch - minDispatch, 0.0) / denom;
    intvCycleBreakdown.tot += intraSMstrcutStall;

    return 0;
}

int GCoMModel::RunMemoryContentionModel(IntervalProfile::IntervalStat intervalStat, unsigned smConcurrentWarps, unsigned numActiveSM,
                CycleBreakdown &intvCycleBreakdown)
{
    // based on MDM memory contention model
    
    double nocSerLat = mHWConfig.coreFreq * mHWConfig.l1DLineSize / mHWConfig.nSectorPerLine / mHWConfig.maxNoCBW; // NoC service latency in core cycles
    double dramSerLat = mHWConfig.coreFreq * mHWConfig.l1DLineSize / mHWConfig.nSectorPerLine / mHWConfig.maxDRAMBW; // DRAM service latency in core cycles
    double l2RMissRate = intervalStat.l2Rmiss + intervalStat.l2Rhit ? (double) intervalStat.l2Rmiss / (intervalStat.l2Rmiss + intervalStat.l2Rhit) : 0.0;
    double l2WMissRate = intervalStat.l2Wmiss + intervalStat.l2Whit ? (double) intervalStat.l2Wmiss / (intervalStat.l2Wmiss + intervalStat.l2Whit) : 0.0;

    bool isMD; // is memory divergence interval
    if (intervalStat.l1Rmiss * smConcurrentWarps > mHWConfig.nMSHR)
        isMD = true;
    else
        isMD = false;
    double rMiss = min(intervalStat.l1Rmiss * smConcurrentWarps, mHWConfig.nMSHR); // concurrent read miss
    double wMiss = intervalStat.l1Wmiss * smConcurrentWarps; // concurrent write miss

    double missDram; // concurrent miss that goes to DRAM
    if ((mHWConfig.L2WritePolicy == EWritePolicy::WRITE_BACK && mHWConfig.L2WriteAllocatePolicy == EWriteAllocatePolicy::LAZY_FETCH_ON_READ)
            || (mHWConfig.L2WritePolicy == EWritePolicy::WRITE_BACK && mHWConfig.L2WriteAllocatePolicy == EWriteAllocatePolicy::FETCH_ON_WRITE))
    {
        // For write-back L2 cache, write misses don't go to DRAM, so consider only reads. 
        // Even with FETCH_ON_WRITE, the model is similar since write acknowledgment is sent once the request reaches L2 cache.
        // The difference is immediate L2->DRAM BW overhead and fewer L2 read misses afterward. Not thoroughly analyzed.
        missDram = rMiss * l2RMissRate;
    }
    else if (mHWConfig.L2WritePolicy == EWritePolicy::LOCAL_WB_GLOBAL_WT && mHWConfig.L2WriteAllocatePolicy == EWriteAllocatePolicy::NO_WRITE_ALLOCATE)
        missDram = rMiss * l2RMissRate + wMiss * l2WMissRate;
    else
    {
        cout << "[Error] Unsupported L2 cache write policy" << endl;
        return -1;
    }
    
    double noContentionLatency = mHWConfig.L2Latency + l2RMissRate * mHWConfig.dramLatency; // TODO: L1 hit latency should be included
    
    double mshrBatchLatency, memDataQueueNoC, memDataQueueDRAM;
    if (isMD && 
            (nocSerLat * (rMiss + wMiss) * numActiveSM > (double) mHWConfig.L2Latency + mHWConfig.dramLatency))
    {
        // memory divergence and NoC saturated interval
        memDataQueueNoC = numActiveSM * (rMiss + wMiss) * nocSerLat;
        memDataQueueDRAM = numActiveSM * missDram * dramSerLat;

        // for each MSHR batching no contention read latency + NoC + DRAM queing delay occurs, 
        // but write queuing delay occur only once.
        mshrBatchLatency = noContentionLatency + numActiveSM * rMiss * nocSerLat + 
                numActiveSM * rMiss * l2RMissRate * dramSerLat;
    }
    else
    {
        memDataQueueNoC = 0.5 * numActiveSM * (rMiss + wMiss) * nocSerLat;
        memDataQueueDRAM = 0.5 * numActiveSM * missDram * dramSerLat;

        // for each MSHR batching no contention read latency + NoC + DRAM queing delay occurs, 
        // but write queuing delay occur only once.
        mshrBatchLatency = noContentionLatency + 0.5 * (numActiveSM * rMiss * nocSerLat + 
                numActiveSM * rMiss * l2RMissRate * dramSerLat);
    }

    double memStructMSHR = 0;
    if (isMD)
    {
        memStructMSHR = mshrBatchLatency *
            (ceil((double) intervalStat.l1Rmiss * smConcurrentWarps / mHWConfig.nMSHR) - 1);
    }

    intvCycleBreakdown.memDataQueueNoC += memDataQueueNoC;
    intvCycleBreakdown.memDataQueueDRAM += memDataQueueDRAM;
    intvCycleBreakdown.memDataQueue += memDataQueueNoC + memDataQueueDRAM;
    intvCycleBreakdown.memData += memDataQueueNoC + memDataQueueDRAM;
    intvCycleBreakdown.memStructMSHR += memStructMSHR;
    intvCycleBreakdown.memStruct += memStructMSHR;
    intvCycleBreakdown.tot += memDataQueueNoC + memDataQueueDRAM + memStructMSHR;

    return 0;
}