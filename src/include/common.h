#ifndef _COMMON_H
#define _COMMON_H

#include <vector>
#include <string>

namespace GCoM
{
    typedef unsigned long long Address; // GPU virtual address type

    enum class EWritePolicy {
        WRITE_BACK,
        WRITE_THROUGH,
        LOCAL_WB_GLOBAL_WT
    };

    enum class EWriteAllocatePolicy {
        NO_WRITE_ALLOCATE,
        FETCH_ON_WRITE,
        LAZY_FETCH_ON_READ
    };

    enum class EHashFunction {
        LINEAR_SET_FUNCTION = 0,
        BITWISE_XORING_FUNCTION,
        HASH_IPOLY_FUNCTION,
        FERMI_HASH_SET_FUNCTION
    };

    // Defualt value is based in RTX2060
    struct HWConfig
    {
        unsigned numSMs = 30;
        unsigned numSubcorePerSM = 4;
        unsigned subcoreIssueRate = 1;
        unsigned maxThreadsPerSM = 1024;
        unsigned maxCTAPerSM = 16;
        unsigned registerPerSM = 65536;
        double coreFreq = 1.365; // GHz

        std::string wrapSchedulePolicy = "gto"; // [gto: greedy then oldest | lrr: loosely round-robin]

        // operation latency
        unsigned intLatency = 2;
        unsigned fpLatency = 2;
        unsigned dpLatency = 64;
        unsigned sfuLatency = 21;
        unsigned specializedUnit1Latency = 4;
        unsigned specializedUnit2Latency = 200;
        unsigned specializedUnit3Latency = 16;
        unsigned specializedUnit4Latency = 4;
        unsigned specializedUnit5Latency;
        unsigned specializedUnit6Latency;
        unsigned specializedUnit7Latency;
        unsigned specializedUnit8Latency;

        // operation initiation interval
        unsigned intInitIntv = 2;
        unsigned fpInitIntv = 2;
        unsigned dpInitIntv = 64;
        unsigned sfuInitIntv = 8;
        unsigned specializedUnit1IntiIntv = 4;
        unsigned specializedUnit2IntiIntv = 4;
        unsigned specializedUnit3IntiIntv = 16;
        unsigned specializedUnit4IntiIntv = 1;
        unsigned specializedUnit5IntiIntv;
        unsigned specializedUnit6IntiIntv;
        unsigned specializedUnit7IntiIntv;
        unsigned specializedUnit8IntiIntv;

        // static latency of pipeline stages. register file read, Bank, FU lane will be modeled separtely
        unsigned computePipeline = 3;
        unsigned loadPipeline = 0;
        unsigned storePipeline = 4;

        unsigned operandCollectorQue = 8; // # of operand collector queue per SM
        unsigned fuQue = 2; // # of queue per functional unitS

        unsigned regFileReadThroughput = 2;

        // OpFUMap // Opcode to FU type mapping is hardcoded in interval_model.h targeting Turing for now

        // L1 cache (per SM)
        unsigned l1DLatency = 32;
        EWritePolicy l1DWritePolicy = EWritePolicy::WRITE_THROUGH;
        EWriteAllocatePolicy l1DWriteAllocatePolicy = EWriteAllocatePolicy::LAZY_FETCH_ON_READ;
        
        unsigned l1DnSet = 4; // Total number of sets (in all banks)
        unsigned l1DAssoc = 192; // Max associativity of L1D cache
        unsigned l1DLineSize = 128; // A line can be splited into multiple sectors
        EHashFunction l1DSetHashFunction = EHashFunction::LINEAR_SET_FUNCTION;
        unsigned l1DnBank = 4; // Number of L1D banks where each bank has part of sets or sectors. Lines (or sectors) in a set cannot be accessed in parallel.
        EHashFunction l1DBankHashFunction = EHashFunction::LINEAR_SET_FUNCTION;
        unsigned l1DBankInterleaveByte = 32; // Could be smaller than line size in sectored cache
        unsigned nSectorPerLine = 4;
        unsigned nMSHR = 256; // Number of MSHR entries
        
        unsigned sharedMemLatency = 30; // L1S
        std::vector<unsigned> sharedMemSizeOption = {32768, 65536}; // Should be ascending order

        // L2 cache
        unsigned L2Latency = 194;
        EWritePolicy L2WritePolicy = EWritePolicy::WRITE_BACK;
        EWriteAllocatePolicy L2WriteAllocatePolicy = EWriteAllocatePolicy::LAZY_FETCH_ON_READ;

        unsigned l2DnSet = 64; // Number of sets in a L2D partition
        unsigned l2DAssoc = 16; // Number of ways in a set
        unsigned l2DLineSize = 128;
        EHashFunction l2DSetHashFunction = EHashFunction::HASH_IPOLY_FUNCTION;
        // unsigned nSectorPerLine = 4; // use the same value as L1D
        unsigned l2DnPerMem = 2; // number of L2D partitions for a memory channel

        // interconnect
        double maxNoCBW = 1048.32; // max NoC bandwidth (GB/s)

        // memory
        unsigned dramLatency = 96;
        double maxDRAMBW = 336; // max DRAM bandwidth (GB/s)
        unsigned numMems = 12; // number of memory channels
    };

    struct KernelInfo 
    {
        unsigned sharedMemSize = 0;
        unsigned numRegs = 0;
        unsigned numCTAs = 0; // total
        unsigned numWarps = 0; // total
        unsigned long long numThreadInsts = 0;
    };

    extern unsigned DEBUG_LEVEL;

    // copied from Accel-sim gpu-misc.cc LOGB2
    unsigned int LogB2(unsigned int v);
}

#endif