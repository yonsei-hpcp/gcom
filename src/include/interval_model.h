#ifndef _INTERVAL_MODEL_H
#define _INTERVAL_MODEL_H

#include <filesystem>
#include <vector>
#include <map>
#include "common.h"
#include "worker.h"
#include "inst_decoder.h"

namespace fs = std::filesystem;

namespace GCoM
{
    enum class EFUType
    {
        SP,
        DP,
        INT,
        SFU,
        MEM,
        SPECIALIZED_UNIT_1, // BRA
        SPECIALIZED_UNIT_2, // TEX
        SPECIALIZED_UNIT_3, // TENSOR
        SPECIALIZED_UNIT_4, // UDP
        SPECIALIZED_UNIT_5,
        SPECIALIZED_UNIT_6,
        SPECIALIZED_UNIT_7,
        SPECIALIZED_UNIT_8
    };

    static const std::map<EUArchOp, EFUType> TuringOpFUMap
    {
        {EUArchOp::SP_OP, EFUType::SP},
        {EUArchOp::DP_OP, EFUType::DP},
        {EUArchOp::ALU_OP, EFUType::INT},
        {EUArchOp::INTP_OP, EFUType::INT},
        {EUArchOp::SFU_OP, EFUType::SFU},
        {EUArchOp::LOAD_OP, EFUType::MEM},
        {EUArchOp::STORE_OP, EFUType::MEM},
        {EUArchOp::BARRIER_OP, EFUType::INT},
        {EUArchOp::MEMORY_BARRIER_OP, EFUType::MEM},
        {EUArchOp::EXIT_OPS, EFUType::INT},
        {EUArchOp::SPECIALIZED_UNIT_1_OP, EFUType::SPECIALIZED_UNIT_1},
        {EUArchOp::SPECIALIZED_UNIT_2_OP, EFUType::SPECIALIZED_UNIT_2},
        {EUArchOp::SPECIALIZED_UNIT_3_OP, EFUType::SPECIALIZED_UNIT_3},
        {EUArchOp::SPECIALIZED_UNIT_4_OP, EFUType::SPECIALIZED_UNIT_4},
        {EUArchOp::SPECIALIZED_UNIT_5_OP, EFUType::SPECIALIZED_UNIT_5},
        {EUArchOp::SPECIALIZED_UNIT_6_OP, EFUType::SPECIALIZED_UNIT_6},
        {EUArchOp::SPECIALIZED_UNIT_7_OP, EFUType::SPECIALIZED_UNIT_7},
        {EUArchOp::SPECIALIZED_UNIT_8_OP, EFUType::SPECIALIZED_UNIT_8}
    };

    class IntervalProfile
    {
    public:
        IntervalProfile() {};

        struct IntervalStat
        {
            unsigned nWarpInst = 0;
            unsigned stallCycle = 0;
            unsigned l1Rmiss = 0; // L1 read miss of rep warp
            unsigned l1Wmiss = 0; // L1 write miss of rep warp
            unsigned l2Rhit = 0; // L2 read hit of all warps
            unsigned l2Rmiss = 0; // L2 read miss of all warps
            unsigned l2Whit = 0; // L2 write hit of all warps
            unsigned l2Wmiss = 0; // L2 write miss of all warps
            unsigned nActiveThreadInst = 0;
            std::map<EFUType,unsigned> dispatch; // cycle to dispatch to each FU including cycle to emptying memaccess queue of a instruction
            bool hasMemInst = false; // includes memory access instruction on which another instruction has dependency
        };
        IntervalStat &operator[](int index)
        {
            return _interval_profile[index];
        }

        void push_back(IntervalStat intervalStat) { _interval_profile.push_back(intervalStat); }

        unsigned len() { return _interval_profile.size(); }
        
        double mIssueProb = 0; // issue probability of warp instruction only considering a single warp and instruction latency. (The same with GPUMech)
        double mAvgIntvWarpInsts = 0; // average number of warp instructions in an interval
    private:
        
        std::vector<IntervalStat> _interval_profile;
    };

    class GCoMModel
    {
    public:
        GCoMModel(HWConfig hwConfig)
        {
            mHWConfig = hwConfig;
        };

        // Data structure shared with external warp selector 
        struct KernelRepWarp
        { 
            int kernelNumber;
            int repWarpIdx;

            template<class Archive>
            void serialize(Archive &ar, const unsigned int version)
            {
                ar & kernelNumber;
                ar & repWarpIdx;
            }
        };
        // Export decoder.mWarps and read representative warp from the external warp selector
        // input: decoded warps
        // outpt: representative warp ID
        int SelectRepresentativeWarp(SASSDecoder *decoder, fs::path warpSelInputPath, fs::path repWarpPath, unsigned kernelIdx, int &representativeWarpID);

        struct DepTableEntry
        {
            unsigned doneCycle;
            WarpInst *warpInst;
        };
        // input: representativeWarp (wirh cache statistics), mHWConfig
        // output: mIntervalProfile (updated representativeWarp is not used anymore)
        int ProfileInterval(Warp &representativeWarp);

        struct Result
        {
            unsigned long long numThreadInst = 0;
            double numCycle = 0;
            double cpi = 0;
            double cpiBase = 0;
            double cpiComData = 0;
            double cpiComStruct = 0;
            double cpiMemData = 0;
            double cpiMemStruct = 0;
            double cpiIdle = 0;

            Result operator+(const Result &other) const;
        };
        // input: mIntervalProfile, mHWConfig, kernelInfo
        // output: result
        int RunPerformanceModel(KernelInfo kernelInfo, Result &result);
        
        struct CycleBreakdown
        {
            double tot = 0;
            double base = 0;
            double comData = 0;
            double comStruct = 0;
            double memData = 0;
            double memStruct = 0;
            double idle = 0;

            // sub-classification
            double memDataScbLast = 0;
            double memDataQueue = 0;
            double memDataQueueNoC = 0;
            double memDataQueueDRAM = 0;
            double memStructMSHR = 0;
            double memStructL1Bank = 0;

            CycleBreakdown operator+(const GCoMModel::CycleBreakdown &other) const;
            CycleBreakdown operator/(const unsigned long long &divisor) const;
            CycleBreakdown operator*(const unsigned &multiplier) const;
        };
    private:
        HWConfig mHWConfig;
        IntervalProfile mIntervalProfile;

        // Check source operand dependency of warpInst and decide issue cycle
        // input: warpInst, instIdx, dependency_table, prevIssueCycle
        // output: curIssueCycle, hasMemDep
        int CheckSrcDep(WarpInst &warpInst, std::map<unsigned, DepTableEntry> &dependency_table, unsigned prevIssueCycle,
                unsigned &curIssueCycle, bool &hasMemDep);

        // input: warpInst, instIdx, mHWConfig
        // output: warpInst.mDecoded.latency, intervalStat
        int DecideInstLatencyCollectIntervalStat(WarpInst &warpInst, unsigned instIdx, IntervalProfile::IntervalStat &intervalStat);

        // input: warpInst.mOpLatencyInitIntvType, mHWConfig
        // output: warpInst.mDecoded.latency, warpInst.mDecoded.initiationInterval
        int ReadLatencyInitIntv(WarpInst &warpInst);

        // Register file read latency (operand collector) model
        // input: warpInst.mDecoded.in.size(), mHWConfig.regFileReadThroughput
        // output: latency
        void ModelRegReadLatency(WarpInst &warpInst, unsigned &latency);

        // input: kernelInfo, mHWConfig
        // output: maxConcurrentCTAPerSM
        int calculateMaxConcurrentCTAPerSM(KernelInfo kernelInfo, unsigned &maxConcurrentCTAPerSM);

        // input: mIntervalProfile, allocedCTAPerSM, maxConcurrentCTAPerSM
        // output: SM cycle breakdown
        int RunPerSMModel(unsigned allocedCTAPerSM, unsigned maxConcurrentCTAPerSM, unsigned numWarpsPerCTA, unsigned numActiveSM,
                CycleBreakdown &smCycleBreakdown);

        // input: intervalIdx, smConcurrentWarps, numActiveSM, mIntervalProfile, mHWConfig
        // output: updated cycleBreakdown
        int RunPerIntervalModel(unsigned intervalIdx, unsigned smConcurrentWarps, unsigned numActiveSM, CycleBreakdown &smCycleBreakdown);

        // input: 
        //  intervalStat
        //  isLast: whether the interval is the last one
        //  numSubcore
        //  issueRate: issue rate of a subcore
        //  smConcurrentWarps, warpSchedulePolicy, mIntervalProfile.mIssueProb, mIntervalProfile.mAvgIntvWarpInsts
        // output: updated intvCycleBreakdown.{base, comData, memData, memDataScbLast, idle(intra-SM)}
        int RunMultiSubCoreModel(IntervalProfile::IntervalStat intervalStat, bool isLast, unsigned numSubcore, 
                double issueRate, unsigned smConcurrentWarps, std::string warpSchedulePolicy,
                CycleBreakdown &intvCycleBreakdown);

        // input:
        //  intervalStat, mIntervalProfile.mIssueProb, mIntervalProfile.mAvgIntvWarpInsts, warpSchedulePolicy
        //  issueRate: issue rate of a subcore
        //  numOthW: number of other warps in a subcore
        // output: 
        //  othBase: base cycle of other warps in a subcore
        //  dataStallSub: data stall cycle of a subcore
        //  issuedInstnInStall: number of issued instructions hiding stall cycles
        int RunSubcoreDataStallModel(IntervalProfile::IntervalStat intervalStat, std::string warpSchedulePolicy,
                double issueRate, double numOthW, double &othBase, double &dataStallSub, double &issuedInstnInStall);

        // input:
        //  intervalStat, numSubcore, smConcurrentWarps
        // output:
        //  updated intvCycleBreakdown.{comStruct, memStruct, memStructL1Bank}
        int RunIntraSMContentionModel(IntervalProfile::IntervalStat intervalStat, unsigned numSubcore, 
                unsigned smConcurrentWarps, unsigned operandCollectorQue, unsigned fuQue, CycleBreakdown &intvCycleBreakdown);

        // input:
        //  intervalStat, smConcurrentWarps, numActiveSM
        //  mHWConfig.{nMSHR, coreFreq, maxNoCBW, maxDRAMBW, l1DLineSize, nSectorPerLine, L2Latency, L2WritePolicy, L2WriteAllocatePolicy, dramLatency}
        // output:
        //  updated intvCycleBreakdown.{memData, memStruct, memDataQueue, memDataQueueNoC, memDataQueueDRAM, memStructMSHR}
        int RunMemoryContentionModel(IntervalProfile::IntervalStat intervalStat, unsigned smConcurrentWarps, unsigned numActiveSM,
                CycleBreakdown &intvCycleBreakdown);
    };
} // namespace GCoM

std::ostream& operator<<(std::ostream &os, const GCoM::GCoMModel::Result &result);

#endif
