#ifndef _INSTRUCTION_H
#define _INSTRUCTION_H

#include <vector>
#include <list>
#include <cstdint>
#include <bitset>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/list.hpp>
#include <boost/serialization/bitset.hpp>
#include "common.h"
#include "trace_parser.h"

#define MAX_WARP_SIZE 32
#define MAX_SECTOR_CHUNCK_SIZE 4

namespace GCoM
{
    typedef std::bitset<MAX_WARP_SIZE> ActiveMask;
    typedef std::bitset<MAX_SECTOR_CHUNCK_SIZE> MemAccessSectorMask;

    enum class EUArchOp
    {
        NO_OP = -1,
        ALU_OP = 1,
        SFU_OP,
        TENSOR_CORE_OP,
        DP_OP,
        SP_OP,
        INTP_OP,
        ALU_SFU_OP,
        LOAD_OP,
        TENSOR_CORE_LOAD_OP, // not used
        TENSOR_CORE_STORE_OP, // not used
        STORE_OP,
        BRANCH_OP,
        BARRIER_OP,
        MEMORY_BARRIER_OP,
        CALL_OPS,
        RET_OPS,
        EXIT_OPS,
        SPECIALIZED_UNIT_1_OP,
        SPECIALIZED_UNIT_2_OP,
        SPECIALIZED_UNIT_3_OP,
        SPECIALIZED_UNIT_4_OP,
        SPECIALIZED_UNIT_5_OP,
        SPECIALIZED_UNIT_6_OP,
        SPECIALIZED_UNIT_7_OP,
        SPECIALIZED_UNIT_8_OP
    };

    enum class EMemorySpace
    {
        UNDEFINED_SPACE = 0,
        GLOBAL_SPACE,
        LOCAL_SPACE,
        SHARED_SPACE,
        CONST_SPACE,
        TEX_SPACE,
        INSTRUCTION_SPACE
    };

    enum class ECacheOpType
    {
        CACHE_UNDEFINED,

        // loads
        CACHE_ALL,       // .ca
        CACHE_LAST_USE,  // .lu
        CACHE_VOLATILE,  // .cv
        CACHE_L1,        // .nc

        // loads and stores
        CACHE_STREAMING,  // .cs
        CACHE_GLOBAL,     // .cg

        // stores
        CACHE_WRITE_BACK,    // .wb
        CACHE_WRITE_THROUGH  // .wt
    };

    enum class EMemAccessType
    {
        GLOBAL_ACC_R,
        LOCAL_ACC_R,
        TEXTURE_ACC_R,
        GLOBAL_ACC_W,
        LOCAL_ACC_W,
        TEXTURE_ACC_w,
        L1_WRBK_ACC,
        L2_WRBK_ACC,
        INST_ACC_R,
        L1_WR_ALLOC_R,
        L2_WR_ALLOC_R
    };

    enum class EOpLatencyInitIntvType
    {
        DEFALUT,
        INT,
        FP,
        FP16,
        DP,
        SFU,
        TENSOR,
        SPECIALIZED_UNIT_1,
        SPECIALIZED_UNIT_2,
        SPECIALIZED_UNIT_3,
        SPECIALIZED_UNIT_4,
        SPECIALIZED_UNIT_5,
        SPECIALIZED_UNIT_6,
        SPECIALIZED_UNIT_7,
        SPECIALIZED_UNIT_8
    };

    class MemAccess
    {
    public:
        MemAccess() {}
        MemAccess(Address addr, bool isWrite, unsigned reqSize, EMemAccessType type, MemAccessSectorMask sectorMask)
            : mAddr(addr), mIsWrite(isWrite), mReqSize(reqSize), mType(type), mSectorMask(sectorMask) {}

        Address mAddr;
        bool mIsWrite;
        unsigned mReqSize; // bytes
        EMemAccessType mType;
        MemAccessSectorMask mSectorMask;

    private:
        friend class boost::serialization::access;
        template<class Archive>
        void serialize(Archive &ar, const unsigned int version)
        {
            ar & mAddr;
            ar & mIsWrite;
            ar & mReqSize;
            ar & mType;
            ar & mSectorMask;
        }
    };

    class Warp; // Forward declaration

    class WarpInst
    {
    public:
        WarpInst() 
        {
            //mWarpPtr = warpPtr;
        }
        ~WarpInst() {}

        //unsigned mInstsId = 3;

        unsigned mDataSize; // (B) need for each thread
        std::vector<Address> mPerThreadMemAddr; // memory address for each thread
        EOpLatencyInitIntvType mOpLatencyInitIntvType; // operation latency and initiation interval type
        struct 
        {
            Address pc = (Address) -1;
            EUArchOp op = EUArchOp::NO_OP;
            std::vector<unsigned> out; // dst register indexs
            std::vector<unsigned> in; // src register indexs
            int pred; // predicate register index
            unsigned latency = 0; // operation latency
            unsigned initiationInterval = 0;
            EMemorySpace space = EMemorySpace::UNDEFINED_SPACE;
            ECacheOpType cacheOp = ECacheOpType::CACHE_UNDEFINED;
            bool isAtomic = false;
            ActiveMask activeMask; // dynamic active mask for timing model (after predication)
            std::list<MemAccess> accessQ; // memory accesses after coalescing

            template<class Archive>
            void serialize(Archive &ar, const unsigned int version)
            {
                ar & pc;
                ar & op;
                ar & out;
                ar & in;
                ar & pred;
                ar & latency;
                ar & initiationInterval;
                ar & space;
                ar & cacheOp;
                ar & isAtomic;
                ar & activeMask;
                ar & accessQ;
            }
        } mDecoded;

        struct MemStat
        {
            uint8_t l1Hit = 0; // # of thread memory request that hit L1
            uint8_t l1Miss = 0; // # of all thread memory request that missed L1
            uint8_t coalescedL1Miss = 0; // L1 <-> L2 mem request #. Miss on the same MSHR entry will be coalesced.
            uint8_t l2Hit = 0; // # of L1 <-> L2 mem request that hit L2
            uint8_t l2Miss = 0; // # of L1 <-> L2 mem request that missed L2

            MemStat& operator=(const struct WarpInstCacheStat& other);

            template<class Archive>
            void serialize(Archive &ar, const unsigned int version)
            {
                ar & l1Hit;
                ar & l1Miss;
                ar & coalescedL1Miss;
                ar & l2Hit;
                ar & l2Miss;
            }
        } mMemStat;

        //inst_trace_t *mInstTracePtr; // trace before decoding // not necessary?
        Warp *mWarpPtr; // backward pointer to warp class
    
    private:
        friend class boost::serialization::access;
        template<class Archive>
        void serialize(Archive &ar, const unsigned int version)
        {
            ar & mDataSize;
            ar & mPerThreadMemAddr;
            ar & mOpLatencyInitIntvType;
            ar & mDecoded;
            ar & mMemStat;
            ar & mWarpPtr;
        }
    };

    class Warp
    {
    public:
        Warp() {}
        ~Warp() {}

        // custom copy constructor
        Warp(const Warp &other)
        {
            mInsts = other.mInsts;
            for (WarpInst &inst : mInsts)
                inst.mWarpPtr = this;
            mGlobalCacheStat = other.mGlobalCacheStat;
        }

        std::vector<WarpInst> mInsts;
        // struct
        // {
        //     std::string kernelName;
        //     unsigned ctaId; // CTA id among all CTAs in a kernel // not neccessary? 
        //     unsigned warpId; // warp id among warps in a CTA // not neccessary?
        // } mInfo;

        struct WarpInstCacheHitMissCount
        {
            unsigned l1HitWarp = 0; // The number of warp instructions that did not miss L1 cache at all
            unsigned l1MissWarp = 0; //The number of warp instructions that missed L1
            unsigned l2HitWarp = 0; // The number of warp instructions that did not miss L2 cache at all
            unsigned l2MissWarp = 0; //The number of warp instructions that missed L2

            template<class Archive>
            void serialize(Archive &ar, const unsigned int version)
            {
                ar & l1HitWarp;
                ar & l1MissWarp;
                ar & l2HitWarp;
                ar & l2MissWarp;
            }
        };
        std::vector<WarpInstCacheHitMissCount> mGlobalCacheStat;
    
    private:
        friend class boost::serialization::access;
        template<class Archive>
        void serialize(Archive &ar, const unsigned int version)
        {
            ar & mInsts;
            ar & mGlobalCacheStat;
        }
    };
}

#endif