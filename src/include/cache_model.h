#ifndef _CACHE_MODEL_H
#define _CACHE_MODEL_H

#include <map>
#include <vector>
#include <filesystem>
#include "instruction.h"
#include <string>
#include <vector>
#include <bitset>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/vector.hpp>
#include "inst_decoder.h"
#include "addrdec.h"

namespace fs = std::filesystem;

namespace GCoM
{
    enum class ECacheModelType
    {
        SIMPLE // Simple wrte-back cache
    };

    class CacheBlock
    {
        public:
        CacheBlock() {
            for (unsigned i = 0; i < MAX_SECTOR_CHUNCK_SIZE; i++)
                m_status.push_back(ECacheBlockState::INVALID);
        }

        enum class ECacheBlockState { INVALID = 0, RESERVED, VALID, MODIFIED };
        ECacheBlockState GetStatus(unsigned sidx);

        void Fill(unsigned sidx);
        
        void Clear();

        private:
        std::vector<ECacheBlockState> m_status;
    };

    class CacheSet
    {
        public:
        CacheSet(unsigned int associativity)
        {
            assoc = associativity;
            tags.resize(assoc);
            vals.resize(assoc);
            lines.resize(assoc);
            order.resize(assoc);
            for (unsigned i = 0; i < assoc; i++)
            {
                vals[i] = false;
                order[i] = assoc - 1;
            }
        }

        std::vector<Address> tags;
        std::vector<bool> vals;
        std::vector<CacheBlock> lines;

        std::vector<int> order;
        
        private:
        unsigned int assoc;
    };

    class SimpleCache
    {
    public:
        SimpleCache() {}
        SimpleCache(unsigned int associativity, unsigned int numSet)
        {
            assoc = associativity;
            nset = numSet;
            
            for (unsigned i = 0; i < nset; i ++)
                sets.push_back(CacheSet(assoc));
        }

        enum class ECacheRequestStatus {
            HIT = 0,
            HIT_RESERVED,
            MISS,
            RESERVATION_FAIL,
            SECTOR_MISS,
            MSHR_HIT,
            HIT_RESERVED_JH
        };
        ECacheRequestStatus access(unsigned set_index, Address tag, MemAccessSectorMask mask);

    private:
        unsigned int assoc;
        unsigned int nset;

        std::vector<CacheSet> sets;

        unsigned mask2idx(MemAccessSectorMask sector_mask);
    };

    // warper class for cache model
    class CacheModel
    {
    public:
        CacheModel() {}
        CacheModel(ECacheModelType type)
        {
            mType = type;
        }

        // Initialization for Simple cache
        CacheModel(ECacheModelType type, HWConfig hwConfig);

        // based on shader_core_config::max_cta of Accel-sim
        // input: kernelInfo, hwConfig
        // output: updated L1D cache size, mMaxCtaPerSm
        int SetSharedMemSizeMaxCta(KernelInfo kernelInfo, HWConfig hwConfig);

        /*
        * Access cache and set mMemStat of WarpInst
        * return 0 for normal exit. -1 for abnormal exit.
        */
        int AccessCache(WarpInst &inst, unsigned smId, HWConfig hwConfig, unsigned instIdx);

        // update global cache statistics of a warp with collected warp instruction statistics
        int UpdateGlobalCacheStat(Warp &warp);
        
        unsigned mMaxCtaPerSm;

        static const unsigned MAX_MEMORY_ACCESS_SIZE = 128; // From Accel-sim abstract_hardware_model.h
        typedef std::bitset<MAX_MEMORY_ACCESS_SIZE> MemAccessByteMask;
        struct transactionInfo {
            std::bitset<4> chunks;  // bitmask: 32-byte chunks accessed
            MemAccessByteMask bytes;
            ActiveMask active;  // threads in this transaction

            bool test_bytes(unsigned start_bit, unsigned end_bit) {
                for (unsigned i = start_bit; i <= end_bit; i++)
                    if (bytes.test(i)) 
                        return true;
                return false;
            }
        };
    private:
        ECacheModelType mType;

        // warp instruction statistic collected during cahce access
        std::map<std::pair<Address, unsigned>, Warp::WarpInstCacheHitMissCount>
                mWarpInstStatistics_PCInstIdx2CacheHitMiss;

        // for simple cache model only
        std::vector<SimpleCache> mL1D;
        std::vector<SimpleCache> mL2D;
        linear_to_raw_address_translation mAddressMapping;

        // do memory coalescing and generate MemAccess
        int GenerateMemAccesses(WarpInst &inst, HWConfig hwConfig);
        int CoalesceMemoryAccess(WarpInst &inst, EMemAccessType accessType, HWConfig hwConfig);
        int CoalesceAtomicMemoryAccess(WarpInst &inst, EMemAccessType accessType, HWConfig hwConfig);
        void ReduceTransactionSize( bool is_write, EMemAccessType access_type, const transactionInfo &info, Address addr, unsigned segmentSize, std::list<MemAccess> &accessQ);
    };

    // shared by CacheModel and GCoMModel
    // input: kernelInfo, hwConfig
    // output: maxConcurrentCTAPerSM
    int CalculateMaxConcurrentCTAPerSM(KernelInfo kernelInfo, HWConfig hwConfig, unsigned &maxConcurrentCTAPerSM);

    // Copide from Accel-sim cache_config::hash_function
    // Currently bank index hashing is used in a warp profiling. Otherwise this may go inside CacheModel.
    // L1, L2 set, bank index hashing function
    unsigned HashAddress(Address addr, unsigned nBins, unsigned offSetBits, EHashFunction hashFunctionType);

    // Emulate TB allocation and warp scheduling, access cache
    int RunSimpleCache(CacheModel *simpleCache, SASSDecoder *decoder, HWConfig hwConfig);
} // namespace GCoM

#endif