#include "cache_model.h"
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <filesystem>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/vector.hpp>

using namespace GCoM;
using namespace std;
namespace fs = std::filesystem;

#define WARP_SIZE 32 // Strictly speaking it is HW specification. But let's consider it HW independent for now since it is implicitly included in SASS trace.

CacheBlock::ECacheBlockState CacheBlock::GetStatus(unsigned sidx)
{
    return m_status[sidx];
}

void CacheBlock::Fill(unsigned sidx)
{
    m_status[sidx] = ECacheBlockState::VALID;
}

void CacheBlock::Clear()
{
    for (unsigned i = 0; i < MAX_SECTOR_CHUNCK_SIZE; ++i) {
        m_status[i] = ECacheBlockState::INVALID;
    }
}

SimpleCache::ECacheRequestStatus SimpleCache::access(unsigned setIndex, Address tag, MemAccessSectorMask mask)
{
    CacheSet &set = sets[setIndex];
    unsigned sidx = mask2idx(mask);
    int emptyIdx = -1;
    int evicCandi = -1;

    for (unsigned i = 0; i < assoc; i++)
    {
        if (set.vals[i])
        {
            if(set.tags[i] == tag)
            {
                /* access order for LRU policy : 0 if hit, ++ oder of others */
                int old = set.order[i];
                for (unsigned j = 0; j < assoc; j++)
                {
                    if (set.order[j] < old)
                        set.order[j] ++;
                }
                set.order[i] = 0;
                if (set.lines[i].GetStatus(sidx) == CacheBlock::ECacheBlockState::VALID)
                    return SimpleCache::ECacheRequestStatus::HIT;
                else //sector miss
                {
                    set.lines[i].Fill(sidx);
                    return SimpleCache::ECacheRequestStatus::MISS;
                }
            }
            else if (set.order[i] == (int) assoc - 1)
            {
                evicCandi = i; // evict block candiadate
            }

            //how to deal with global mem cache coherence?
        }
        else
        {
            emptyIdx = i; // keep empty block.
        }
    }

    /* when miss */
    int fillingIdx;
    if (emptyIdx != -1)
        fillingIdx = emptyIdx;
    else //evict
    {
        fillingIdx = evicCandi;
        set.lines[fillingIdx].Clear();
    }
    set.tags[fillingIdx] = tag;
    set.vals[fillingIdx] = true;
    set.lines[fillingIdx].Fill(sidx);
    int old = set.order[fillingIdx];
    for (unsigned j = 0; j < assoc; j++)
    {
        if (set.order[j] < old)
            set.order[j]++;
    }
    set.order[fillingIdx] = 0;
    return SimpleCache::ECacheRequestStatus::MISS;
    
}

unsigned SimpleCache::mask2idx(MemAccessSectorMask sectorMask)
{
    assert(sectorMask.count() == 1);
    for (unsigned i = 0; i < MAX_SECTOR_CHUNCK_SIZE; ++i) { // to further optimize, this code could go to SimpleCache::access and change sectorMask to sector index
        if (sectorMask.to_ulong() & (1 << i)) 
            return i;
    }
    return -1; // cannot reach here
}

// Initialization for Simple cache mode
CacheModel::CacheModel(ECacheModelType type, HWConfig hwConfig)
{
    mType = type;
    assert(type == ECacheModelType::SIMPLE);

    // L1D cache will be initialized in SetSharedMemSizeMaxCta
    
    for (unsigned i = 0; i < hwConfig.numMems * hwConfig.l2DnPerMem; i++)
        mL2D.push_back(SimpleCache(hwConfig.l2DAssoc, hwConfig.l2DnSet));

    mAddressMapping.addrdec_setoption();
    mAddressMapping.init(hwConfig.numMems, hwConfig.l2DnPerMem);
}

// based on shader_core_config::max_cta of Accel-sim
int CacheModel::SetSharedMemSizeMaxCta(KernelInfo kernelInfo, HWConfig hwConfig)
{
    if (CalculateMaxConcurrentCTAPerSM(kernelInfo, hwConfig, mMaxCtaPerSm) != 0)
        return -1;

    // Set L1D cache size according to shared memory requirement (adaptive cache)
    unsigned totSharedMem = kernelInfo.sharedMemSize * mMaxCtaPerSm;
    if (totSharedMem < 0 || totSharedMem > hwConfig.sharedMemSizeOption.back())
    {
        cout << "[Error] Unsupported shared memory requirement" << endl;
        return -1;
    }

    unsigned totalUnifiedL1 = hwConfig.l1DAssoc * hwConfig.l1DnSet * hwConfig.l1DLineSize;
    unsigned newL1DAssoc;
    for (unsigned sharedMemSizeOpt : hwConfig.sharedMemSizeOption)
    {
        if (totSharedMem <= sharedMemSizeOpt)
        {

            newL1DAssoc = hwConfig.l1DAssoc * (totalUnifiedL1 - sharedMemSizeOpt);
            if (newL1DAssoc % totalUnifiedL1 != 0 || totalUnifiedL1 <= sharedMemSizeOpt)
            {
                cout << "[Error] Invalid shared memory size option" << endl;
                return -1;
            }
            newL1DAssoc = newL1DAssoc / totalUnifiedL1;
            break;
        }
    }

    for (unsigned i = 0; i < hwConfig.numSMs; i++)
        mL1D.push_back(SimpleCache(newL1DAssoc, hwConfig.l1DnSet));

    return 0;
}

inline Address blockAddr(Address addr, unsigned lineSize)
{
    return addr & ~(Address)(lineSize - 1);
}

inline Address tagAddr(Address addr, unsigned lineSize)
{
    // For generality, the tag includes both index and tag. This allows for more
    // complex set index calculations that can result in different indexes
    // mapping to the same set, thus the full tag + index is required to check
    // for hit/miss. Tag is now identical to the block address.

    return addr & ~(Address)(lineSize - 1);
}

// for Simple cache mode
int CacheModel::AccessCache(WarpInst &inst, unsigned smId, HWConfig hwConfig, unsigned instIdx)
{
    assert(hwConfig.nSectorPerLine == MAX_SECTOR_CHUNCK_SIZE || hwConfig.nSectorPerLine == 1);
    if (GenerateMemAccesses(inst, hwConfig) != 0)
        return -1;

    for (auto access : inst.mDecoded.accessQ)
    {
        // if mem_space == local, global
        Address addr = access.mAddr;
        Address baddr = blockAddr(addr, hwConfig.l1DLineSize);
        unsigned set_index = HashAddress(baddr, hwConfig.l1DnSet, LogB2(hwConfig.l1DLineSize), hwConfig.l1DSetHashFunction);
        Address tag = tagAddr(baddr, hwConfig.l1DLineSize);
        MemAccessSectorMask mask = MemAccessSectorMask(1);
        if (hwConfig.nSectorPerLine == MAX_SECTOR_CHUNCK_SIZE)
            mask = access.mSectorMask;
        SimpleCache::ECacheRequestStatus l1Status = mL1D[smId].access(set_index, tag, mask);
        if (l1Status == SimpleCache::ECacheRequestStatus::MISS)
        {
            inst.mMemStat.l1Miss++;
            inst.mMemStat.coalescedL1Miss++;
            baddr = blockAddr(addr, hwConfig.l2DLineSize);
            Address part_addr = mAddressMapping.partition_address(baddr);
            set_index = HashAddress(part_addr, hwConfig.l2DnSet, LogB2(hwConfig.l2DLineSize), hwConfig.l2DSetHashFunction);
            tag = tagAddr(baddr, hwConfig.l2DLineSize);
            addrdec_t m_raw_addr;
            mAddressMapping.addrdec_tlx(addr, &m_raw_addr);
            SimpleCache::ECacheRequestStatus l2Status = mL2D[m_raw_addr.sub_partition].access(set_index, tag, mask);
            if (l2Status == SimpleCache::ECacheRequestStatus::MISS)
                inst.mMemStat.l2Miss++;
            else
                inst.mMemStat.l2Hit++;
        }
        else
            inst.mMemStat.l1Hit++;
    }

    // Collect warp instruction cache statistics
    auto key = make_pair(inst.mDecoded.pc, instIdx);
    // make element in mWarpInstStatistics_PCInstIdx2CacheHitMiss if not exist, else reference it
    Warp::WarpInstCacheHitMissCount &warpInstCStat = mWarpInstStatistics_PCInstIdx2CacheHitMiss[key];

    if (inst.mMemStat.l1Miss > 0)
        warpInstCStat.l1MissWarp += 1;
    else
        warpInstCStat.l1HitWarp += 1;

    if (inst.mMemStat.l2Miss > 0)
        warpInstCStat.l2MissWarp += 1;
    else
        warpInstCStat.l2HitWarp += 1;

    return 0;
}

int CacheModel::UpdateGlobalCacheStat(Warp &warp)
{
    for (unsigned instIdx = 0; instIdx < warp.mInsts.size(); instIdx++)
    {
        WarpInst &warpInst = warp.mInsts[instIdx];

        // find corresponding global cache statistics of a warp instruction
        // if not found, set to 0
        auto key = make_pair(warpInst.mDecoded.pc, instIdx);
        auto it = mWarpInstStatistics_PCInstIdx2CacheHitMiss.find(key);
        if (it != mWarpInstStatistics_PCInstIdx2CacheHitMiss.end())
            warp.mGlobalCacheStat.push_back(it->second);
        else
            warp.mGlobalCacheStat.push_back(Warp::WarpInstCacheHitMissCount());
    }
    return 0;
}

// Based on Accel-sim void warp_inst_t::generate_mem_accesses()
int CacheModel::GenerateMemAccesses(WarpInst &inst, HWConfig hwConfig)
{
    if (inst.mDecoded.activeMask.count() == 0)
        return 0; // predicated off
    if (inst.mPerThreadMemAddr.size() == 0)
    {
        cout << "[Error] load/store with no memory address" << endl;
        return -1; // no memory access
    }

    EMemAccessType accessType;
    if (inst.mDecoded.space == EMemorySpace::GLOBAL_SPACE)
        accessType = inst.mDecoded.op == EUArchOp::STORE_OP ? EMemAccessType::GLOBAL_ACC_W : EMemAccessType::GLOBAL_ACC_R;
    else if (inst.mDecoded.space == EMemorySpace::LOCAL_SPACE)
        accessType = inst.mDecoded.op == EUArchOp::STORE_OP ? EMemAccessType::LOCAL_ACC_W : EMemAccessType::LOCAL_ACC_R;
    else if (inst.mDecoded.space == EMemorySpace::SHARED_SPACE)
        return 0; // currently GCoM ignores bank conflict, multiple access for simplicity
    else
    {
        cout << "[Error] Invalid memory space" << endl;
        if (inst.mDecoded.space == EMemorySpace::CONST_SPACE)
            return -1; // not used
        else if (inst.mDecoded.space == EMemorySpace::TEX_SPACE)
            return -1; // Accel-sim SASS mode ignores texture loads, consider it as ALU_OP
        else
            return -1; // should not reach here
    }
    
    // Calculate memory accesses generated by this warp
    if (inst.mDecoded.space == EMemorySpace::GLOBAL_SPACE || inst.mDecoded.space == EMemorySpace::LOCAL_SPACE)
    {
        if (inst.mDecoded.isAtomic)
            return CoalesceAtomicMemoryAccess(inst, accessType, hwConfig);
        else
            return CoalesceMemoryAccess(inst, accessType, hwConfig);
    } else
        return -1; // should not reach here
    
}

// Based on Accel-sim Address line_size_based_tag_func(...)
inline Address GetLineSizeBasedTag(Address address,
                                       Address line_size) {
  // gives the tag for an address based on a given line size
  return address & ~(line_size - 1);
};

// Based on Accel-sim void warp_inst_t::memory_coalescing_arch(...)
int CacheModel::CoalesceMemoryAccess(WarpInst &inst, EMemAccessType accessType, HWConfig hwConfig)
{
    // see the CUDA manual where it discusses coalescing rules before reading this
    unsigned segmentSize = 0;
    if (hwConfig.nSectorPerLine != 1 && hwConfig.nSectorPerLine != 4)
    {
        cout << "[Error] Unsupported sector size" << endl;
        return -1;
    }
    if (hwConfig.l1DLineSize != 128)
    {
        cout << "[Error] Unsupported line size" << endl;
        return -1;
    }
    bool sectorSegmentSize = (hwConfig.nSectorPerLine != 1);

    switch (inst.mDataSize)
    {
    case 1:
        segmentSize = 32;
        break;
    case 2:
        segmentSize = sectorSegmentSize ? 32 : 64;
        break;
    case 4:
    case 8:
    case 16:
        segmentSize = sectorSegmentSize ? 32 : 128;
        break;
    }
    std::map<Address, transactionInfo> warpTransactions;

    // step 1: find all transactions generated by this warp
    for (unsigned thread = 0;
            thread < WARP_SIZE; thread++)
    {
        if (inst.mDecoded.activeMask.test(thread) == false)
            continue;

        unsigned dataSizeCoales = inst.mDataSize;

        if (inst.mDecoded.space == EMemorySpace::LOCAL_SPACE)
        {
            // Local memory accesses >4B were split into 4B chunks
            if (inst.mDataSize >= 4)
                dataSizeCoales = 4;
            // Otherwise keep the same inst.mDataSize for sub-4B access to local memory
        }

        Address addr = inst.mPerThreadMemAddr[thread];
        Address blockAddress =
            GetLineSizeBasedTag(addr, segmentSize);
        unsigned chunk =
            (addr & 127) / 32; // which 32-byte chunk within in a 128-byte
                                // chunk does this thread access?
        transactionInfo &info = warpTransactions[blockAddress];

        // can only write to one segment
        // it seems like in trace driven, a thread can write to more than one
        // segment assert(blockAddress ==
        // GetLineSizeBasedTag(addr+dataSizeCoales-1,segmentSize));

        info.chunks.set(chunk);
        info.active.set(thread);
        unsigned idx = (addr & 127);
        for (unsigned i = 0; i < dataSizeCoales; i++)
            if ((idx + i) < MAX_MEMORY_ACCESS_SIZE)
                info.bytes.set(idx + i);

        // it seems like in trace driven, a thread can write to more than one
        // segment handle this special case
        if (blockAddress != GetLineSizeBasedTag(
                                    addr + dataSizeCoales - 1, segmentSize))
        {
            addr = addr + dataSizeCoales - 1;
            Address blockAddress =
                GetLineSizeBasedTag(addr, segmentSize);
            unsigned chunk = (addr & 127) / 32;
            transactionInfo &info = warpTransactions[blockAddress];
            info.chunks.set(chunk);
            info.active.set(thread);
            unsigned idx = (addr & 127);
            for (unsigned i = 0; i < dataSizeCoales; i++)
                if ((idx + i) < MAX_MEMORY_ACCESS_SIZE)
                    info.bytes.set(idx + i);
        }
    
    }

    // step 2: reduce each transaction size, if possible
    std::map<Address, transactionInfo>::iterator t;
    for (t = warpTransactions.begin(); t != warpTransactions.end();
            t++)
    {
        Address addr = t->first;
        const transactionInfo &info = t->second;

        ReduceTransactionSize((inst.mDecoded.op == EUArchOp::STORE_OP), accessType, info, addr,
                                                segmentSize, inst.mDecoded.accessQ);
    }
    return 0;
}

// Based on Accel-sim void warp_inst_t::memory_coalescing_arch_atomic(...)
int CacheModel::CoalesceAtomicMemoryAccess(WarpInst &inst, EMemAccessType accessType, HWConfig hwConfig)
{
    if(inst.mDecoded.space != EMemorySpace::GLOBAL_SPACE)
    {
        cout << "[Error] Atomic operations allowed only for global memory" << endl;
        return -1; // Atomics allowed only for global memory
    }   

    // see the CUDA manual where it discusses coalescing rules before reading this
    unsigned segmentSize = 0;
    if (hwConfig.nSectorPerLine != 1 && hwConfig.nSectorPerLine != 4)
    {
        cout << "[Error] Unsupported sector size" << endl;
        return -1;
    }
    if (hwConfig.l1DLineSize != 128)
    {
        cout << "[Error] Unsupported line size" << endl;
        return -1;
    }
    bool sectorSegmentSize = (hwConfig.nSectorPerLine != 1);

    switch (inst.mDataSize)
    {
    case 1:
        segmentSize = 32;
        break;
    case 2:
        segmentSize = sectorSegmentSize ? 32 : 64;
        break;
    case 4:
    case 8:
    case 16:
        segmentSize = sectorSegmentSize ? 32 : 128;
        break;
    }

    std::map<Address, std::list<transactionInfo> > warpTransactions; // each block addr maps to a list of transactions

    // step 1: find all transactions generated by this warp
    for (unsigned thread = 0;
            thread < WARP_SIZE; thread++)
    {
        if (inst.mDecoded.activeMask.test(thread) == false)
            continue;

        Address addr = inst.mPerThreadMemAddr[thread];
        Address blockAddress =
            GetLineSizeBasedTag(addr, segmentSize);
        unsigned chunk =
            (addr & 127) / 32; // which 32-byte chunk within in a 128-byte chunk
                                // does this thread access?

        // can only write to one segment
        assert(blockAddress ==
                GetLineSizeBasedTag(addr + inst.mDataSize - 1, segmentSize));
        
        // Find a transaction that does not conflict with this thread's accesses
        bool new_transaction = true;
        std::list<transactionInfo>::iterator it;
        transactionInfo *info;
        for (it = warpTransactions[blockAddress].begin();
                it != warpTransactions[blockAddress].end(); it++)
        {
            unsigned idx = (addr & 127);
            if (not it->test_bytes(idx, idx + inst.mDataSize - 1))
            {
                new_transaction = false;
                info = &(*it);
                break;
            }
        }
        if (new_transaction)
        {
            // Need a new transaction
            warpTransactions[blockAddress].push_back(transactionInfo());
            info = &warpTransactions[blockAddress].back();
        }
        assert(info);

        info->chunks.set(chunk);
        info->active.set(thread);
        unsigned idx = (addr & 127);
        for (unsigned i = 0; i < inst.mDataSize; i++)
        {
            assert(!info->bytes.test(idx + i));
            info->bytes.set(idx + i);
        }
    }

    // step 2: reduce each transaction size, if possible
    std::map<Address, std::list<transactionInfo>>::iterator t_list;
    for (t_list = warpTransactions.begin();
            t_list != warpTransactions.end(); t_list++)
    {
        // For each block addr
        Address addr = t_list->first;
        const std::list<transactionInfo> &transactionList = t_list->second;

        std::list<transactionInfo>::const_iterator t;
        for (t = transactionList.begin(); t != transactionList.end(); t++)
        {
            // For each transaction
            const transactionInfo &info = *t;
            ReduceTransactionSize((inst.mDecoded.op == EUArchOp::STORE_OP), accessType, info,
                                                    addr, segmentSize, inst.mDecoded.accessQ);
        }
    }
    return 0;
}

// Based on Accel-sim void warp_inst_t::memory_coalescing_arch_reduce_and_send(...)
void CacheModel::ReduceTransactionSize(
    bool isWrite, EMemAccessType accessType, const transactionInfo &info,
    Address addr, unsigned segmentSize, std::list<MemAccess> &accessQ)
{
    assert((addr & (segmentSize - 1)) == 0);

    const std::bitset<4> &q = info.chunks;
    assert(q.count() >= 1);
    std::bitset<2> h; // halves (used to check if 64 byte segment can be
                      // compressed into a single 32 byte segment)

    unsigned size = segmentSize;
    if (segmentSize == 128)
    {
        bool lowerHalfUsed = q[0] || q[1];
        bool upperHalfUsed = q[2] || q[3];
        if (lowerHalfUsed && !upperHalfUsed)
        {
            // only lower 64 bytes used
            size = 64;
            if (q[0])
                h.set(0);
            if (q[1])
                h.set(1);
        }
        else if ((!lowerHalfUsed) && upperHalfUsed)
        {
            // only upper 64 bytes used
            addr = addr + 64;
            size = 64;
            if (q[2])
                h.set(0);
            if (q[3])
                h.set(1);
        }
        else
        {
            assert(lowerHalfUsed && upperHalfUsed);
        }
    }
    else if (segmentSize == 64)
    {
        // need to set halves
        if ((addr % 128) == 0)
        {
            if (q[0])
                h.set(0);
            if (q[1])
                h.set(1);
        }
        else
        {
            assert((addr % 128) == 64);
            if (q[2])
                h.set(0);
            if (q[3])
                h.set(1);
        }
    }
    if (size == 64)
    {
        bool lowerHalfUsed = h[0];
        bool upperHalfUsed = h[1];
        if (lowerHalfUsed && !upperHalfUsed)
        {
            size = 32;
        }
        else if ((!lowerHalfUsed) && upperHalfUsed)
        {
            addr = addr + 32;
            size = 32;
        }
        else
        {
            assert(lowerHalfUsed && upperHalfUsed);
        }
    }
    accessQ.push_back(MemAccess(addr, isWrite, size, accessType, info.chunks));
}

int GCoM::RunSimpleCache(CacheModel *simpleCache, SASSDecoder *decoder, HWConfig hwConfig)
{
    unsigned numWarpPerCTA = decoder->mKernelInfo.numWarps / decoder->mKernelInfo.numCTAs;

    struct CTACtx
    {
        unsigned numDoneWarp = 0;
        unsigned warpInstIdx = 0;
    };
    struct SMCtx
    {
        std::map<int, struct CTACtx> mCtaCtx;
    };
    // smCtxs = 
    // {sm0 : cta0, done warp #, warpInstIdx
    //      cta1, done warp #, warpInstIdx
    //      .., 
    // sm1:  }
    std::vector<SMCtx> smCtxs(hwConfig.numSMs);

    // warpidx2cta (warpIdx) {return warpIdx / numWarpPerCTA }
    // cta2warpidxs (ctaid) {return ctaid*numWarpPerCTA<=  <(ctaid+1)*numWarpPerCTA}
    
    unsigned numFinishedCta = 0, numAssignedCta = 0, nextSmId = 0;
    while (numFinishedCta < decoder->mKernelInfo.numCTAs)
    {
        // assign cta to sm in RR fashion until m
        for (unsigned i = numAssignedCta; i < decoder->mKernelInfo.numCTAs; i++)
        {
            bool assigned = false;
            for (unsigned j = 0; j < hwConfig.numSMs; j++)
            {
                unsigned smId = (nextSmId + j) % hwConfig.numSMs;
                if (smCtxs[smId].mCtaCtx.size() < simpleCache->mMaxCtaPerSm)
                {
                    struct CTACtx new_cta;
                    smCtxs[smId].mCtaCtx.insert({i, new_cta});
                    nextSmId = (smId + 1) % hwConfig.numSMs;
                    numAssignedCta++;
                    assigned = true;
                    break;
                }
            }
            if (!assigned)
                break;
        }

        bool cta_done = false;
        while(!cta_done)
        {
            // for each sm proceed warp inst, access L1 cache
            for (unsigned smId = 0; smId < hwConfig.numSMs; smId++)
            {
                std::vector<int> doneCtaIds;
                for (auto cit = smCtxs[smId].mCtaCtx.begin(); cit != smCtxs[smId].mCtaCtx.end(); cit++)
                {
                    unsigned ctaid = cit->first;
                    unsigned warpInstIdx = cit->second.warpInstIdx;
                    cit->second.numDoneWarp = 0;
                    for (unsigned warpIdx = ctaid * numWarpPerCTA; warpIdx < (ctaid+1) * numWarpPerCTA; warpIdx++)
                    {
                        // if warp is not done
                        if (warpInstIdx < decoder->GetWarp(warpIdx).mInsts.size())
                        {
                            WarpInst &inst = decoder->GetWarpInst(warpIdx, warpInstIdx);
                            if ((inst.mDecoded.op == EUArchOp::LOAD_OP || inst.mDecoded.op == EUArchOp::STORE_OP) && 
                                (inst.mDecoded.space == EMemorySpace::GLOBAL_SPACE || inst.mDecoded.space == EMemorySpace::LOCAL_SPACE))
                            {
                                if (simpleCache->AccessCache(inst, smId, hwConfig, warpInstIdx) != 0)
                                    return -1;
                            }
                        }
                        else
                            cit->second.numDoneWarp++;
                    }
                    cit->second.warpInstIdx++;

                    // when done warp # > numWarpPerCTA break, to assign new cta
                    if (cit->second.numDoneWarp >= numWarpPerCTA)
                    {
                        doneCtaIds.push_back(cit->first);
                        numFinishedCta++;
                        cta_done = true;
                    }
                }
                for (std::vector<int>::iterator it = doneCtaIds.begin(); it != doneCtaIds.end(); it++)
                    smCtxs[smId].mCtaCtx.erase(*it);
            }
        }
    }

    return 0;
}

// Copide from Accel-sim hashing.cc
unsigned HashAdressWithIpolyFunction(Address higherBits, unsigned index, unsigned nBins)
{
    /*
     * Set Indexing function from "Pseudo-randomly interleaved memory."
     * Rau, B. R et al.
     * ISCA 1991
     * http://citeseerx.ist.psu.edu/viewdoc/download;jsessionid=348DEA37A3E440473B3C075EAABC63B6?doi=10.1.1.12.7149&rep=rep1&type=pdf
     *
     * equations are corresponding to IPOLY(37) and are adopted from:
     * "Sacat: streaming-aware conflict-avoiding thrashing-resistant gpgpu
     * cache management scheme." Khairy et al. IEEE TPDS 2017.
     *
     * equations for 16 banks are corresponding to IPOLY(5)
     * equations for 32 banks are corresponding to IPOLY(37)
     * equations for 64 banks are corresponding to IPOLY(67)
     * To see all the IPOLY equations for all the degrees, see
     * http://wireless-systems.ece.gatech.edu/6604/handouts/Peterson's%20Table.pdf
     *
     * We generate these equations using GF(2) arithmetic:
     * http://www.ee.unb.ca/cgi-bin/tervo/calc.pl?num=&den=&f=d&e=1&m=1
     *
     * We go through all the strides 128 (10000000), 256 (100000000),...  and
     * do modular arithmetic in GF(2) Then, we create the H-matrix and group
     * each bit together, for more info read the ISCA 1991 paper
     *
     * IPOLY hashing guarantees conflict-free for all 2^n strides which widely
     * exit in GPGPU applications and also show good performance for other
     * strides.
     */
    if (nBins == 16)
    {
        std::bitset<64> a(higherBits);
        std::bitset<4> b(index);
        std::bitset<4> newIndex(index);

        newIndex[0] =
            a[11] ^ a[10] ^ a[9] ^ a[8] ^ a[6] ^ a[4] ^ a[3] ^ a[0] ^ b[0];
        newIndex[1] =
            a[12] ^ a[8] ^ a[7] ^ a[6] ^ a[5] ^ a[3] ^ a[1] ^ a[0] ^ b[1];
        newIndex[2] = a[9] ^ a[8] ^ a[7] ^ a[6] ^ a[4] ^ a[2] ^ a[1] ^ b[2];
        newIndex[3] = a[10] ^ a[9] ^ a[8] ^ a[7] ^ a[5] ^ a[3] ^ a[2] ^ b[3];

        return newIndex.to_ulong();
    }
    else if (nBins == 32)
    {
        std::bitset<64> a(higherBits);
        std::bitset<5> b(index);
        std::bitset<5> newIndex(index);

        newIndex[0] =
            a[13] ^ a[12] ^ a[11] ^ a[10] ^ a[9] ^ a[6] ^ a[5] ^ a[3] ^ a[0] ^ b[0];
        newIndex[1] = a[14] ^ a[13] ^ a[12] ^ a[11] ^ a[10] ^ a[7] ^ a[6] ^ a[4] ^
                       a[1] ^ b[1];
        newIndex[2] =
            a[14] ^ a[10] ^ a[9] ^ a[8] ^ a[7] ^ a[6] ^ a[3] ^ a[2] ^ a[0] ^ b[2];
        newIndex[3] =
            a[11] ^ a[10] ^ a[9] ^ a[8] ^ a[7] ^ a[4] ^ a[3] ^ a[1] ^ b[3];
        newIndex[4] =
            a[12] ^ a[11] ^ a[10] ^ a[9] ^ a[8] ^ a[5] ^ a[4] ^ a[2] ^ b[4];
        return newIndex.to_ulong();
    }
    else if (nBins == 64)
    {
        std::bitset<64> a(higherBits);
        std::bitset<6> b(index);
        std::bitset<6> newIndex(index);

        newIndex[0] = a[18] ^ a[17] ^ a[16] ^ a[15] ^ a[12] ^ a[10] ^ a[6] ^ a[5] ^
                       a[0] ^ b[0];
        newIndex[1] = a[15] ^ a[13] ^ a[12] ^ a[11] ^ a[10] ^ a[7] ^ a[5] ^ a[1] ^
                       a[0] ^ b[1];
        newIndex[2] = a[16] ^ a[14] ^ a[13] ^ a[12] ^ a[11] ^ a[8] ^ a[6] ^ a[2] ^
                       a[1] ^ b[2];
        newIndex[3] = a[17] ^ a[15] ^ a[14] ^ a[13] ^ a[12] ^ a[9] ^ a[7] ^ a[3] ^
                       a[2] ^ b[3];
        newIndex[4] = a[18] ^ a[16] ^ a[15] ^ a[14] ^ a[13] ^ a[10] ^ a[8] ^ a[4] ^
                       a[3] ^ b[4];
        newIndex[5] =
            a[17] ^ a[16] ^ a[15] ^ a[14] ^ a[11] ^ a[9] ^ a[5] ^ a[4] ^ b[5];
        return newIndex.to_ulong();
    }
    else
    { /* Else incorrect number of bins for the hashing function */
        assert(
            "\nmemory_partition_indexing error: The number of "
            "bins should be "
            "16, 32 or 64 for the hashing IPOLY index function. other bin "
            "numbers are not supported. Generate it by yourself! \n" &&
            0);

        return 0;
    }
}

int GCoM::CalculateMaxConcurrentCTAPerSM(KernelInfo kernelInfo, HWConfig hwConfig, unsigned &maxConcurrentCTAPerSM)
{

    unsigned paddedThreadsPerCTA = kernelInfo.numWarps / kernelInfo.numCTAs;
    paddedThreadsPerCTA += (kernelInfo.numWarps % kernelInfo.numCTAs) ? 1 : 0;
    paddedThreadsPerCTA *= MAX_WARP_SIZE;
    
    // Limit by n_threads / SM
    unsigned resultThread = hwConfig.maxThreadsPerSM / paddedThreadsPerCTA;

    // Limit by shared memory / SM
    unsigned resultSharedMem = (unsigned) -1;
    if (kernelInfo.sharedMemSize > 0)
    {
        // last element of hwConfig.sharedMemSizeOption is the maximum shared memory size
        unsigned maxSharedMemSize = hwConfig.sharedMemSizeOption.back();
        resultSharedMem = maxSharedMemSize / kernelInfo.sharedMemSize;
    }

    // Limit by register count, rounded up to multiple of 4.
    unsigned resultReg = (unsigned) -1;
    if (kernelInfo.numRegs > 0)
        resultReg = hwConfig.registerPerSM / (paddedThreadsPerCTA * ((kernelInfo.numRegs + 3) & ~3));

    // Limit by CTA
    unsigned resultCTA = hwConfig.maxCTAPerSM;

    maxConcurrentCTAPerSM = min({resultThread, resultSharedMem, resultReg, resultCTA});

    if (maxConcurrentCTAPerSM < 1)
    {
        cout << "[Error] Kernel requires more resources than SM" << endl;
        return -1;
    }
    else
        return 0;
}

// Copide from Accel-sim cache_config::hash_function
unsigned GCoM::HashAddress(Address addr, unsigned nBins, unsigned offSetBits, EHashFunction hashFunctionType)
{
    unsigned log2nBins = LogB2(nBins);

    unsigned binIndex = 0;

    switch (hashFunctionType)
    {
    case EHashFunction::FERMI_HASH_SET_FUNCTION:
    {
        /*
         * Set Indexing function from "A Detailed GPU Cache Model Based on Reuse
         * Distance Theory" Cedric Nugteren et al. HPCA 2014
         */
        unsigned lowerXor = 0;
        unsigned upperXor = 0;

        if (nBins == 32 || nBins == 64)
        {
            // Lower xor value is bits 7-11
            lowerXor = (addr >> offSetBits) & 0x1F;

            // Upper xor value is bits 13, 14, 15, 17, and 19
            upperXor = (addr & 0xE000) >> 13;   // Bits 13, 14, 15
            upperXor |= (addr & 0x20000) >> 14; // Bit 17
            upperXor |= (addr & 0x80000) >> 15; // Bit 19

            binIndex = (lowerXor ^ upperXor);

            // 48KB cache prepends the binIndex with bit 12
            if (nBins == 64)
                binIndex |= (addr & 0x1000) >> 7;
        }
        else
        {
            assert(false && "Incorrect number of bins for Fermi hashing function.\n");
        }
        break;
    }
    case EHashFunction::BITWISE_XORING_FUNCTION:
    {
        Address higherBits = addr >> (offSetBits + log2nBins);
        unsigned index = (addr >> offSetBits) & (nBins - 1);
        binIndex = (index) ^ (higherBits & (nBins - 1));
        break;
    }
    case EHashFunction::HASH_IPOLY_FUNCTION:
    {
        Address higherBits = addr >> (offSetBits + log2nBins);
        unsigned index = (addr >> offSetBits) & (nBins - 1);
        binIndex = HashAdressWithIpolyFunction(higherBits, index, nBins);
        break;
    }

    case EHashFunction::LINEAR_SET_FUNCTION:
    {
        binIndex = (addr >> offSetBits) & (nBins - 1);
        break;
    }

    default:
    {
        assert(false && "\nUndefined hash function.\n");
        break;
    }
    }

    assert(binIndex < nBins);

    return binIndex;
}