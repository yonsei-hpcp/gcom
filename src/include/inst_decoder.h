#ifndef _INST_DECODER_H
#define _INST_DECODER_H

#include <vector>
#include <unordered_map>
#include <filesystem>
#include "worker.h"
#include "instruction.h"
#include "trace_parser.h"

namespace GCoM
{
    // SASS decoder dependent on SASS trace, but independent of hardware configuration
    class SASSDecoder
    {
    public:
        SASSDecoder(OptionParser opp);
        ~SASSDecoder();

        // input: mAccelsimTracer
        // output: mWarps, kernelInfo
        int ParseKernalTrace();

        WarpInst &GetWarpInst(int warpId, int instIdx); // This will also used to update cache stats in warp inst
        Warp &GetWarp(int warpId);
        int FreeWarps(); // free warps after cache simulation and representative warp selection of a kernel

        void ExportWarps(std::filesystem::path exportPath);

        KernelInfo mKernelInfo;

    private:
        std::vector<Warp> mWarps; // decoded warps for a kernel
        std::filesystem::path mSASSTracePath;

        trace_parser *mAccelsimTracer;
        std::vector<trace_command> mCommandlist;
        std::vector<trace_command>::iterator mCurrentGPUCommand;

        // Input: instTrace, kernel_trace_info
        // Output: warpInst with decoded information, mKernelInfo.numThreadInsts
        int DecocdeInst(inst_trace_t &instTrace, kernel_trace_t kernel_trace_info, WarpInst &warpInst);
        int ClassifyLatencyInitIntv(std::string instTraceOpcode, WarpInst &warpInst);

        
    }; // class SASSDecoder

    // Modified based on Accel-sim gpu-simulator/ISA_Def/turing_opcode.h
    static const std::unordered_map<std::string, EUArchOp> TuringOpcodeMap = {
        // Floating Point 32 Instructions
        {"FADD", EUArchOp::SP_OP},
        {"FADD32I", EUArchOp::SP_OP},
        {"FCHK", EUArchOp::SP_OP},
        {"FFMA32I", EUArchOp::SP_OP},
        {"FFMA", EUArchOp::SP_OP},
        {"FMNMX", EUArchOp::SP_OP},
        {"FMUL", EUArchOp::SP_OP},
        {"FMUL32I", EUArchOp::SP_OP},
        {"FSEL", EUArchOp::SP_OP},
        {"FSET", EUArchOp::SP_OP},
        {"FSETP", EUArchOp::SP_OP},
        {"FSWZADD", EUArchOp::SP_OP},
        // SFU
        {"MUFU", EUArchOp::SFU_OP},

        // Floating Point 16 Instructions
        {"HADD2", EUArchOp::SP_OP},
        {"HADD2_32I", EUArchOp::SP_OP},
        {"HFMA2", EUArchOp::SP_OP},
        {"HFMA2_32I", EUArchOp::SP_OP},
        {"HMUL2", EUArchOp::SP_OP},
        {"HMUL2_32I", EUArchOp::SP_OP},
        {"HSET2", EUArchOp::SP_OP},
        {"HSETP2", EUArchOp::SP_OP},

        // Tensor Core Instructions
        // Execute Tensor Core Instructions on SPECIALIZED_UNIT_3
        {"HMMA", EUArchOp::SPECIALIZED_UNIT_3_OP},
        {"BMMA", EUArchOp::SPECIALIZED_UNIT_3_OP},
        {"IMMA", EUArchOp::SPECIALIZED_UNIT_3_OP},

        // Double Point Instructions
        {"DADD", EUArchOp::DP_OP},
        {"DFMA", EUArchOp::DP_OP},
        {"DMUL", EUArchOp::DP_OP},
        {"DSETP", EUArchOp::DP_OP},

        // Integer Instructions
        {"BMSK", EUArchOp::INTP_OP},
        {"BREV", EUArchOp::INTP_OP},
        {"FLO", EUArchOp::INTP_OP},
        {"IABS", EUArchOp::INTP_OP},
        {"IADD", EUArchOp::INTP_OP},
        {"IADD3", EUArchOp::INTP_OP},
        {"IADD32I", EUArchOp::INTP_OP},
        {"IDP", EUArchOp::INTP_OP},
        {"IDP4A", EUArchOp::INTP_OP},
        {"IMAD", EUArchOp::INTP_OP},
        {"IMNMX", EUArchOp::INTP_OP},
        {"IMUL", EUArchOp::INTP_OP},
        {"IMUL32I", EUArchOp::INTP_OP},
        {"ISCADD", EUArchOp::INTP_OP},
        {"ISCADD32I", EUArchOp::INTP_OP},
        {"ISETP", EUArchOp::INTP_OP},
        {"LEA", EUArchOp::INTP_OP},
        {"LOP", EUArchOp::INTP_OP},
        {"LOP3", EUArchOp::INTP_OP},
        {"LOP32I", EUArchOp::INTP_OP},
        {"POPC", EUArchOp::INTP_OP},
        {"SHF", EUArchOp::INTP_OP},
        {"SHL", EUArchOp::INTP_OP},  //////////
        {"SHR", EUArchOp::INTP_OP},
        {"VABSDIFF", EUArchOp::INTP_OP},
        {"VABSDIFF4", EUArchOp::INTP_OP},

        // Conversion Instructions
        {"F2F", EUArchOp::ALU_OP},
        {"F2I", EUArchOp::ALU_OP},
        {"I2F", EUArchOp::ALU_OP},
        {"I2I", EUArchOp::ALU_OP},
        {"I2IP", EUArchOp::ALU_OP},
        {"FRND", EUArchOp::ALU_OP},

        // Movement Instructions
        {"MOV", EUArchOp::ALU_OP},
        {"MOV32I", EUArchOp::ALU_OP},
        {"MOVM", EUArchOp::ALU_OP},  // move matrix
        {"PRMT", EUArchOp::ALU_OP},
        {"SEL", EUArchOp::ALU_OP},
        {"SGXT", EUArchOp::ALU_OP},
        {"SHFL", EUArchOp::ALU_OP},

        // Predicate Instructions
        {"PLOP3", EUArchOp::ALU_OP},
        {"PSETP", EUArchOp::ALU_OP},
        {"P2R", EUArchOp::ALU_OP},
        {"R2P", EUArchOp::ALU_OP},

        // Load/Store Instructions
        {"LD", EUArchOp::LOAD_OP},
        // For now, we ignore constant loads, consider it as ALU_OP, TO DO
        {"LDC", EUArchOp::ALU_OP},
        {"LDG", EUArchOp::LOAD_OP},
        {"LDL", EUArchOp::LOAD_OP},
        {"LDS", EUArchOp::LOAD_OP},
        {"LDSM", EUArchOp::LOAD_OP},  //
        {"ST", EUArchOp::STORE_OP},
        {"STG", EUArchOp::STORE_OP},
        {"STL", EUArchOp::STORE_OP},
        {"STS", EUArchOp::STORE_OP},
        {"MATCH", EUArchOp::ALU_OP},
        {"QSPC", EUArchOp::ALU_OP},
        {"ATOM", EUArchOp::LOAD_OP}, // Reference: Accelsim bool trace_warp_inst_t::parse_from_trace_struct
        {"ATOMS", EUArchOp::STORE_OP},
        {"ATOMG", EUArchOp::LOAD_OP}, // Reference: Accelsim bool trace_warp_inst_t::parse_from_trace_struct
        {"RED", EUArchOp::LOAD_OP}, // Reference: Accelsim bool trace_warp_inst_t::parse_from_trace_struct
        {"CCTL", EUArchOp::ALU_OP},
        {"CCTLL", EUArchOp::ALU_OP},
        {"ERRBAR", EUArchOp::ALU_OP},
        {"MEMBAR", EUArchOp::MEMORY_BARRIER_OP},
        {"CCTLT", EUArchOp::ALU_OP},

        // Uniform Datapath Instruction
        // UDP unit
        // for more info about UDP, see
        // https://www.hotchips.org/hc31/HC31_2.12_NVIDIA_final.pdf
        {"R2UR", EUArchOp::SPECIALIZED_UNIT_4_OP},
        {"S2UR", EUArchOp::SPECIALIZED_UNIT_4_OP},
        {"UBMSK", EUArchOp::SPECIALIZED_UNIT_4_OP},
        {"UBREV", EUArchOp::SPECIALIZED_UNIT_4_OP},
        {"UCLEA", EUArchOp::SPECIALIZED_UNIT_4_OP},
        {"UFLO", EUArchOp::SPECIALIZED_UNIT_4_OP},
        {"UIADD3", EUArchOp::SPECIALIZED_UNIT_4_OP},
        {"UIMAD", EUArchOp::SPECIALIZED_UNIT_4_OP},
        {"UISETP", EUArchOp::SPECIALIZED_UNIT_4_OP},
        {"ULDC", EUArchOp::SPECIALIZED_UNIT_4_OP},
        {"ULEA", EUArchOp::SPECIALIZED_UNIT_4_OP},
        {"ULOP", EUArchOp::SPECIALIZED_UNIT_4_OP},
        {"ULOP3", EUArchOp::SPECIALIZED_UNIT_4_OP},
        {"ULOP32I", EUArchOp::SPECIALIZED_UNIT_4_OP},
        {"UMOV", EUArchOp::SPECIALIZED_UNIT_4_OP},
        {"UP2UR", EUArchOp::SPECIALIZED_UNIT_4_OP},
        {"UPLOP3", EUArchOp::SPECIALIZED_UNIT_4_OP},
        {"UPOPC", EUArchOp::SPECIALIZED_UNIT_4_OP},
        {"UPRMT", EUArchOp::SPECIALIZED_UNIT_4_OP},
        {"UPSETP", EUArchOp::SPECIALIZED_UNIT_4_OP},
        {"UR2UP", EUArchOp::SPECIALIZED_UNIT_4_OP},
        {"USEL", EUArchOp::SPECIALIZED_UNIT_4_OP},
        {"USGXT", EUArchOp::SPECIALIZED_UNIT_4_OP},
        {"USHF", EUArchOp::SPECIALIZED_UNIT_4_OP},
        {"USHL", EUArchOp::SPECIALIZED_UNIT_4_OP},
        {"USHR", EUArchOp::SPECIALIZED_UNIT_4_OP},
        {"VOTEU", EUArchOp::SPECIALIZED_UNIT_4_OP},

        // Texture Instructions
        // For now, we ignore texture loads, consider it as ALU_OP
        {"TEX", EUArchOp::SPECIALIZED_UNIT_2_OP},
        {"TLD", EUArchOp::SPECIALIZED_UNIT_2_OP},
        {"TLD4", EUArchOp::SPECIALIZED_UNIT_2_OP},
        {"TMML", EUArchOp::SPECIALIZED_UNIT_2_OP},
        {"TXD", EUArchOp::SPECIALIZED_UNIT_2_OP},
        {"TXQ", EUArchOp::SPECIALIZED_UNIT_2_OP},

        // Surface Instructions //
        {"SUATOM", EUArchOp::ALU_OP},
        {"SULD", EUArchOp::ALU_OP},
        {"SURED", EUArchOp::ALU_OP},
        {"SUST", EUArchOp::ALU_OP},

        // Control Instructions
        // execute branch insts on a dedicated branch unit (SPECIALIZED_UNIT_1)
        {"BMOV", EUArchOp::SPECIALIZED_UNIT_1_OP},
        {"BPT", EUArchOp::SPECIALIZED_UNIT_1_OP},
        {"BRA", EUArchOp::SPECIALIZED_UNIT_1_OP},
        {"BREAK", EUArchOp::SPECIALIZED_UNIT_1_OP},
        {"BRX", EUArchOp::SPECIALIZED_UNIT_1_OP},
        {"BRXU", EUArchOp::SPECIALIZED_UNIT_1_OP},  //
        {"BSSY", EUArchOp::SPECIALIZED_UNIT_1_OP},
        {"BSYNC", EUArchOp::SPECIALIZED_UNIT_1_OP},
        {"CALL", EUArchOp::SPECIALIZED_UNIT_1_OP},
        {"EXIT", EUArchOp::EXIT_OPS},
        {"JMP", EUArchOp::SPECIALIZED_UNIT_1_OP},
        {"JMX", EUArchOp::SPECIALIZED_UNIT_1_OP},
        {"JMXU", EUArchOp::SPECIALIZED_UNIT_1_OP},  ///
        {"KILL", EUArchOp::SPECIALIZED_UNIT_3_OP},
        {"NANOSLEEP", EUArchOp::SPECIALIZED_UNIT_1_OP},
        {"RET", EUArchOp::SPECIALIZED_UNIT_1_OP},
        {"RPCMOV", EUArchOp::SPECIALIZED_UNIT_1_OP},
        {"RTT", EUArchOp::SPECIALIZED_UNIT_1_OP},
        {"WARPSYNC", EUArchOp::SPECIALIZED_UNIT_1_OP},
        {"YIELD", EUArchOp::SPECIALIZED_UNIT_1_OP},

        // Miscellaneous Instructions
        {"B2R", EUArchOp::ALU_OP},
        {"BAR", EUArchOp::BARRIER_OP},
        {"CS2R", EUArchOp::ALU_OP},
        {"CSMTEST", EUArchOp::ALU_OP},
        {"DEPBAR", EUArchOp::ALU_OP},
        {"GETLMEMBASE", EUArchOp::ALU_OP},
        {"LEPC", EUArchOp::ALU_OP},
        {"NOP", EUArchOp::ALU_OP},
        {"PMTRIG", EUArchOp::ALU_OP},
        {"R2B", EUArchOp::ALU_OP},
        {"S2R", EUArchOp::ALU_OP},
        {"SETCTAID", EUArchOp::ALU_OP},
        {"SETLMEMBASE", EUArchOp::ALU_OP},
        {"VOTE", EUArchOp::ALU_OP},
        {"VOTE_VTG", EUArchOp::ALU_OP}
    };

    // Modified based on Accel-sim gpu-simulator/ISA_Def/volta_opcode.h
    static const std::unordered_map<std::string, EUArchOp> VoltaOpcodeMap = {
        // Floating Point 32 Instructions
        {"FADD", EUArchOp::SP_OP},
        {"FADD32I", EUArchOp::SP_OP},
        {"FCHK", EUArchOp::SP_OP},
        {"FFMA32I", EUArchOp::SP_OP},
        {"FFMA", EUArchOp::SP_OP},
        {"FMNMX", EUArchOp::SP_OP},
        {"FMUL", EUArchOp::SP_OP},
        {"FMUL32I", EUArchOp::SP_OP},
        {"FSEL", EUArchOp::SP_OP},
        {"FSET", EUArchOp::SP_OP},
        {"FSETP", EUArchOp::SP_OP},
        {"FSWZADD", EUArchOp::SP_OP},
        // SFU
        {"MUFU", EUArchOp::SFU_OP},

        // Floating Point 16 Instructions
        {"HADD2", EUArchOp::SP_OP},
        {"HADD2_32I", EUArchOp::SP_OP},
        {"HFMA2", EUArchOp::SP_OP},
        {"HFMA2_32I", EUArchOp::SP_OP},
        {"HMUL2", EUArchOp::SP_OP},
        {"HMUL2_32I", EUArchOp::SP_OP},
        {"HSET2", EUArchOp::SP_OP},
        {"HSETP2", EUArchOp::SP_OP},

        // Tensor Core Instructions
        // Execute Tensor Core Instructions on SPECIALIZED_UNIT_3
        {"HMMA", EUArchOp::SPECIALIZED_UNIT_3_OP},

        // Double Point Instructions
        {"DADD", EUArchOp::DP_OP},
        {"DFMA", EUArchOp::DP_OP},
        {"DMUL", EUArchOp::DP_OP},
        {"DSETP", EUArchOp::DP_OP},

        // Integer Instructions
        {"BMSK", EUArchOp::INTP_OP},
        {"BREV", EUArchOp::INTP_OP},
        {"FLO", EUArchOp::INTP_OP},
        {"IABS", EUArchOp::INTP_OP},
        {"IADD", EUArchOp::INTP_OP},
        {"IADD3", EUArchOp::INTP_OP},
        {"IADD32I", EUArchOp::INTP_OP},
        {"IDP", EUArchOp::INTP_OP},
        {"IDP4A", EUArchOp::INTP_OP},
        {"IMAD", EUArchOp::INTP_OP},
        {"IMMA", EUArchOp::INTP_OP},
        {"IMNMX", EUArchOp::INTP_OP},
        {"IMUL", EUArchOp::INTP_OP},
        {"IMUL32I", EUArchOp::INTP_OP},
        {"ISCADD", EUArchOp::INTP_OP},
        {"ISCADD32I", EUArchOp::INTP_OP},
        {"ISETP", EUArchOp::INTP_OP},
        {"LEA", EUArchOp::INTP_OP},
        {"LOP", EUArchOp::INTP_OP},
        {"LOP3", EUArchOp::INTP_OP},
        {"LOP32I", EUArchOp::INTP_OP},
        {"POPC", EUArchOp::INTP_OP},
        {"SHF", EUArchOp::INTP_OP},
        {"SHR", EUArchOp::INTP_OP},
        {"VABSDIFF", EUArchOp::INTP_OP},
        {"VABSDIFF4", EUArchOp::INTP_OP},

        // Conversion Instructions
        {"F2F", EUArchOp::ALU_OP},
        {"F2I", EUArchOp::ALU_OP},
        {"I2F", EUArchOp::ALU_OP},
        {"I2I", EUArchOp::ALU_OP},
        {"I2IP", EUArchOp::ALU_OP},
        {"FRND", EUArchOp::ALU_OP},

        // Movement Instructions
        {"MOV", EUArchOp::ALU_OP},
        {"MOV32I", EUArchOp::ALU_OP},
        {"PRMT", EUArchOp::ALU_OP},
        {"SEL", EUArchOp::ALU_OP},
        {"SGXT", EUArchOp::ALU_OP},
        {"SHFL", EUArchOp::ALU_OP},

        // Predicate Instructions
        {"PLOP3", EUArchOp::ALU_OP},
        {"PSETP", EUArchOp::ALU_OP},
        {"P2R", EUArchOp::ALU_OP},
        {"R2P", EUArchOp::ALU_OP},

        // Load/Store Instructions
        {"LD", EUArchOp::LOAD_OP},
        // For now, we ignore constant loads, consider it as ALU_OP, TO DO
        {"LDC", EUArchOp::ALU_OP},
        {"LDG", EUArchOp::LOAD_OP},
        {"LDL", EUArchOp::LOAD_OP},
        {"LDS", EUArchOp::LOAD_OP},
        {"ST", EUArchOp::STORE_OP},
        {"STG", EUArchOp::STORE_OP},
        {"STL", EUArchOp::STORE_OP},
        {"STS", EUArchOp::STORE_OP},
        {"MATCH", EUArchOp::ALU_OP},
        {"QSPC", EUArchOp::ALU_OP},
        {"ATOM", EUArchOp::LOAD_OP}, // Reference: Accelsim bool trace_warp_inst_t::parse_from_trace_struct
        {"ATOMS", EUArchOp::STORE_OP},
        {"ATOMG", EUArchOp::LOAD_OP}, // Reference: Accelsim bool trace_warp_inst_t::parse_from_trace_struct
        {"RED", EUArchOp::LOAD_OP}, // Reference: Accelsim bool trace_warp_inst_t::parse_from_trace_struct
        {"CCTL", EUArchOp::ALU_OP},
        {"CCTLL", EUArchOp::ALU_OP},
        {"ERRBAR", EUArchOp::ALU_OP},
        {"MEMBAR", EUArchOp::MEMORY_BARRIER_OP},
        {"CCTLT", EUArchOp::ALU_OP},

        // Texture Instructions
        // For now, we ignore texture loads, consider it as ALU_OP
        {"TEX", EUArchOp::SPECIALIZED_UNIT_2_OP},
        {"TLD", EUArchOp::SPECIALIZED_UNIT_2_OP},
        {"TLD4", EUArchOp::SPECIALIZED_UNIT_2_OP},
        {"TMML", EUArchOp::SPECIALIZED_UNIT_2_OP},
        {"TXD", EUArchOp::SPECIALIZED_UNIT_2_OP},
        {"TXQ", EUArchOp::SPECIALIZED_UNIT_2_OP},

        // Control Instructions
        // execute branch insts on a dedicated branch unit (SPECIALIZED_UNIT_1)
        {"BMOV", EUArchOp::SPECIALIZED_UNIT_1_OP},
        {"BPT", EUArchOp::SPECIALIZED_UNIT_1_OP},
        {"BRA", EUArchOp::SPECIALIZED_UNIT_1_OP},
        {"BREAK", EUArchOp::SPECIALIZED_UNIT_1_OP},
        {"BRX", EUArchOp::SPECIALIZED_UNIT_1_OP},
        {"BSSY", EUArchOp::SPECIALIZED_UNIT_1_OP},
        {"BSYNC", EUArchOp::SPECIALIZED_UNIT_1_OP},
        {"CALL", EUArchOp::SPECIALIZED_UNIT_1_OP},
        {"EXIT", EUArchOp::EXIT_OPS},
        {"JMP", EUArchOp::SPECIALIZED_UNIT_1_OP},
        {"JMX", EUArchOp::SPECIALIZED_UNIT_1_OP},
        {"KILL", EUArchOp::SPECIALIZED_UNIT_1_OP},
        {"NANOSLEEP", EUArchOp::SPECIALIZED_UNIT_1_OP},
        {"RET", EUArchOp::SPECIALIZED_UNIT_1_OP},
        {"RPCMOV", EUArchOp::SPECIALIZED_UNIT_1_OP},
        {"RTT", EUArchOp::SPECIALIZED_UNIT_1_OP},
        {"WARPSYNC", EUArchOp::SPECIALIZED_UNIT_1_OP},
        {"YIELD", EUArchOp::SPECIALIZED_UNIT_1_OP},

        // Miscellaneous Instructions
        {"B2R", EUArchOp::ALU_OP},
        {"BAR", EUArchOp::BARRIER_OP},
        {"CS2R", EUArchOp::ALU_OP},
        {"CSMTEST", EUArchOp::ALU_OP},
        {"DEPBAR", EUArchOp::ALU_OP},
        {"GETLMEMBASE", EUArchOp::ALU_OP},
        {"LEPC", EUArchOp::ALU_OP},
        {"NOP", EUArchOp::ALU_OP},
        {"PMTRIG", EUArchOp::ALU_OP},
        {"R2B", EUArchOp::ALU_OP},
        {"S2R", EUArchOp::ALU_OP},
        {"SETCTAID", EUArchOp::ALU_OP},
        {"SETLMEMBASE", EUArchOp::ALU_OP},
        {"VOTE", EUArchOp::ALU_OP},
        {"VOTE_VTG", EUArchOp::ALU_OP}
    };
}

#endif