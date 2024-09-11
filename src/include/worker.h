#ifndef _WORKER_H
#define _WORKER_H

#include <stdlib.h>
#include <iostream>
#include <filesystem>
#include "common.h"

namespace GCoM
{
    class OptionParser
    {
    public:
        void ParseCommandLine(int argc, char **argv);
        void ShowHelp(char *argStr);

        // input: mConfigPath
        // output: hwConfig
        int ReadHWConfig(HWConfig &hwConfig);
    
        // Command line options - general
        int mWorkType;
        std::filesystem::path mRepWarpPath;
        std::filesystem::path mConfigPath;
        std::filesystem::path mSASSTracePath;

        // Command line options - RunGCoMWithSimpleCacheModel
        std::filesystem::path mExportWarpSelInputPath;
    };

    class Worker
    {
    public:
        Worker();
        int DoWork(int argc, char **argv);
        void ShowHelp(char *argStr);
    
    private:
        int (*WorkCases[20])(OptionParser opp);
    };

    int RunGCoMWithSimpleCacheModel(OptionParser opp);
}

#endif