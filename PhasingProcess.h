#ifndef PHASINGPROCESS_H
#define PHASINGPROCESS_H
#include "Util.h"



struct PhasingParameters
{
    int numThreads;
    int distance;
    std::string snpFile;
    std::string svFile;
    std::vector<std::string> bamFile;
    std::string modFile="";
    std::string ponFile="";
    std::string strictPonFile="";
    std::string fastaFile;
    std::string resultPrefix;
    std::string callerStr;
    Caller caller;
    bool generateDot;
    bool phaseIndel;
    bool disablePonTag=true;
    bool disableCalling=false;
    bool disableRefineSomatic=false;
    bool outputLOH = false;  // Whether to output LOH results
    bool outputSGE = false;  // Whether to output SmallGenomicEvent results
    bool outputLGE = false;  // Whether to output LargeGenomicEvent results
    bool outputGE = false;   // Whether to output GenomicEvent results
    
    int connectAdjacent;
    int mappingQuality;
    double mismatchRate;
    
    int baseQuality;
    double edgeWeight;
    
    double snpConfidence;
    double readConfidence;
    
    double edgeThreshold;
    double overlapThreshold;
    
    int somaticConnectAdjacent;
    
    std::string version;
    std::string command;
    
    // If negative, purity is not provided by user and should be estimated
    double purity = -1.0;
};

class PhasingProcess
{

    public:
        PhasingProcess(PhasingParameters params);
        ~PhasingProcess();

};


#endif
