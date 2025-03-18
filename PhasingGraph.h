#ifndef PHASINGGRAPH_H
#define PHASINGGRAPH_H

#include "Util.h"
#include "ParsingBam.h"
#include "PhasingProcess.h"


typedef std::pair<int, int> PosAllele;
typedef std::map<std::string, int> ReadBaseMap;

class SubEdge{
    
    private:
        int readCount;
        // Edge information. The vector store next pos
        // < next position, read name >
        std::map<int, std::vector<std::string> > *refRead;
        std::map<int, std::vector<std::string> > *altRead;
        // sum of edge pair quality, pos1 quality + pos2 quality
        // < next position, quality sum >
        std::map<int, int> *refQuality;
        std::map<int, int> *altQuality;
        // < next position, read count >
        std::map<int, float> *refReadCount;
        std::map<int, float> *altReadCount;

    public:
        
        SubEdge();
        ~SubEdge();
        
        void destroy();
        
        void addSubEdge(Variant &currentNode, Variant &connectNode, std::string readName, int baseQuality, double edgeWeight, bool fakeRead);
        std::pair<float,float> BestPair(int targetPos);
        float getRefReadCount(int targetPos);
        float getAltReadCount(int targetPos);        
        
        std::vector<std::string> showEdge(std::string message);
        std::vector<std::pair<int,int>> getConnectPos();
 
        int getQuality(PosAllele targetPos);
        int getAvgQuality(PosAllele targetPos);

};

//use to store the voting info from the previous variants
struct VoteResult{
    int Pos;		//who votes
    float para;		//rr+aa
    float cross;	//ra+ar
    float weight;	//how much weight
    int hap;		//which haplotype 
    double ESR;		//similarity of para and cross

    VoteResult( int currPos, float weight ) ;
};

struct VariantInfo {
    VariantType type;
    // bool homozygous; // 0 heterozygous 1 homozygous
    // bool somatic; // 0 germline 1 somatic
};

struct VariantEdge{
    int currPos;
    SubEdge* alt;
    SubEdge* ref;
    int refcnt ; // count the ref base amount
    int altcnt ; // count the alt base amount
    double vaf ; // count the vaf of the left snp
    int coverage ; // count the coverge on the snp
    
    VariantEdge(int currPos);
    // node pair 
    std::pair<PosAllele,PosAllele> findBestEdgePair(std::map<int, VariantInfo>::iterator currNodeIter, std::map<int, VariantInfo>::iterator nextNodeIter, bool isONT, double diffRatioThreshold, VoteResult &vote, bool debug);
    // number of read of two node. AA and AB combination
    std::pair<float,float> findNumberOfRead(int targetPos);
    bool get_fakeSnp();
};


struct BlockRead{
    std::map<std::string,int> readVec;
    
    void recordRead(std::string readName);
};

struct EdgeResult{
    int rr;
    int ra;
    int ar;
    int aa;
};

class VairiantGraph{
    
    private:
        PhasingParameters *params;
        std::string *ref;
        std::vector<std::string> dotResult;
        std::vector<ReadVariant> *readVariant;
        
        // By default, a Map in C++ is sorted in increasing order based on its key.
        // position, edge
        std::map<int,VariantEdge*> *edgeList;
        // position, type
        std::map<int,VariantInfo> *variantPosType;

        // position phasing result
        PosPhasingResult *posPhasingResult;
        // store phased read and read's haplotype
        std::map<std::string,int> *readHpMap;

        // produce PS tag and determine phased GT tag
        void storeResultPath();
        
        void readCorrection();

        void edgeConnectResult();
        
        std::pair<float,float> Onelongcase( std::vector<VoteResult> vote ) ;

    public:
    
        VairiantGraph(std::string &ref, PhasingParameters &params);
        ~VairiantGraph();
    
        void addEdge(std::vector<ReadVariant> &in_readVariant);
        
        void phasingProcess(PosPhasingResult &posPhasingResult);

        void exportPhasingResult(PosPhasingResult &posPhasingResult);
        
        void writingDotFile(std::string dotPrefix);
        std::map<std::string,int>* getReadHP();
        int totalNode();

        void destroy();
        
};

class Roller {
public:
    Roller();
    Roller(size_t windowSize);
    Roller(size_t windowSize, const std::vector<double>& input);
    
    Roller& setValues(const std::vector<double>& input);
    Roller& smooth();
    Roller& forward();
    Roller& reverse();
    Roller& downSampling(int distance);
    Roller& directionalDifference(size_t distance = 10);
    Roller& opposite();
    Roller& reverseNetValues();
    // Get the final result
    std::vector<double> getResult() const;
    
    std::vector<double> smooth(const std::vector<double>& values);
    std::vector<double> forward(const std::vector<double>& values);
    std::vector<double> reverse(const std::vector<double>& values);
    std::vector<double> downSampling(const std::vector<double>& values, int distance);
    std::vector<double> directionalDifference(const std::vector<double>& values, size_t distance = 10);
    std::vector<double> opposite(const std::vector<double>& values);
    std::vector<double> reverseNetValues(const std::vector<double>& values);

    static std::vector<double> backValues(const std::vector<size_t>& idxs, const std::vector<double>& values);
    static std::vector<int> backPosition(const std::vector<size_t>& idxs, const std::vector<int>& values);
    static void replaceValue(ClipCount &clipCount, std::vector<int>& keys, std::vector<size_t>& idxs, std::vector<double>& replaceValues);

private:
    size_t window;
    std::vector<double> values;
    void updateWindow(double& rollingSum, std::deque<double>& windowValues, double newValue);
};

class PointFinder {
public:
    PointFinder(size_t peakDistance, int calibrationThreshold = 5, size_t calibrationDistance = 0);
    
    std::vector<size_t> findAllPeaks(const std::vector<double>& values, double peakThreshold = 0);
    std::vector<size_t> findPeaks(const std::vector<double>& values, double peakThreshold = 0);
    std::vector<size_t> getForwardGentle(const std::vector<double>& values, std::vector<size_t> peaks);
    std::vector<size_t> getReverseGentle(const std::vector<double>& values, std::vector<size_t> peaks);

private:
    size_t peakDistance;
    int calibrationThreshold;
    size_t calibrationDistance;

    size_t start(size_t target, size_t distance, size_t distanceDivisor = 1);
    size_t end(size_t target, size_t distance, size_t size, size_t distanceDivisor = 1);
};

class Clip{
    private:
        std::string chr;
        std::vector<int> *largeGenomicEventInterval;
        std::vector<std::pair<int, int>> *smallGenomicEventRegion;
        void amplifyGentleClip(ClipCount &clipCount);
        void detectInterval(ClipCount &clipCount);
    public:
        Clip(std::string &chr);
        ~Clip();
        void detectGenomicEventInterval(ClipCount &clipCount, std::vector<int> &largeGenomicEventInterval, std::vector<std::pair<int, int>> &smallGenomicEventRegion);
};

#endif
