#include "PhasingGraph.h"

//SubEdge

SubEdge::SubEdge():readCount(0){ 
    refRead = new std::map<int, std::vector<std::string> >;
    altRead = new std::map<int, std::vector<std::string> >;
    refQuality = new std::map<int, int>;
    altQuality = new std::map<int, int>;
    refReadCount = new std::map<int, float>;
    altReadCount = new std::map<int, float>;
}

SubEdge::~SubEdge(){ 
}

void SubEdge::destroy(){
    delete refRead;
    delete altRead;
    delete refQuality;
    delete altQuality;
    delete refReadCount;
    delete altReadCount;
}

void SubEdge::addSubEdge(Variant &currentNode, Variant &connectNode, std::string readName, int conditionQuality, double lowQualityWeight, bool fakeRead){
    double edgeWeight = 1;
    //if the base quality on both snps is high enough and didn't be marked as fakeRead, the edge has normal weight 
    if(fakeRead){
        edgeWeight = 0.01;
    }
    else if( (currentNode.quality > VARIANT_UNDEFINED && currentNode.quality < conditionQuality ) || (connectNode.quality > VARIANT_UNDEFINED && connectNode.quality < conditionQuality ) ){
        edgeWeight = lowQualityWeight;
    }

    // target noded is REF allele
    if(connectNode.allele == 0 ){
        // debug, this parameter will record the names of all reads between two podoubles
        //(*refRead)[connectNode.position].push_back(readName);
        // quality sum
        // (*refQuality)[connectNode.position] += currentQuality + connectNode.quality;
        (*refReadCount)[connectNode.position] += edgeWeight;
    }
    // target noded is ALT allele
    else if(connectNode.allele == 1 ){
        // debug, this parameter will record the names of all reads between two points
        // (*altRead)[connectNode.position].push_back(readName);
        // quality sum
        // (*refQuality)[connectNode.position] += currentQuality + connectNode.quality;
        (*altReadCount)[connectNode.position] += edgeWeight;
    }
    // readCount++;
}

std::pair<float,float> SubEdge::BestPair(int targetPos){
    return std::make_pair( getRefReadCount(targetPos), getAltReadCount(targetPos) );
}

float SubEdge::getRefReadCount(int targetPos){
    std::map<int, float>::iterator posIter = refReadCount->find(targetPos);
    if( posIter != refReadCount->end() ){
        return posIter->second;
    }
    return 0.0f;
}

float SubEdge::getAltReadCount(int targetPos){
    std::map<int, float>::iterator posIter = altReadCount->find(targetPos);
    if( posIter != altReadCount->end() ){
        return posIter->second;
    }
    return 0.0f;
}

std::vector<std::string> SubEdge::showEdge(std::string message){
    std::vector<std::string> result;
    for(std::map<int, float >::iterator edgeIter = refReadCount->begin() ; edgeIter != refReadCount->end() ; edgeIter++ ){
        result.push_back(message +" -> ref_" + std::to_string((*edgeIter).first) + "[label=" + std::to_string((*edgeIter).second) + "];");
    }
    for(std::map<int, float >::iterator edgeIter = altReadCount->begin() ; edgeIter != altReadCount->end() ; edgeIter++ ){
        result.push_back(message +" -> alt_" + std::to_string((*edgeIter).first) + "[label=" + std::to_string((*edgeIter).second) + "];");
    }
    return result;
}

std::vector<std::pair<int,int>> SubEdge::getConnectPos(){
    std::vector<std::pair<int,int>> result;
    for(std::map<int, float >::iterator edgeIter = refReadCount->begin() ; edgeIter != refReadCount->end() ; edgeIter++ ){
        result.push_back( std::make_pair( (*edgeIter).first, 0 ) );
    }
    for(std::map<int, float >::iterator edgeIter = altReadCount->begin() ; edgeIter != altReadCount->end() ; edgeIter++ ){
        result.push_back( std::make_pair( (*edgeIter).first, 1 ) );
    }
    return result;
}

int SubEdge::getQuality(PosAllele targetPos){
    // target is Ref allele
    if( targetPos.second == 1 ){
        std::map<int, int>::iterator qIter = refQuality->find(targetPos.first);
        if( qIter == refQuality->end() )
            return 0;
        else
            return (*refQuality)[targetPos.first];
    }
    // target is Alt allele
    if( targetPos.second == 2 ){
        std::map<int, int>::iterator qIter = altQuality->find(targetPos.first);
        if( qIter == altQuality->end() )
            return 0;
        else
            return (*altQuality)[targetPos.first];
    }
    return 0;
}

int SubEdge::getAvgQuality(PosAllele targetPos){
    // target is Ref allele
    if( targetPos.second == 1 ){
        std::map<int, int>::iterator qIter = refQuality->find(targetPos.first);
        if( qIter == refQuality->end() )
            return 0;
        else
            return (*refQuality)[targetPos.first]/(*refReadCount)[targetPos.first];
    }
    // target is Alt allele
    if( targetPos.second == 2 ){
        std::map<int, int>::iterator qIter = altQuality->find(targetPos.first);
        if( qIter == altQuality->end() )
            return 0;
        else
            return (*altQuality)[targetPos.first]/(*altReadCount)[targetPos.first];
    }
    return 0;
}

VoteResult::VoteResult( int currPos, float variantweight ) {
    Pos = currPos ;
    weight = variantweight ;
}

VariantEdge::VariantEdge(int inCurrPos){
    currPos = inCurrPos;
    alt = new SubEdge();
    ref = new SubEdge();
    edgeCount = new std::map<int, std::map<int, ThreePointEdge>>;
    refcnt = 0; 
    altcnt = 0; 
    coverage = 0;  
}

void VariantEdge::addSegmentEdge(Variant &currentNode, Variant &connectNode, Variant &connectSecondNode){
    int alleleCombination = (currentNode.allele << 2) | (connectNode.allele << 1) | connectSecondNode.allele;
    (*edgeCount)[connectNode.position][connectSecondNode.position][alleleCombination]++;
}

ThreePointEdge VariantEdge::findThreePointEdge(int pos, int nextPos){
    auto iter = edgeCount->find(pos);
    if(iter != edgeCount->end()){
        auto innerIter = iter->second.find(nextPos);
        if(innerIter != iter->second.end()){
            return innerIter->second;
        }
    }
    return ThreePointEdge();
}

//to get the value of fakeSnp
bool VariantEdge::get_fakeSnp(){
    bool fakeSnp;
    if(vaf == 0 || vaf == 1)
      fakeSnp = true;
    else
      fakeSnp = false;
    
    return fakeSnp;  
}  

//VariantEdge
std::pair<PosAllele,PosAllele> VariantEdge::findBestEdgePair(std::map<int, VariantInfo>::iterator currNodeIter, std::map<int, VariantInfo>::iterator nextNodeIter, bool isONT, double edgeThreshold, VoteResult &vote, bool debug){
    int targetPos = nextNodeIter->first;
    std::pair<float,float> refBestPair  = ref->BestPair(targetPos);
    std::pair<float,float> altBestPair  = alt->BestPair(targetPos);
    // get the weight of each pair
    float rr = refBestPair.first;
    float ra = refBestPair.second;
    float ar = altBestPair.first;
    float aa = altBestPair.second;
    float para;
    float cross;
    if(currNodeIter->second.origin != SOMATIC && nextNodeIter->second.origin != SOMATIC){
        para = rr + aa;
        cross = ra + ar;
    }else if(currNodeIter->second.origin == SOMATIC && nextNodeIter->second.origin != SOMATIC){
        para = aa;
        cross = ar;
    }else if(currNodeIter->second.origin != SOMATIC && nextNodeIter->second.origin == SOMATIC){
        para = aa;
        cross = ra;
    }else if(currNodeIter->second.origin == SOMATIC && nextNodeIter->second.origin == SOMATIC){
        para = aa;
        cross = 0;
    }
    // initialize the edge connection
    // -1 : not connect
    int refAllele = -1;
    int altAllele = -1;
    
    double edgeSimilarRatio = (double)std::min((para),(cross)) / (double)std::max((para),(cross));
    vaf = (float)altcnt/(refcnt+altcnt);
    
    //std::cout << currPos+1 << "\t" << altcnt << "\t" << refcnt << "\t" << vaf << "\n" ;
    if( para > cross ){
        // RR conect
        refAllele = 1;
        altAllele = 2;
    }
    else if( para < cross ){
        // RA connect
        refAllele = 2;
        altAllele = 1;
    }
    else if( para == cross ){
        // no connect 
        // not sure which is better
    }

    if((currNodeIter->second.type == SNP && (nextNodeIter->second.type == MOD_FORWARD_STRAND || nextNodeIter->second.type == MOD_REVERSE_STRAND)) ||
       ((currNodeIter->second.type == MOD_FORWARD_STRAND || currNodeIter->second.type == MOD_REVERSE_STRAND) && nextNodeIter->second.type == SNP)){
        edgeThreshold = 0.3;
        if((para+cross) < 1){
            edgeThreshold = -1;
        }
    }

    if( edgeSimilarRatio > edgeThreshold ){
        refAllele = -1;
        altAllele = -1;
    }
    

    if(debug){
        std::cout << currPos << "\t->\t" << targetPos << "\t|rr aa | ra ar\t" << "\t" << rr << "\t" << aa << "\t" << ra << "\t" << ar  << "\n";
    }

    // if the vaf is 0 or 1, we think this variant is a fake variant and we lower its weight
    if ( vaf == 0 || vaf == 1 ) {
        vote.weight = 0.01 ;
        //std::cout<< "fakesnp\t" << currPos+1 << "->" << targetPos+1 << "\t" << vote.weight << "\n";
    }
    // the lower the edgeSimilarRatio means the higher reads consistency, and we will make the weight bigger if the reads consistency is high enough
    else if ( (edgeSimilarRatio <= 0.1 && (para + cross) >= 1)  || ((para)<1&&(cross)>=1) || ((para)>=1&&(cross)<1) ) {
        vote.weight = 20 ;
    }

    vote.para = para ;
    vote.cross = cross ;
    vote.ESR = edgeSimilarRatio ;

    // create edge pairs
    PosAllele refEdge = std::make_pair( targetPos, refAllele );
    PosAllele altEdge = std::make_pair( targetPos, altAllele );
    // return edge pair
    return std::make_pair( refEdge, altEdge );
}

int VairiantGraph::patternMining(ThreePointEdge threePointEdge){
    constexpr float condition = 2;
    SomaticVote vote = VOTE_UNDEFINED;

    float largest = -1;
    float secondLargest = -1;
    float thirdLargest = -1;

    for (float value : threePointEdge) {
        if (value > largest) {
            thirdLargest = secondLargest;
            secondLargest = largest;
            largest = value;
        } else if (value > secondLargest) {
            thirdLargest = secondLargest;
            secondLargest = value;
        } else if (value > thirdLargest) {
            thirdLargest = value;
        }
    }

    if (secondLargest == 0.0f) {
        return VOTE_UNDEFINED;
    }
    // find the vote of the third path
    if (thirdLargest > 0.0f) {
        for (auto &pattern : highSomaticPatterns) {
            const auto &edge = pattern.edges;
            if (threePointEdge[edge[0]] >= thirdLargest &&
                threePointEdge[edge[1]] >= thirdLargest &&
                threePointEdge[edge[0]] >= condition &&
                threePointEdge[edge[1]] >= condition &&
                threePointEdge[edge[2]] >= thirdLargest) {
                vote = pattern.vote;
                break;
            }
        }
    }
    // find the vote of the second path
    if(vote == VOTE_UNDEFINED && secondLargest/2 > thirdLargest){
        for(auto &pattern : lowSomaticPatterns){
            const auto &edge = pattern.edges;
            if (threePointEdge[edge[1]] >= secondLargest && 
                threePointEdge[edge[2]] >= secondLargest) {
                vote = pattern.vote;
                break;
            }
        }
    }
    // else if the vote is still undefined, we think the third path is disagree
    if(vote == VOTE_UNDEFINED && secondLargest >= condition){
        vote = DISAGREE;
    }

    return vote;
}

std::pair<float,float> VariantEdge::findNumberOfRead(int targetPos){
    std::pair<float,float> refBestPair  = ref->BestPair(targetPos);
    std::pair<float,float> altBestPair  = alt->BestPair(targetPos);
    // get the weight of each pair
    float rr = refBestPair.first;
    float ra = refBestPair.second;
    float ar = altBestPair.first;
    float aa = altBestPair.second;
    return std::make_pair( rr + aa , ra +ar );
}

//BlockRead
void BlockRead::recordRead(std::string readName){
    std::map<std::string,int>::iterator readIter = readVec.find(readName);
    if( readIter == readVec.end() )
        readVec[readName] = 1;
    else
        readVec[readName]++;
}

//Handle the special case which One Long Read provides wrong info repeatedly
std::pair<float,float> VairiantGraph::Onelongcase( std::vector<VoteResult> vote ){

    int counter = 0 ;
    float h1 = 0 ;
    float h2 = 0 ;

    // iterate all the voting that previous variants provide
    for ( size_t i = 0 ; i < vote.size() ; i++ ) {

	    // count the votes that refer to only one read
        if ( (vote[i].para+vote[i].cross) <= 1 ) {
            counter++ ;
        }
	    // we will only count the votes that is not INDEL and have lower ESR beacause the INDEL is the variant has higher error rate and the lower ESR means higher reads consistency,
        else if ( vote[i].ESR < 0.2 && vote[i].weight >= 1 && (*variantPosType)[vote[i].Pos].type != INDEL ) {
            if ( vote[i].hap == 1 ) {
                h1+=vote[i].weight ;
            }
            else if ( vote[i].hap == 2 ) {
                h2+=vote[i].weight ;
            }
        }
    }

    //if there has less than three variants use one read to vote we cancel the mechanism
    if ( counter <= 3 || (h1==0&&h2==0) ) {
        return std::make_pair( -1 , -1 ) ;
    }
    else {
        return std::make_pair( h1 , h2 ) ;
    }

}

//VairiantGraph
void VairiantGraph::edgeConnectResult(std::vector<LOHSegment> &LOHSegments){
    // current snp, haplotype (1 or 2), support snp
    // std::map<int, std::map<int,std::vector<int> > > *hpCountMap = new std::map<int, std::map<int,std::vector<int> > >;
    // current variant position, haplotype (1 or 2), previous variants' voting result
    std::map<int, std::map<int,float> > *hpCountMap2 = new std::map<int, std::map<int,float> > ;
    //current variant position, haplotype (1 or 2), previous variants' voting information 
    std::map<int, std::vector<VoteResult> > *hpCountMap3 = new std::map<int, std::vector<VoteResult> > ;
    
    int blockStart = -1;
    int currPos = -1;
    int nextPos = -1;
    int lastConnectPos = -1;
    auto prevIter = variantPosType->end();

    auto lohIter = LOHSegments.begin();
    bool inLOHRegion = false;
    bool genomicEventChange = false;
    int prevPhasedNode = -1;
    int lohStart = -1;
    int lohEnd = -1;
    Haplotype connectHP = HAPLOTYPE1;

    // Visit all position and assign SNPs to haplotype.
    // Avoid recording duplicate information,
    // only one of the two alleles needs to be used for each SNP
    for(auto variantIter = variantPosType->begin() ; variantIter != variantPosType->end() ; variantIter++ ){
        // check next position
        auto nextNodeIter = std::next(variantIter, 1);
        if( nextNodeIter == variantPosType->end() ){
            break;
        }
        
        currPos = variantIter->first;
        nextPos = nextNodeIter->first;
        
        // There should not be a large distance between any two variants, 
        // with the default being a distance of 300000bp, equivalent to one centromere length.
        if(std::abs(nextPos-currPos) > params->distance ){
            continue;
        }
        
        Haplotype currHP = HAPLOTYPE1;
        while (lohIter != LOHSegments.end() && currPos > lohIter->end){
            lohIter++;
            if(lohIter != LOHSegments.end()){
                lohStart = lohIter->start;
                lohEnd = lohIter->end;
            }
        }
        bool backLOHRegion = inLOHRegion;
        inLOHRegion = lohStart <= currPos && currPos <= lohEnd;
        if(backLOHRegion != inLOHRegion){
            genomicEventChange = true;
        }
        if(genomicEventChange && variantIter->second.homozygous == inLOHRegion){
            std::map<int,VariantEdge*>::iterator edgeIter = edgeList->find(prevPhasedNode);
            if(edgeIter!=edgeList->end()){
                float ra = edgeIter->second->ref->getAltReadCount(currPos);
                float ar = edgeIter->second->alt->getRefReadCount(currPos);
                float aa = edgeIter->second->alt->getAltReadCount(currPos);
                float connectRef = inLOHRegion ? ra : ar;
                float connectAlt = aa;
                if(connectRef + connectAlt > 0 && std::max(connectRef, connectAlt) / (connectRef + connectAlt) >= 0.8){
                    bool isRefDominant = (connectRef > connectAlt);
                    if (inLOHRegion) {
                        lohIter->startAllele = isRefDominant ? REF_ALLELE : ALT_ALLELE;
                        if (isRefDominant) {
                            connectHP = ((*posPhasingResult)[prevPhasedNode].refHaplotype == HAPLOTYPE1) ? HAPLOTYPE2 : HAPLOTYPE1;
                        } else {
                            connectHP = ((*posPhasingResult)[prevPhasedNode].refHaplotype == HAPLOTYPE2) ? HAPLOTYPE2 : HAPLOTYPE1;
                        }
                    } else {
                        (lohIter - 1)->endAllele = isRefDominant ? REF_ALLELE : ALT_ALLELE;
                        currHP = isRefDominant ? (connectHP == HAPLOTYPE1 ? HAPLOTYPE2 : HAPLOTYPE1) : connectHP;
                    }
                    hpCountMap2->clear();
                    hpCountMap3->clear();
                    lastConnectPos = currPos;
                }
            }
            genomicEventChange = false;
        }
        if(variantIter->second.homozygous){
            if(inLOHRegion){
                prevPhasedNode = currPos;
            }
            continue;
        }
        // get the number of HP1 and HP2 supported reference allele
        //int h1 = (*hpCountMap)[currPos][HAPLOTYPE1].size();
        //int h2 = (*hpCountMap)[currPos][HAPLOTYPE2].size();
        float h1 = (*hpCountMap2)[currPos][HAPLOTYPE1] ;
        float h2 = (*hpCountMap2)[currPos][HAPLOTYPE2] ;


        //Handle the special case which One Long Read provides wrong info repeatedly
        std::pair<float, float> special = Onelongcase( (*hpCountMap3)[currPos] ) ;
        if ( special.first != -1 ) {
            h1 = special.first ;
            h2 = special.second ;
        }

        // new block, set this position as block start 
        if( h1 == h2 ){
            // No new blocks should be created if the next SNP has already been picked up
            if( currPos < lastConnectPos ){
                if( variantIter->second.origin == SOMATIC ){
                    PhasingResult phasingResult(HAPLOTYPE_UNDEFINED, blockStart, variantIter->second.type, variantIter->second.origin == SOMATIC);
                    posPhasingResult->emplace(currPos, phasingResult);
                }
                continue;
            }
            blockStart = currPos;
        }
        else{
            currHP = ( h1 > h2 ? HAPLOTYPE1 : HAPLOTYPE2 );
        }
        // If 'continue' is not executed, a phasing result is created for the current position
        PhasingResult phasingResult(currHP, blockStart, variantIter->second.type, variantIter->second.origin == SOMATIC);
        posPhasingResult->emplace(currPos, phasingResult);
        if(!inLOHRegion){
            prevPhasedNode = currPos;
        }
        
        // Check if there is no edge from current node
        std::map<int,VariantEdge*>::iterator edgeIter = edgeList->find( currPos );
        if( edgeIter==edgeList->end() ){
            continue;
        }
        // check connect between surrent SNP and next n SNPs
        for(int i = 0 ; i < params->connectAdjacent;){
            if(nextNodeIter->second.homozygous == false){
                i++;
                VoteResult vote(currPos, 1); //used to store previous 20 variants' voting information

                // consider reads from the currnt SNP and the next (i+1)'s SNP
                std::pair<PosAllele,PosAllele> tmp = edgeIter->second->findBestEdgePair(variantIter, nextNodeIter, params->isONT, params->edgeThreshold, vote, false);
                
                // if the target is a danger indel change its weight to 0.1
                if ( (*variantPosType)[currPos].type == DANGER_INDEL ) {
                    vote.weight = 0.1 ;
                }
                // -1 : no connect  
                //  1 : the haplotype of next (i+1)'s SNP are same as previous
                //  2 : the haplotype of next (i+1)'s SNP are different as previous
                if( tmp.first.second != -1 ){
                    // record the haplotype resut of next (i+1)'s SNP
                    if( (*posPhasingResult)[currPos].refHaplotype == HAPLOTYPE1 ){
                        if( tmp.first.second == 1 ){
                            // (*hpCountMap)[nextNodeIter->first][HAPLOTYPE1].push_back(currPos);
                            (*hpCountMap2)[nextNodeIter->first][HAPLOTYPE1] += vote.weight;
                            vote.hap = 1 ;
                        }
                        if( tmp.first.second == 2 ){
                            // (*hpCountMap)[nextNodeIter->first][HAPLOTYPE2].push_back(currPos);
                            (*hpCountMap2)[nextNodeIter->first][HAPLOTYPE2] += vote.weight;
                            vote.hap = 2 ;
                        }
                    }
                    if( (*posPhasingResult)[currPos].refHaplotype == HAPLOTYPE2 ){
                        if( tmp.first.second == 1 ){
                            // (*hpCountMap)[nextNodeIter->first][HAPLOTYPE2].push_back(currPos);
                            (*hpCountMap2)[nextNodeIter->first][HAPLOTYPE2] += vote.weight;
                            vote.hap = 2 ;
                        }
                        if( tmp.first.second == 2 ){
                            // (*hpCountMap)[nextNodeIter->first][HAPLOTYPE1].push_back(currPos);
                            (*hpCountMap2)[nextNodeIter->first][HAPLOTYPE1] += vote.weight;
                            vote.hap = 1 ;
                        }
                    }

                    (*hpCountMap3)[nextNodeIter->first].push_back( vote );

                    if( params->generateDot ){
                        std::string e1 = std::to_string(currPos+1) + ".1\t->\t" + std::to_string(tmp.first.first+1) + "." + std::to_string(tmp.first.second);
                        std::string e2 = std::to_string(currPos+1) + ".2\t->\t" + std::to_string(tmp.second.first+1) + "." + std::to_string(tmp.second.second);

                        dotResult.push_back(e1);
                        dotResult.push_back(e2);
                    }
                    
                    lastConnectPos = nextNodeIter->first;
                }
            }
            nextNodeIter++;
            if( nextNodeIter == variantPosType->end() ){
                break;
            }
        }
        prevIter = variantIter;
    }

    // delete hpCountMap;
    delete hpCountMap2;
    delete hpCountMap3;
}

VairiantGraph::VairiantGraph(std::string &in_ref, PhasingParameters &in_params, std::string &in_chr){
    params=&in_params;
    ref=&in_ref;
    chr = in_chr;
    
    variantPosType = new std::map<int,VariantInfo>;
    edgeList = new std::map<int,VariantEdge*>;
    readHpMap = new std::map<std::string,int>;
}

VairiantGraph::~VairiantGraph(){
}

void VairiantGraph::destroy(){
    dotResult.clear();
    dotResult.shrink_to_fit();
    if(readVariant != nullptr){
        delete readVariant;
        readVariant = nullptr;
    }

    for( auto edgeIter = edgeList->begin() ; edgeIter != edgeList->end() ; edgeIter++ ){
        edgeIter->second->ref->destroy();
        edgeIter->second->alt->destroy();
        edgeIter->second->edgeCount->clear();
        delete edgeIter->second->ref;
        delete edgeIter->second->alt;
        delete edgeIter->second->edgeCount;
        delete edgeIter->second;
    }
    edgeList->clear();
    
    delete variantPosType;
    delete edgeList;
    delete readHpMap;
}
    
void VairiantGraph::addEdge(std::vector<ReadVariant> *in_readVariant){
    readVariant = in_readVariant;
    std::map<std::string,ReadVariant> mergeReadMap;

    // each read will record fist and list variant posistion
    std::map<std::string, std::pair<int,int>> alignRange;
    // record an iterator for all alignments of a read.
    std::map<std::string, std::vector<int>> readIdxVec;
    // record need del read index
    std::vector<int> delReadIdx;

    // Check for overlaps among different alignments of a read and filter out the shorter overlapping alignments.
    for (int readIter = 0; readIter < (int)in_readVariant->size(); readIter++) {
        int is_toDelete = 0;
        std::string readName = (*in_readVariant)[readIter].read_name;
        if((*in_readVariant)[readIter].variantVec.empty())continue;
        int firstVariantPos = (*in_readVariant)[readIter].variantVec.front().position;
        int lastVariantPos = (*in_readVariant)[readIter].variantVec.back().position;
        auto& readRange = alignRange[readName];
        auto& readIdxVecRef = readIdxVec[readName];

        // Initialize readRange if it's the first appearance
        if (alignRange.find(readName) == alignRange.end()) {
            readRange = {firstVariantPos,lastVariantPos};
        } else {
            // Check for overlaps
            while (readRange.first <= firstVariantPos && firstVariantPos <= readRange.second) {
                if (lastVariantPos < readRange.second) {
                    is_toDelete = 1;
                    delReadIdx.push_back(readIter);
                    break;
                }

                int preAlignIdx = readIdxVecRef.size() - 1;
                if (preAlignIdx < 0 ) break;

                const auto& previousAlignment = (*in_readVariant)[readIdxVecRef[preAlignIdx]];
                const auto& prevVariantVec = previousAlignment.variantVec;
                int prevStart = prevVariantVec.front().position;
                int prevEnd = prevVariantVec.back().position;

                double overlapStart = std::max(prevStart, firstVariantPos);
                double overlapEnd = std::min(prevEnd, lastVariantPos);
                if (overlapStart > overlapEnd) break; // No overlap
                double overlapLen = overlapEnd - overlapStart + 1;

                double alignStart = std::max(prevEnd, lastVariantPos);
                double alignEnd = std::min(prevStart, firstVariantPos);
                double alignSpan = alignStart - alignEnd + 1;
                double overlapRatio = overlapLen / alignSpan;

                // Filtering highly overlapping alignments
                if (overlapRatio >= params->overlapThreshold) {
                    int alignLen1 = prevEnd - prevStart + 1;
                    int alignLen2 = lastVariantPos - firstVariantPos + 1;

                    if (alignLen2 <= alignLen1) {
                        is_toDelete = 1;
                        delReadIdx.push_back(readIter); // Current alignment is shorter
                        break;
                    } else {
                        delReadIdx.push_back(readIdxVecRef[preAlignIdx]); // Previous alignment is shorter
                        readIdxVecRef.pop_back();
                        readRange.second = (preAlignIdx > 0) ? (*in_readVariant)[readIdxVecRef[preAlignIdx - 1]].variantVec.back().position : firstVariantPos;
                    }
                } else {
                    break;
                }
            }
            // update range
            readRange.second = lastVariantPos;
        }
        if (is_toDelete == 0 )
            readIdxVecRef.push_back(readIter);
    }

    // sort read index
    std::sort(delReadIdx.begin(), delReadIdx.end());
    // remove overlap alignment
    delReadIdx.push_back((int)in_readVariant->size());
    int saveIter = *(delReadIdx.begin());
    for (auto delIter = delReadIdx.begin(), nextdelIter = std::next(delReadIdx.begin(), 1); nextdelIter != delReadIdx.end(); delIter++ , nextdelIter++) {
        auto nowDelIter = *delIter+1;
        while (nowDelIter<*nextdelIter){
            (*in_readVariant)[saveIter++]=(*in_readVariant)[nowDelIter++];
        }
    }
    in_readVariant->erase( std::next(in_readVariant->begin(), saveIter), in_readVariant->end());

    //position, allele, VariantBases
    std::map<int, std::array<VariantBases, 2>> posAlleleCount;
    for(std::vector<ReadVariant>::iterator readIter = in_readVariant->begin() ; readIter != in_readVariant->end() ; readIter++ ){
        for( auto variant : (*readIter).variantVec ){
            posAlleleCount[variant.position][variant.allele].targetCount++;
            for(auto offsetBaseIter : variant.offsetBase){
                posAlleleCount[variant.position][variant.allele].offsetDiffRefCount[offsetBaseIter.first]++;
            }
        }
    }
    std::ofstream clusterFile(params->resultPrefix + chr + ".cluster");
    std::vector<int> delPos;
    for(const auto& posAlleleCountIter : posAlleleCount) {
        const auto& alleleIter = posAlleleCountIter.second;
        const auto& refOffsetMap = alleleIter[REF_ALLELE].offsetDiffRefCount;
        const auto& altOffsetMap = alleleIter[ALT_ALLELE].offsetDiffRefCount;
        int targetAltCount = alleleIter[ALT_ALLELE].targetCount;

        int sameCount = 0;
        for (const auto& altOffsetMapIter : altOffsetMap) {
            int ra = 0;
            if (refOffsetMap.count(altOffsetMapIter.first)) {
                ra = refOffsetMap.at(altOffsetMapIter.first);
            }
            int aa = altOffsetMapIter.second;
            double condition1 = (double)aa / targetAltCount;
            double condition2 = (double)aa / (ra + aa);
            
            if(condition1 >= 0.5 && condition2 >= 0.6){
                sameCount++;
                // if(sameCount == 2){
                //     break;
                // }
            }
        }
        if(sameCount >= 2){
            delPos.push_back(posAlleleCountIter.first);
            clusterFile << chr << "\t" << posAlleleCountIter.first << "\t" << sameCount << "\n";
        }
    }
    clusterFile.close();

    int readCount=0;
    // merge alignment
    for(std::vector<ReadVariant>::iterator readIter = in_readVariant->begin() ; readIter != in_readVariant->end() ; readIter++ ){

        std::map<std::string,ReadVariant>::iterator posIter = mergeReadMap.find((*readIter).read_name) ;

        // fakeRead is initialize as fake
        if ( posIter == mergeReadMap.end() ) {
            mergeReadMap[(*readIter).read_name].fakeRead = false ;
        }
        //std::cout << (*readIter).mm_rate << "\n";

	    //if the mmrate too high we think it's a fake read
        if( (*readIter).fakeRead == true ){
          mergeReadMap[(*readIter).read_name].fakeRead = true ;
        }

        // Creating a pseudo read which allows filtering out variants that should not be phased
        //ReadVariant tmpRead;
        // Visiting all the variants on the read
        for( auto variant : (*readIter).variantVec ){
            readCount++;
            if(std::binary_search(delPos.begin(), delPos.end(), variant.position)){
                // (*variantPosType)[variant.position].origin = GERMLINE;
                continue;
            }
            if( variant.quality <= VARIANT_UNDEFINED ){
                (*variantPosType)[variant.position].type = static_cast<VariantType>(variant.quality);
            }
            else{
                (*variantPosType)[variant.position].type = SNP;
            }
            (*variantPosType)[variant.position].homozygous = variant.homozygous;

            mergeReadMap[(*readIter).read_name].variantVec.push_back(variant);
        }
    }    
    
    for(std::map<std::string,ReadVariant>::iterator readIter = mergeReadMap.begin() ; readIter != mergeReadMap.end() ; readIter++){
        (*readIter).second.sort();
        
        // iter all pair of snp and construct initial graph
        std::vector<Variant>::iterator variant1Iter = readIter->second.variantVec.begin();
        std::vector<Variant>::iterator variant2Iter = std::next(variant1Iter,1);
        std::vector<Variant>::iterator variant3Iter = std::next(variant2Iter,1);
        
        while(variant1Iter != readIter->second.variantVec.end() && variant2Iter != readIter->second.variantVec.end() ){
            // create new edge if not exist
            std::map<int,VariantEdge*>::iterator posIter = edgeList->find(variant1Iter->position);
            if( posIter == edgeList->end() )
                (*edgeList)[variant1Iter->position] = new VariantEdge(variant1Iter->position);

            //count the ref and alt base amount on the variant
            if( (*variant1Iter).allele == 0 && (*readIter).second.fakeRead == false ) {
                (*edgeList)[(*variant1Iter).position]->refcnt++ ;
            }
            if( (*variant1Iter).allele == 1 && (*readIter).second.fakeRead == false ) {
                (*edgeList)[(*variant1Iter).position]->altcnt++ ;
            }
            (*edgeList)[(*variant1Iter).position]->coverage++;
            auto node = (*edgeList)[(*variant1Iter).position];
            int addHomCount = 0;
            // add edge process
            for(int nextNode = 0 ; nextNode < params->connectAdjacent;){
                if(variant1Iter->homozygous == true || variant2Iter->homozygous == true){
                    addHomCount++;
                }
                if(variant1Iter->homozygous == false || variant2Iter->homozygous == false){
                    nextNode++;
                    for(int nextNextNode = 0 ; nextNextNode < params->somaticConnectAdjacent && nextNode < params->somaticConnectAdjacent;){
                        if( variant3Iter == readIter->second.variantVec.end() ){
                            break;
                        }
                        if(variant3Iter->homozygous == false){
                            nextNextNode++;
                            node->addSegmentEdge((*variant1Iter), (*variant2Iter), (*variant3Iter));
                        }
                        variant3Iter++;
                    }
                }
                if(addHomCount <= 6){
                    // this allele support ref
                    if( variant1Iter->allele == 0 )
                        node->ref->addSubEdge((*variant1Iter), (*variant2Iter), readIter->first, params->baseQuality, params->edgeWeight, (*readIter).second.fakeRead);
                    // this allele support alt
                    if( (*variant1Iter).allele == 1 )
                        node->alt->addSubEdge((*variant1Iter), (*variant2Iter), readIter->first, params->baseQuality, params->edgeWeight, (*readIter).second.fakeRead);
                }
                // next snp
                variant2Iter++;
                variant3Iter = std::next(variant2Iter,1);
                if( variant2Iter == readIter->second.variantVec.end() ){
                    break;
                }
            }

            variant1Iter++;
            variant2Iter = std::next(variant1Iter,1);
            variant3Iter = std::next(variant2Iter,1);
        }

        //count the ref and alt base amount of the last variant on the read
	    if ( variant1Iter != (*readIter).second.variantVec.end() && variant2Iter == (*readIter).second.variantVec.end() ) {
            std::map<int,VariantEdge*>::iterator posIter = edgeList->find((*variant1Iter).position);
            if( posIter == edgeList->end() ) {
                (*edgeList)[(*variant1Iter).position] = new VariantEdge((*variant1Iter).position);
                //(*edgeList)[(*variant1Iter).position]->vaf = (*currentVariants)[(*variant1Iter).position].vaf ;
            }

            if( (*variant1Iter).allele == 0 && (*readIter).second.fakeRead == false)
                (*edgeList)[(*variant1Iter).position]->refcnt++ ;
            if( (*variant1Iter).allele == 1 && (*readIter).second.fakeRead == false)
                (*edgeList)[(*variant1Iter).position]->altcnt++ ;
            (*edgeList)[(*variant1Iter).position]->coverage++;
        }
    }
}

void VairiantGraph::somaticCalling(std::map<int, RefAlt>* variants){
    std::string logFilePath = params->resultPrefix + chr + ".somatic";
    std::ofstream logFile(logFilePath, std::ofstream::out | std::ofstream::trunc);
    if (!logFile.is_open()) {
        std::cerr << "failed to open log file: " << logFilePath << std::endl;
    }
    std::map<int, std::array<int, 18>> voteResult;
    std::map<int, double> totalArtifactPathRatio;
    std::map<int, VariantInfo>::iterator nodeIter;
    std::map<int, VariantInfo>::iterator nextNodeIter;
    std::map<int, VariantInfo>::iterator nextNextNodeIter;

    for(nodeIter = variantPosType->begin() ; nodeIter != variantPosType->end() ; nodeIter++ ){
        auto variantIter = variants->find(nodeIter->first);
        if(variantIter->second.germline){
            nodeIter->second.origin = GERMLINE;
        }
        if (nodeIter->second.homozygous) {
            continue;
        }
        // check next position
        nextNodeIter = std::next(nodeIter, 1);
        if( nextNodeIter == variantPosType->end() ){
            break;
        }

        // There should not be a large distance between any two variants, 
        // with the default being a distance of 300000bp, equivalent to one centromere length.
        if(std::abs(nextNodeIter->first - nodeIter->first) > params->distance ){
            continue;
        }

        // Check if there is no edge from current node
        std::map<int,VariantEdge*>::iterator edgeIter = edgeList->find(nodeIter->first);
        if( edgeIter == edgeList->end() ){
            continue;
        }
        VariantEdge *edge = edgeIter->second;

        // check connect between surrent SNP and next n SNPs
        for(int i = 0 ; i < params->somaticConnectAdjacent;){
            if( nextNodeIter == variantPosType->end() ){
                break;
            }
            if (!nextNodeIter->second.homozygous) {
                i++;
                nextNextNodeIter = std::next(nextNodeIter, 1);
                for(int j = 0 ; j < params->somaticConnectAdjacent;){
                    if( nextNextNodeIter == variantPosType->end() ){
                        break;
                    }
                    if (!nextNextNodeIter->second.homozygous) {
                        j++;

                        ThreePointEdge threePointEdge = edge->findThreePointEdge(nextNodeIter->first,nextNextNodeIter->first);
                        int voteTmp = patternMining(threePointEdge);
                        
                        if (voteTmp == LEFT_HIGH_SA_PAR || voteTmp == LEFT_HIGH_SR_PAR || voteTmp == LEFT_HIGH_SA_CR || voteTmp == LEFT_HIGH_SR_CR || 
                            voteTmp == LEFT_LOW_SA_PAR){
                            voteResult[nodeIter->first][voteTmp]++;
                            voteResult[nextNodeIter->first][DISAGREE]++;
                            voteResult[nextNextNodeIter->first][DISAGREE]++;
                        } else if (voteTmp == RIGHT_HIGH_SA_PAR || voteTmp == RIGHT_HIGH_SR_PAR || voteTmp == RIGHT_HIGH_SA_CR || voteTmp == RIGHT_HIGH_SR_CR || 
                                voteTmp == RIGHT_LOW_SA_PAR){
                            voteResult[nextNextNodeIter->first][voteTmp]++;
                            voteResult[nextNodeIter->first][DISAGREE]++;
                            voteResult[nodeIter->first][DISAGREE]++;
                        } else if (voteTmp == MID_HIGH_SA_PAR || voteTmp == MID_HIGH_SR_PAR || voteTmp == MID_HIGH_SA_CR || voteTmp == MID_HIGH_SR_CR || 
                                   voteTmp == MID_LOW_SA_PAR || voteTmp == DISAGREE){
                            voteResult[nextNodeIter->first][voteTmp]++;
                            voteResult[nextNextNodeIter->first][DISAGREE]++;
                            voteResult[nodeIter->first][DISAGREE]++;
                        }

                        EdgeType edgeType = EDGE_RRR;
                        if (voteTmp == LEFT_HIGH_SA_PAR || voteTmp == LEFT_HIGH_SR_PAR || voteTmp == LEFT_HIGH_SA_CR || voteTmp == LEFT_HIGH_SR_CR){
                            if (voteTmp == LEFT_HIGH_SA_PAR){
                                edgeType = EDGE_AAA;
                            } else if (voteTmp == LEFT_HIGH_SR_PAR){
                                edgeType = EDGE_ARR;
                            } else if (voteTmp == LEFT_HIGH_SA_CR){
                                edgeType = EDGE_AAR;
                            } else if (voteTmp == LEFT_HIGH_SR_CR){
                                edgeType = EDGE_ARA;
                            }
                            totalArtifactPathRatio[nodeIter->first] += (threePointEdge[edgeType] / (threePointEdge[EDGE_ARR]+threePointEdge[EDGE_ARA]+threePointEdge[EDGE_AAR]+threePointEdge[EDGE_AAA]));
                        } else if (voteTmp == RIGHT_HIGH_SA_PAR || voteTmp == RIGHT_HIGH_SR_PAR || voteTmp == RIGHT_HIGH_SA_CR || voteTmp == RIGHT_HIGH_SR_CR){
                            if (voteTmp == RIGHT_HIGH_SA_PAR){
                                edgeType = EDGE_RRA;
                            } else if (voteTmp == RIGHT_HIGH_SR_PAR){
                                edgeType = EDGE_AAA;
                            } else if (voteTmp == RIGHT_HIGH_SA_CR){
                                edgeType = EDGE_ARA;
                            } else if (voteTmp == RIGHT_HIGH_SR_CR){
                                edgeType = EDGE_RAA;
                            }
                            totalArtifactPathRatio[nextNextNodeIter->first] += (threePointEdge[edgeType] / (threePointEdge[EDGE_RRA]+threePointEdge[EDGE_RAA]+threePointEdge[EDGE_ARA]+threePointEdge[EDGE_AAA]));
                        } else if (voteTmp == MID_HIGH_SA_PAR || voteTmp == MID_HIGH_SR_PAR || voteTmp == MID_HIGH_SA_CR || voteTmp == MID_HIGH_SR_CR){
                            if (voteTmp == MID_HIGH_SA_PAR){
                                edgeType = EDGE_AAA;
                            } else if (voteTmp == MID_HIGH_SR_PAR){
                                edgeType = EDGE_RAR;
                            } else if (voteTmp == MID_HIGH_SA_CR){
                                edgeType = EDGE_AAR;
                            } else if (voteTmp == MID_HIGH_SR_CR){
                                edgeType = EDGE_RAA;
                            }
                            totalArtifactPathRatio[nextNodeIter->first] += (threePointEdge[edgeType] / (threePointEdge[EDGE_RAR]+threePointEdge[EDGE_RAA]+threePointEdge[EDGE_AAR]+threePointEdge[EDGE_AAA]));
                        }
                    }
                    nextNextNodeIter++;
                }
            }
            nextNodeIter++;
        }
        if(nodeIter->second.origin == GERMLINE){
            continue;
        }
        auto voteResultIter = voteResult.find(nodeIter->first);
        if(voteResultIter != voteResult.end()){
            const auto& voteResultArray = voteResultIter->second;
            int highLeftVote = voteResultArray[LEFT_HIGH_SR_CR] + voteResultArray[LEFT_HIGH_SA_PAR] + voteResultArray[LEFT_HIGH_SR_PAR] + voteResultArray[LEFT_HIGH_SA_CR];
            int highRightVote = voteResultArray[RIGHT_HIGH_SR_PAR] + voteResultArray[RIGHT_HIGH_SA_CR] + voteResultArray[RIGHT_HIGH_SR_CR] + voteResultArray[RIGHT_HIGH_SA_PAR];
            int highMidVote = voteResultArray[MID_HIGH_SR_CR] + voteResultArray[MID_HIGH_SA_PAR] + voteResultArray[MID_HIGH_SR_PAR] + voteResultArray[MID_HIGH_SA_CR];
            int highAllVote = highLeftVote + highRightVote + highMidVote;
            int lowVote = voteResultArray[MID_LOW_SA_PAR]+voteResultArray[LEFT_LOW_SA_PAR]+voteResultArray[RIGHT_LOW_SA_PAR];
            float allLowVote = voteResultArray[DISAGREE] + lowVote;
            double totalRaito = totalArtifactPathRatio[nodeIter->first];
            logFile<< chr << "\t" << nodeIter->first << "\t" << highAllVote << "\t" << lowVote << "\t" << voteResultArray[DISAGREE] << "\t" << totalRaito << "\n";
            if((highAllVote > 0 && (totalArtifactPathRatio[nodeIter->first]/highAllVote) > 0.8) || (lowVote > 0 && (lowVote/allLowVote) >= 0.2)){
                nodeIter->second.origin = SOMATIC;
            }
        }
    }
    logFile.close();
}

void VairiantGraph::readCorrection(std::map<double, int> *ploidyRatioMap){
    // Allocate a map for storing haplotype allele counts
    // haplotype, <position <allele, base count>>
    std::map<int,std::map<int,std::map<double,double>>> *hpAlleleCountMap = new std::map<int,std::map<int,std::map<double,double>>>;

    // Process reads to compute haplotype allele counts
    processReadVariants(hpAlleleCountMap);

    // Reassign allele results based on the computed haplotype counts
    reassignAlleleResult(hpAlleleCountMap, ploidyRatioMap);

    delete hpAlleleCountMap;
}

void VairiantGraph::processReadVariants(std::map<int,std::map<int,std::map<double,double>>> *hpAlleleCountMap) {
    // phasing result and variant allele mapping
    const int variantHaplotype[2][2] = {
        {HAPLOTYPE1, HAPLOTYPE2},
        {HAPLOTYPE2, HAPLOTYPE1}
    };

    // iter all read, determine the haplotype of the read
    for(std::vector<ReadVariant>::iterator readIter = (*readVariant).begin() ; readIter != (*readVariant).end() ; readIter++ ){
        double haplotype1Count = 0;
        double haplotype2Count = 0;
        //int block;
        
        // loop all variant 
        for( auto variant : (*readIter).variantVec ){
            const auto& posPhasingResultIter = posPhasingResult->find(variant.position);
            if( posPhasingResultIter != posPhasingResult->end() && posPhasingResultIter->second.refHaplotype != HAPLOTYPE_UNDEFINED ){
                const Haplotype refHaplotype = posPhasingResultIter->second.refHaplotype;
                std::map<int,VariantEdge*>::iterator edgeIter = edgeList->find( variant.position );
                //when vaf is 0 or 1, the fakeSnp will be true
                bool fakeSnp = edgeIter->second->get_fakeSnp();
                double edgeWeight = 1;
                if(fakeSnp){
                    edgeWeight = 0.01;
                }
                else if (variant.quality == DANGER_INDEL || variant.quality == INDEL) {
                    edgeWeight = 0.1;
                }else if (variant.quality == MOD_FORWARD_STRAND || variant.quality == MOD_REVERSE_STRAND) {
                    continue;
                }
                // If variant is not somatic or is somatic with ALT allele, add edgeWeight to the corresponding haplotype count
                if (!posPhasingResultIter->second.somatic || variant.allele == ALT_ALLELE) {
                    double &hapCount = (variantHaplotype[refHaplotype][variant.allele] == HAPLOTYPE1) ? haplotype1Count : haplotype2Count;
                    hapCount += edgeWeight;
                }
            }
            
        }
        
        // tag high confident reads
        if( std::max(haplotype1Count,haplotype2Count)/(haplotype1Count+haplotype2Count) > params->readConfidence && (haplotype1Count + haplotype2Count) > 1 ){
            // tag read with the corresponding haplotype
            int belongHP = ( haplotype1Count > haplotype2Count ? 0 : 1 );
            (*readHpMap)[(*readIter).read_name] = belongHP;
            
            //readBlockHP[(*readIter).read_name][(*readIter).reference_start][block]=belongHP;
            //readBlockHPcount[(*readIter).read_name][block][belongHP]++;
            
            for(auto variantIter = (*readIter).variantVec.begin() ; variantIter != (*readIter).variantVec.end() ; variantIter++ ){
                if( (*variantIter).allele == 0 || (*variantIter).allele == 1){
                    // when the mmrate is too high, we think it's a fakeRead
                    // (*hpAlleleCountMap)[belongHP][(*variantIter).position][(*variantIter).allele]++;
                    if( (*readIter).fakeRead == true )
                      (*hpAlleleCountMap)[belongHP][(*variantIter).position][(*variantIter).allele]+=0.01;
                    else
                      (*hpAlleleCountMap)[belongHP][(*variantIter).position][(*variantIter).allele]++;
                }
            }
        }
        else{
            (*readHpMap)[(*readIter).read_name] = -1;
        }
    }
}

void VairiantGraph::calculatePloidyRatio(double hp1Ref, double hp2Ref, std::map<double, int> *ploidyRatioMap) {
    if (hp1Ref + hp2Ref <= 0) return;
    double ploidyRatio = std::max(hp1Ref, hp2Ref) / (hp1Ref + hp2Ref);
    ploidyRatioMap->emplace(ploidyRatio, 0).first->second++;
}

void VairiantGraph::reassignAlleleResult(std::map<int,std::map<int,std::map<double,double>>> *hpAlleleCountMap, std::map<double, int> *ploidyRatioMap) {
    double snpConfidenceThreshold = params->snpConfidence;
    std::map<int,std::map<int,int>> hpAllele;
    std::string logFilePath = params->resultPrefix + chr + ".purity";
    std::ofstream logFile(logFilePath, std::ofstream::out | std::ofstream::trunc);
    if (!logFile.is_open()) {
        std::cerr << "failed to open log file: " << logFilePath << std::endl;
    }
    // reassign allele result
    for(auto variantIter = variantPosType->begin() ; variantIter != variantPosType->end() ; variantIter++ ){
        int position = variantIter->first;
        auto posPhasingResultIter = posPhasingResult->find(position);
        if (posPhasingResultIter == posPhasingResult->end() || posPhasingResultIter->second.refHaplotype == HAPLOTYPE_UNDEFINED) {
            continue;
        }
        
        double hp1Ref = (*hpAlleleCountMap)[0][position][0];
        double hp1Alt = (*hpAlleleCountMap)[0][position][1];
        double hp2Ref = (*hpAlleleCountMap)[1][position][0];
        double hp2Alt = (*hpAlleleCountMap)[1][position][1];
        if(ploidyRatioMap != nullptr && posPhasingResultIter->second.somatic){
            logFile << chr << "\t" << position << "\t" << 1 << "\t" << hp1Ref << "\t" << hp1Alt << "\t" << hp2Ref << "\t" << hp2Alt << "\t" << 1 << std::endl;

            calculatePloidyRatio(hp1Ref, hp2Ref, ploidyRatioMap);
        }
        double result1reads = hp1Ref + hp2Alt;
        double result2reads = hp2Ref + hp1Alt;
        if(posPhasingResultIter->second.somatic){
            result1reads = hp2Alt;
            result2reads = hp1Alt;
        }
        double resultConfidence = std::max(result1reads, result2reads) / (result1reads + result2reads);
        
        Haplotype refHaplotypeResult = HAPLOTYPE_UNDEFINED;
        Haplotype altHaplotypeResult = HAPLOTYPE_UNDEFINED;
        //std::cout << "RC\t" << position+1 << "\t" << result1reads << "\t" << result2reads << "\t" << resultConfidence << "\n" ;
        
        if( resultConfidence > snpConfidenceThreshold ){
            if( result1reads > result2reads ){
                refHaplotypeResult = HAPLOTYPE1;
                altHaplotypeResult = HAPLOTYPE2;
            }
            else if( result1reads < result2reads ){
                refHaplotypeResult = HAPLOTYPE2;
                altHaplotypeResult = HAPLOTYPE1;
            }
        }

        if( refHaplotypeResult != HAPLOTYPE_UNDEFINED && altHaplotypeResult != HAPLOTYPE_UNDEFINED ){
            posPhasingResultIter->second.refHaplotype = refHaplotypeResult;
        }
        else{
            if(!posPhasingResultIter->second.somatic){
                posPhasingResult->erase(posPhasingResultIter);
            }else{
                posPhasingResultIter->second.refHaplotype = HAPLOTYPE_UNDEFINED;
            }
        }
    }
    logFile.close();
}

void VairiantGraph::convertNonGermlineToSomatic() {
    for(auto variantIter = variantPosType->begin(); variantIter != variantPosType->end(); variantIter++) {
        if(variantIter->second.origin != GERMLINE) {
            variantIter->second.origin = SOMATIC;
        }
    }
}

void VairiantGraph::exportPhasingResult(PosPhasingResult &posPhasingResult, std::vector<LOHSegment> &LOHSegments) {
    bool genomicEventChange = false;
    int lastPhaseSet = -1;
    int assignPhaseSet = -1;
    Haplotype lastHP = HAPLOTYPE_UNDEFINED;
    Haplotype connectedHP = HAPLOTYPE_UNDEFINED;
    Allele nextConnectedAllele = Allele_UNDEFINED;
    auto lohIter = LOHSegments.begin();
    bool inLOHRegion = false;
    for(auto variantIter = variantPosType->begin() ; variantIter != variantPosType->end() ; variantIter++ ){
        while(lohIter != LOHSegments.end() && variantIter->first > lohIter->end){
            lohIter++;
            connectedHP = HAPLOTYPE_UNDEFINED;
            nextConnectedAllele = Allele_UNDEFINED;
            assignPhaseSet = -1;
        }
        bool backLOHRegion = inLOHRegion;
        inLOHRegion = (lohIter != LOHSegments.end() && variantIter->first >= lohIter->start && variantIter->first <= lohIter->end);
        if(backLOHRegion != inLOHRegion){
            genomicEventChange = true;
        }
        if(genomicEventChange && inLOHRegion){
            genomicEventChange = false;
            Allele connectedAllele = lohIter->startAllele;
            if(connectedAllele != Allele_UNDEFINED){
                connectedHP = (connectedAllele == REF_ALLELE) ? 
                            (lastHP == HAPLOTYPE1 ? HAPLOTYPE2 : HAPLOTYPE1) :
                            lastHP;
                // lastPhaseSet = lastPhaseSet;
            }else{
                lastPhaseSet = variantIter->first;
            }
            nextConnectedAllele = lohIter->endAllele;
        }
        if(variantIter->second.homozygous){
            if(inLOHRegion){
                std::string genotype = (connectedHP == HAPLOTYPE1) ? ".|1" : "1|.";
                PhasingResult phasingResult(lastPhaseSet, variantIter->second.type, genotype);
                if(variantIter->second.origin == SOMATIC){
                    std::replace(genotype.begin(), genotype.end(), '1', '0');
                    phasingResult.genotype.insert(phasingResult.genotype.begin(), genotype);
                    phasingResult.phaseSet.push_back(-1);
                }
                posPhasingResult.emplace(variantIter->first, phasingResult);
            }
        }else{
            auto posPhasingResultIter = posPhasingResult.find(variantIter->first);
            if(posPhasingResultIter != posPhasingResult.end()){
                auto &result = posPhasingResultIter->second;
                if(inLOHRegion){
                    if(result.somatic){
                        if(connectedHP == HAPLOTYPE1){
                            if(result.refHaplotype == HAPLOTYPE_UNDEFINED){
                                result.genotype = {".|0", ".|1", ".|1"};
                            }else if(result.refHaplotype == HAPLOTYPE1){
                                result.genotype = {".|0", ".|0", ".|1"};
                            }else if(result.refHaplotype == HAPLOTYPE2){
                                result.genotype = {".|0", ".|1", ".|0"};
                            }
                        }else{
                            if(result.refHaplotype == HAPLOTYPE_UNDEFINED){
                                result.genotype = {"0|.", "1|.", "1|."};
                            }else if(result.refHaplotype == HAPLOTYPE1){
                                result.genotype = {"0|.", "0|.", "1|."};
                            }else if(result.refHaplotype == HAPLOTYPE2){
                                result.genotype = {"0|.", "1|.", "0|."};
                            }
                        }
                        int phaseSet = result.phaseSet.front();
                        result.phaseSet={lastPhaseSet, phaseSet, phaseSet};
                    }else{
                        posPhasingResult.erase(posPhasingResultIter);
                    }
                }else{
                    if(genomicEventChange){
                        genomicEventChange = false;
                        if(nextConnectedAllele != Allele_UNDEFINED){
                            assignPhaseSet = result.phaseSet.front();
                        }
                    }
                    if(result.somatic){
                        if(result.refHaplotype == HAPLOTYPE_UNDEFINED){
                            result.genotype = {"0|0", "1|1"};
                        }
                        else if(result.refHaplotype == HAPLOTYPE1){
                            result.genotype = {"0|0", ".|1"};
                        }
                        else if(result.refHaplotype == HAPLOTYPE2){
                            result.genotype = {"0|0", "1|."};
                        }
                        result.phaseSet.push_back(-1);
                    }else{
                        if(result.refHaplotype == HAPLOTYPE1){
                            result.genotype = {"0|1"};
                        }
                        else if(result.refHaplotype == HAPLOTYPE2){
                            result.genotype = {"1|0"};
                        }
                    }
                    if(assignPhaseSet == result.phaseSet.front()){
                        result.phaseSet.front() = lastPhaseSet;
                    }
                    lastHP = result.refHaplotype;
                    lastPhaseSet = result.phaseSet.front();
                }
            }
        }
    }
}

void VairiantGraph::writingDotFile(std::string dotPrefix){
    
    std::ofstream resultVcf(dotPrefix+".dot");

    if(!resultVcf.is_open()){
        std::cerr<< "Fail to open write file: " << dotPrefix+".vcf" << "\n";
    }
    else{
        resultVcf << "digraph G {\n";

        for(auto edge : dotResult){
            resultVcf << edge << "\n";
        }
        resultVcf << "}\n";
    }
    return;
}

std::map<std::string,int>* VairiantGraph::getReadHP(){
    return readHpMap;
}

int VairiantGraph::totalNode(){
    return variantPosType->size();
}

void VairiantGraph::phasingProcess(PosPhasingResult &inPosPhasingResult, std::vector<LOHSegment> &LOHSegments, std::map<double, int> *ploidyRatioMap){
    posPhasingResult = &inPosPhasingResult;

    // This step involves converting all reads into a graph structure, which will be stored as an edge list
    // in a two-layer map. The first layer of the map uses the starting coordinate as the key and contains
    // a second layer map as the value. The second layer map uses the destination coordinate as the key and
    // stores the number of support read as values. (There is another map used for debugging purposes that
    // treats the read name vector as a value.) The method begins by visiting the coordinates covered by each 
    // read and recording this information in 'variantPosType.' Subsequently, it connects the coordinates contained 
    // in each read on the graph. Specifically, each coordinate is connected to the next N coordinates in a 
    // linear fashion.
    this->edgeConnectResult(LOHSegments);

    // This step will utilize the results of graph phasing to attempt to separate all the reads into two 
    // haplotypes and then identify high-confidence SNPs using reads from the two distinct haplotypes.
    this->readCorrection(ploidyRatioMap);  
}

Roller::Roller() {}

Roller::Roller(size_t windowSize) : window(windowSize) {}

Roller::Roller(size_t windowSize, const std::vector<double>& input) : window(windowSize), values(input) {}

Roller& Roller::setValues(const std::vector<double>& input) {
    values = input;
    return *this;
}
    
Roller& Roller::smooth() {
    values = smooth(values);
    return *this;
}
    
Roller& Roller::forward() {
    values = forward(values);
    return *this;
}
    
Roller& Roller::reverse() {
    values = reverse(values);
    return *this;
}
    
Roller& Roller::downSampling(int distance) {
    values = downSampling(values, distance);
    return *this;
}
    
Roller& Roller::directionalDifference(size_t distance) {
    values = directionalDifference(values, distance);
    return *this;
}
    
Roller& Roller::opposite() {
    values = opposite(values);
    return *this;
}
    
Roller& Roller::reverseNetValues() {
    values = reverseNetValues(values);
    return *this;
}

// Get the final result
std::vector<double> Roller::getResult() const {
    return values;
}

std::vector<double> Roller::smooth(const std::vector<double>& values) {
    std::vector<double> results(values.size());
    std::deque<double> windowValues;
    double rollingSum = 0.0;
    size_t halfWindow = window / 2;

    // initialize the first window
    for (size_t i = 0; i <= std::min(halfWindow, values.size() - 1); ++i) {
        windowValues.push_back(values[i]);
        rollingSum += values[i];
    }
    results[0] = (rollingSum / windowValues.size());
    
    // start to process the rest of the values
    for (size_t i = 1; i < values.size(); ++i) {
        updateWindow(rollingSum, windowValues, values[i]);
        results[i] = (rollingSum / windowValues.size());
    }

    return results;
}

std::vector<double> Roller::forward(const std::vector<double>& values) {
    std::vector<double> results(values.size());
    std::deque<double> windowValues;
    double rollingSum = 0;

    for (int i = values.size() - 1; i >= 0; --i) {
        updateWindow(rollingSum, windowValues, values[i]);
        results[i] = (double)rollingSum / window;
    }
    return results;
}

std::vector<double> Roller::reverse(const std::vector<double>& values) {
    std::vector<double> results(values.size());
    std::deque<double> windowValues;
    double rollingSum = 0;

    for (size_t i = 0; i < values.size(); ++i) {
        updateWindow(rollingSum, windowValues, values[i]);
        results[i] = (double)rollingSum / window;
    }
    return results;
}

std::vector<double> Roller::downSampling(const std::vector<double>& values, int distance) {
    std::vector<double> results(values.size() / distance + 1);
    for (size_t i = 0; i < values.size(); i+=distance) {
        results[i / distance] = values[i];
    }
    return results;
}

std::vector<double> Roller::directionalDifference(const std::vector<double>& values, size_t distance) {
    std::vector<double> results(values.size());

    for (size_t i = 0; i < values.size(); ++i) {
        double forwardShift = (i + distance < values.size()) ? values[i + distance] : 0.0;
        double backwardShift = (i >= distance) ? values[i - distance] : 0.0;
        results[i] = forwardShift - backwardShift;
    }
    return results;
}

std::vector<double> Roller::opposite(const std::vector<double>& values) {
    std::vector<double> results(values.size());
    for (size_t i = 0; i < values.size(); ++i) {
        results[i] = -values[i];
    }
    return results;
}

std::vector<double> Roller::reverseNetValues(const std::vector<double>& values) {
    std::vector<double> results(values.size());
    results[0] = values[0];
    for (size_t i = 1; i < values.size(); ++i) {
        results[i] = values[i] - values[i-1];
    }
    return results;
}

std::vector<double> Roller::backValues(const std::vector<size_t>& idxs, const std::vector<double>& values) {
    std::vector<double> results(idxs.size());
    for (size_t i = 0; i < idxs.size(); ++i) {
        results[i] = values[idxs[i]];
    }
    return results;
}

std::vector<int> Roller::backPosition(const std::vector<size_t>& idxs, const std::vector<int>& values) {
    std::vector<int> results(idxs.size());
    for (size_t i = 0; i < idxs.size(); ++i) {
        results[i] = values[idxs[i]];
    }
    return results;
}

void Roller::replaceValue(ClipCount &clipCount, std::vector<int>& keys, std::vector<size_t>& idxs, std::vector<double>& replaceValues) {
    for (size_t i = 0; i < idxs.size(); ++i) {
        if(replaceValues[i] > 0){
            clipCount[keys[idxs[i]]][FRONT] = replaceValues[i];
        }
        else{
            clipCount[keys[idxs[i]]][BACK] = -replaceValues[i];
        }
    }
}

void Roller::updateWindow(double& rollingSum, std::deque<double>& windowValues, double newValue) {
    // add the new value
    rollingSum += newValue;
    windowValues.push_back(newValue);
    // remove the old value
    if (windowValues.size() > window) {
        rollingSum -= windowValues.front();
        windowValues.pop_front();
    }
}

// PeakFinder Implementation
PointFinder::PointFinder(size_t peakDistance, int calibrationThreshold, size_t calibrationDistance)
    : peakDistance(peakDistance), calibrationThreshold(calibrationThreshold) {
    this->calibrationDistance = (calibrationDistance == 0) ? peakDistance * 2 : calibrationDistance;
}

std::vector<size_t> PointFinder::findAllPeaks(const std::vector<double>& values, double peakThreshold) {
    std::vector<size_t> peaks;
    bool isPeak = true;
    for (size_t i = 1; i < values.size() - 1; ++i) {
        isPeak = true;
        if (values[i] >= peakThreshold) {
            size_t start = this->start(i, peakDistance, 2);
            size_t end = this->end(i, peakDistance, values.size() - 1, 2);
            for (size_t j = start; j <= end; ++j) {
                if (values[j] > values[i]){
                    isPeak = false;
                    break;
                }
            }
            if (isPeak){
                peaks.push_back(i);
            }
        }
    }
    return peaks;
}

std::vector<size_t> PointFinder::findPeaks(const std::vector<double>& values, double peakThreshold) {
    std::vector<size_t> peakKeys;
    if (values.empty()){
        return peakKeys;
    }
    size_t valuesEnd = values.size() - 1;

    for (size_t i = 1; i < valuesEnd; ++i) {
        if (values[i] >= peakThreshold) {
            size_t start = this->start(i, peakDistance, 2);
            size_t end = this->end(i, peakDistance, valuesEnd, 2);
            double start_value = start == 0 ? 0 : values[start];
            double end_value = end == valuesEnd ? 0 : values[end];
            if (values[i] >= start_value && values[i] >= end_value) {
                if (!peakKeys.empty() && i - peakKeys.back() < peakDistance){
                    if (values[i] > values[peakKeys.back()]){
                        peakKeys.back() = i;
                    }
                }
                else{
                    peakKeys.push_back(i);
                }
            }
        }
    }
    return peakKeys;
}

std::vector<size_t> PointFinder::getForwardGentle(const std::vector<double>& values, std::vector<size_t> peaks) {
    std::vector<size_t> gentlePeaks;
    for (size_t& peak : peaks) {
        bool isGentle = true;
        size_t end = this->end(peak, calibrationDistance, values.size());
        for (size_t j = peak; j < end; ++j) {
            if (values[j + 1] >= calibrationThreshold) {
                isGentle = false;
                break;
            }
        }
        if (isGentle){
            gentlePeaks.push_back(peak);
        }
    }
    return gentlePeaks;
}

std::vector<size_t> PointFinder::getReverseGentle(const std::vector<double>& values, std::vector<size_t> peaks) {
    std::vector<size_t> gentlePeaks;
    for (size_t& peak : peaks) {
        bool isGentle = true;
        size_t start = this->start(peak, calibrationDistance);
        for (size_t j = peak; j > start; --j) {
            if (values[j - 1] <= -calibrationThreshold) {
                isGentle = false;
                break;
            }
        }
        if (isGentle){
            gentlePeaks.push_back(peak);
        }
    }
    return gentlePeaks;
}

size_t PointFinder::start(size_t target, size_t distance, size_t distanceDivisor) {
    return (target > distance) ? target - (distance / distanceDivisor) : 0;
}

size_t PointFinder::end(size_t target, size_t distance, size_t size, size_t distanceDivisor) {
    return std::min(target + (distance / distanceDivisor), size);
}

Clip::Clip(std::string &chr){
    this->chr = chr;
    this->largeGenomicEventInterval = nullptr;
    this->smallGenomicEventRegion = nullptr;
}

Clip::~Clip(){
}

void Clip::detectGenomicEventInterval(ClipCount &clipCount, std::vector<int> &largeGenomicEventInterval, std::vector<std::pair<int, int>> &smallGenomicEventRegion){
    this->largeGenomicEventInterval = &largeGenomicEventInterval;
    this->smallGenomicEventRegion = &smallGenomicEventRegion;
    // Check if there are any clip counts to process
    if (clipCount.size() > 0){
        // Aggregate gentle clip values into a consolidated point
        this->amplifyGentleClip(clipCount);
        // Calculate and store the genomic event intervals, including both large and small event regions
        this->detectInterval(clipCount);
    }
}

void Clip::amplifyGentleClip(ClipCount &clipCount){
    size_t window = 100;
    size_t clipCountSize = clipCount.size();
    if (clipCountSize <= window) return;
    std::vector<int> positions(clipCountSize);
    std::vector<double> netValues(clipCountSize);
    std::vector<double> cumulativeSum(clipCountSize);

    int totalNetValue = 0;
    size_t i = 0;
    for (auto& clipIter : clipCount) {
        auto& counts = clipIter.second;
        int netValue = counts[FRONT] - counts[BACK];
        totalNetValue += netValue;
        positions[i] = clipIter.first;
        netValues[i] = netValue;
        cumulativeSum[i] = totalNetValue;
        i++;
    }

    std::vector<double> directionalDifference = Roller(window, cumulativeSum).smooth().directionalDifference().getResult();
    std::vector<double> smoothedForwardSum = Roller(window, netValues).forward().forward().getResult();
    std::vector<double> smoothedReverseOppositeSum = Roller(window, netValues).reverse().reverse().opposite().getResult();

    PointFinder pointFinder(window);
    std::vector<size_t> forwardPeaksIndex = pointFinder.findPeaks(smoothedForwardSum, 0.25);
    std::vector<size_t> reversePeaksIndex = pointFinder.findPeaks(smoothedReverseOppositeSum, 0.25);
    std::vector<size_t> gentleRiseIndex = pointFinder.getForwardGentle(netValues, forwardPeaksIndex);
    std::vector<size_t> gentleFallIndex = pointFinder.getReverseGentle(netValues, reversePeaksIndex);

    Roller::replaceValue(clipCount, positions, gentleRiseIndex, directionalDifference);
    Roller::replaceValue(clipCount, positions, gentleFallIndex, directionalDifference);
}

void Clip::detectInterval(ClipCount &clipCount){
    int upCount = 0;
    int downCount = 0;
    int smallRegion = 10000;
    int lastPos = -smallRegion;

    bool convert = 0;
    bool push = 0;
    bool status = clipCount.begin()->second[FRONT] > clipCount.begin()->second[BACK]; // 1: up, 0: down
    int candidatePos = -1;
    int candidateCount = 0;

    int rejectCount = 0;
    int lastRejectPos = -1;

    bool artificialStatus = false;
    int artificialCount = 0;
    int artificialPos = -1;

    clipCount[clipCount.rbegin()->first + smallRegion] = clipCount.rbegin()->second;

    for(auto posIter = clipCount.begin(); posIter != clipCount.end() ; posIter++ ){
        upCount = posIter->second[FRONT];
        downCount = posIter->second[BACK];
        if (posIter->first - lastPos >= smallRegion) {
            if (candidatePos != -1 && push) {
                largeGenomicEventInterval->push_back(candidatePos);
                candidatePos = -1;
            }else if (candidatePos != -1 && !push){
                smallGenomicEventRegion->push_back(candidatePos > lastRejectPos ? std::make_pair(lastRejectPos, candidatePos) : std::make_pair(candidatePos, lastRejectPos));
                candidatePos = -1;
            }
            convert = 0;
            lastRejectPos = -1;
        }

        if (upCount >= 5 || downCount >= 5) {
            if (candidatePos == -1) {
                candidatePos = posIter->first;
                status = posIter->second[FRONT] > posIter->second[BACK];
                candidateCount = status ? upCount : downCount;
                rejectCount = candidateCount / 4;
                push = 1;
                if (status != artificialStatus && artificialCount > rejectCount && candidatePos - artificialPos < 10000) {
                    push = 0;
                    lastRejectPos = artificialPos;
                }
            }
            else {
                if (!convert){
                    if (status && downCount >= rejectCount) {
                        status = 0;
                        convert = 1;
                        push = 0;
                        lastRejectPos = posIter->first;
                    }
                    else if (!status && upCount >= rejectCount) {
                        status = 1;
                        convert = 1;
                        push = 0;
                        lastRejectPos = posIter->first;
                    }
                    else {
                        candidateCount = std::max(candidateCount, status ? upCount : downCount);
                        rejectCount = candidateCount / 4;
                        if (!push && rejectCount > artificialCount ){
                            push = 0;
                            candidatePos = posIter->first;
                        }
                        if (!status){
                            candidatePos = posIter->first;
                        }
                    }
                }
                else {
                    lastRejectPos = posIter->first;
                }
            }
            lastPos = posIter->first;
        }
        else {
            if (artificialPos == -1) {
                artificialPos = posIter->first;
                artificialCount = std::max(upCount, downCount);
                artificialStatus = (upCount > downCount);
            }
            else {
                int newMax = std::max(upCount, downCount);
                if (newMax >= artificialCount || posIter->first - artificialPos > smallRegion) {
                    artificialStatus = (upCount > downCount);
                    if (push){
                        if (status != artificialStatus){
                            artificialCount = newMax;
                        }
                    }
                    else {
                        artificialCount = newMax;
                    }
                    if (push && candidatePos - posIter->first < smallRegion && status != artificialStatus && candidatePos != -1 && artificialCount > rejectCount) {
                        push = 0;
                        lastRejectPos = posIter->first;
                    }
                    else if (!push && posIter->first - candidatePos < std::abs(artificialPos - candidatePos)){
                        lastRejectPos = posIter->first;
                    }
                    artificialPos = posIter->first;
                }
            }
        }
    }
    clipCount.erase(--clipCount.end());
    if((!largeGenomicEventInterval->empty() && largeGenomicEventInterval->front() != clipCount.begin()->first) ||
        (largeGenomicEventInterval->empty() && !clipCount.empty())){
        largeGenomicEventInterval->insert(largeGenomicEventInterval->begin(), clipCount.begin()->first);
    }
    if((!largeGenomicEventInterval->empty() && largeGenomicEventInterval->back() != clipCount.rbegin()->first) ||
        (largeGenomicEventInterval->empty() && !clipCount.empty())){
        largeGenomicEventInterval->push_back(clipCount.rbegin()->first);
    }
}

void Clip::detectLOHRegion(SnpParser &snpMap, std::vector<LOHSegment> &LOHSegments){
    std::map<int, RefAlt> *currentVariants = snpMap.getVariants(chr);

    int hetCountNum = 0;
    int homCountNum = 0;
    int filterStart = 0;
    int filterEnd = 0;
    auto intervalIter = largeGenomicEventInterval->begin();
    auto regionIter = smallGenomicEventRegion->begin();

    for (auto posIter = currentVariants->begin(); posIter != currentVariants->end(); ++posIter) {
        if (intervalIter != largeGenomicEventInterval->end() && posIter->first >= *intervalIter) {
            double genotypeRatio = (hetCountNum + homCountNum > 0) ? static_cast<double>(hetCountNum) / (homCountNum + hetCountNum) : 0.0;
            if (hetCountNum + homCountNum > 0){
                if (genotypeRatio < 0.09){
                    LOHSegments.push_back(LOHSegment(*(intervalIter-1), *intervalIter, genotypeRatio));
                }
            }
            
            // reset count
            homCountNum = 0;
            hetCountNum = 0;
            // move to next interval
            while (intervalIter != largeGenomicEventInterval->end() && posIter->first >= *intervalIter) {
                intervalIter++;
            }
        }
        while (regionIter != smallGenomicEventRegion->end() && posIter->first > regionIter->second){
            regionIter++;
        }
        if (regionIter != smallGenomicEventRegion->end()){
            filterStart = regionIter->first;
            filterEnd = regionIter->second;
        }
        if (posIter->first < filterStart || posIter->first > filterEnd){
            if (posIter->second.vaf >= 0.8 || posIter->second.homozygous){
                homCountNum++;
            }else {
                hetCountNum++;
            }
        }
    }
}
// Implementation of PurityCalculator class methods
std::map<double, int> PurityCalculator::mergeDistributionMap(const std::map<std::string, std::map<double, int>>& data) {
    std::map<double, int> distributionSumMap;
    for (const auto& chrPair : data) {
        for (const auto& valuePair : chrPair.second) {
            distributionSumMap[valuePair.first] += valuePair.second;
        }
    }
    return distributionSumMap;
}

int PurityCalculator::getTotalCount(const std::map<double, int>& data) {
    int totalCount = 0;
    for (const auto& pair : data) {
        totalCount += pair.second;
    }
    return totalCount;
}

double PurityCalculator::findQuartile(const std::map<double, int>& data, double targetPos) {
    int cumulativeCount = 0;
    double prevKey = data.begin()->first;
    int prevCumulativeCount = 0;

    for (const auto& dataIter : data) {
        cumulativeCount += dataIter.second;

        if (cumulativeCount >= targetPos) {
            if (cumulativeCount == targetPos) {
                return dataIter.first;
            }
            double fraction = (targetPos - prevCumulativeCount) / (cumulativeCount - prevCumulativeCount);
            return prevKey + fraction * (dataIter.first - prevKey);
        }

        prevCumulativeCount = cumulativeCount;
        prevKey = dataIter.first;
    }
    return prevKey;
}


double PurityCalculator::getPurity(std::map<std::string, std::map<double, int>> &inChrDistributionMap, std::string &output_root_path) {
    std::map<double, int> distributionSumMap = mergeDistributionMap(inChrDistributionMap);
    int totalCount = getTotalCount(distributionSumMap);
    double q1 = findQuartile(distributionSumMap, 0.25 * (totalCount + 1));
    double q3 = findQuartile(distributionSumMap, 0.75 * (totalCount + 1));
    std::ofstream outputFile(output_root_path + ".q1_q3");
    if (outputFile.is_open()) {
        outputFile << q1 << "\t" << q3 << "\n";
        outputFile.close();
    }
    double purity = -5.3134 + 11.5568*q1 + 2.1985*q3 - 39.4693*q1*q1 + 47.8858*q1*q3 - 18.2485*q3*q3;
    // Clamp the purity value between 0 and 1
    return std::max(0.0, std::min(purity, 1.0));
}
