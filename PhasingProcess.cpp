#include "PhasingProcess.h"
#include "PhasingGraph.h"
#include "ParsingBam.h"

PhasingProcess::PhasingProcess(PhasingParameters params)
{
    std::cerr<< "LongPhase Ver " << params.version << "\n";
    std::cerr<< "\n";
    std::cerr<< "--- File Parameter --- \n";
    std::cerr<< "SNP File       : " << params.snpFile       << "\n";
    std::cerr<< "SV  File       : " << params.svFile        << "\n";
    std::cerr<< "MOD File       : " << params.modFile       << "\n";
    std::cerr<< "PON File       : " << params.ponFile       << "\n";
    std::cerr<< "Strict PON File: " << params.strictPonFile << "\n";
    std::cerr<< "REF File       : " << params.fastaFile     << "\n";
    std::cerr<< "Output Prefix  : " << params.resultPrefix  << "\n";
    std::cerr<< "Generate Dot   : " << ( params.generateDot ? "True" : "False" ) << "\n";
    std::cerr<< "Output LOH     : " << ( params.outputLOH ? "True" : "False" ) << "\n";
    std::cerr<< "Output SGE     : " << ( params.outputSGE ? "True" : "False" ) << "\n";
    std::cerr<< "Output LGE     : " << ( params.outputLGE ? "True" : "False" ) << "\n";
    std::cerr<< "Output GE      : " << ( params.outputGE ? "True" : "False" ) << "\n";
    std::cerr<< "BAM File       : ";
    for( auto file : params.bamFile){
        std::cerr<< file <<" " ;   
    }
    std::cerr << "\n";
    
    std::cerr<< "\n";
    std::cerr<< "--- Phasing Parameter --- \n";
    std::cerr<< "Seq Platform       : " << ( params.isONT ? "ONT" : "PB" ) << "\n";
    std::cerr<< "Phase Indel        : " << ( params.phaseIndel ? "True" : "False" )  << "\n";
    std::cerr<< "Distance Threshold : " << params.distance        << "\n";
    std::cerr<< "Connect Adjacent   : " << params.connectAdjacent << "\n";
    std::cerr<< "Somatic Connect Adjacent   : " << params.somaticConnectAdjacent << "\n";
    std::cerr<< "Edge Threshold     : " << params.edgeThreshold   << "\n";
    std::cerr<< "Overlap Threshold  : " << params.overlapThreshold   << "\n";
    std::cerr<< "Mapping Quality    : " << params.mappingQuality  << "\n";
    std::cerr<< "Mismatch Rate      : " << params.mismatchRate  << "\n";
    std::cerr<< "Variant Confidence : " << params.snpConfidence   << "\n";
    std::cerr<< "ReadTag Confidence : " << params.readConfidence  << "\n";
    std::cerr<< "\n";
    
    std::time_t processBegin = time(NULL);
        
    // load SNP vcf file
    std::time_t begin = time(NULL);
    std::cerr<< "parsing VCF ... ";
    SnpParser snpFile(params);
    std::cerr<< difftime(time(NULL), begin) << "s\n";

    // load PON vcf and strict PON vcf file
    begin = time(NULL);
    std::cerr<< "parsing PON VCF ... ";
    snpFile.setGermline(params.ponFile, params.strictPonFile);
    std::cerr<< difftime(time(NULL), begin) << "s\n";

    // load SV vcf file
    begin = time(NULL);
    std::cerr<< "parsing SV VCF ... ";
    SVParser svFile(params, snpFile);
    std::cerr<< difftime(time(NULL), begin) << "s\n";
 
    //Parse mod vcf file
	begin = time(NULL);
	std::cerr<< "parsing Meth VCF ... ";
    METHParser modFile(params, snpFile, svFile);
	std::cerr<< difftime(time(NULL), begin) << "s\n";
 
    // parsing ref fasta 
    begin = time(NULL);
    std::cerr<< "reading reference ... ";
    std::vector<int> last_pos;
    for(auto chr :snpFile.getChrVec()){
        last_pos.push_back(snpFile.getLastSNP(chr));
    }
    FastaParser fastaParser(params.fastaFile, snpFile.getChrVec(), last_pos, params.numThreads);
    std::cerr<< difftime(time(NULL), begin) << "s\n";

    // get all detected chromosome
    std::vector<std::string> chrName = snpFile.getChrVec();

    // record chromosome info
    std::map<std::string, ChrInfo> chrInfoMap;
    // Initialize chrInfoMap entries for each chromosome name in advance
    // to avoid potential issues during parallel processing later
    for (std::vector<std::string>::iterator chrIter = chrName.begin(); chrIter != chrName.end(); chrIter++)    {
        chrInfoMap[*chrIter] = ChrInfo();
    }

    // init data structure and get core n
    htsThreadPool threadPool = {NULL, 0};

    // creat thread pool
    if (!(threadPool.pool = hts_tpool_init(params.numThreads))) {
        fprintf(stderr, "Error creating thread pool\n");
    }

    std::cerr << "parsing BAM, somatic calling, and phasing" << std::endl;
    begin = time(NULL);
    
    // loop all chromosome
    #pragma omp parallel for schedule(dynamic) num_threads(params.numThreads)
    for(std::vector<std::string>::iterator chrIter = chrName.begin(); chrIter != chrName.end() ; chrIter++ ){
        
        std::time_t chrbegin = time(NULL);
        ChrInfo &chrInfo = chrInfoMap[*chrIter];
        // get lase SNP variant position
        int lastSNPpos = snpFile.getLastSNP((*chrIter));
        // therer is no variant on SNP file. 
        if( lastSNPpos == -1 ){
            continue;
        }

	    // fetch chromosome string
        std::string chr_reference = fastaParser.chrString.at(*chrIter);
        // create a bam parser object and prepare to fetch varint from each vcf file
	    BamParser *bamParser = new BamParser((*chrIter), params.bamFile, snpFile, svFile, modFile, chr_reference);
        // use to store variant
        std::vector<ReadVariant> *readVariantVec = new std::vector<ReadVariant>();
        // run fetch variant process
        bamParser->direct_detect_alleles(lastSNPpos, threadPool, params, *readVariantVec, chrInfo.clipCount, chr_reference);
        // free memory
        delete bamParser;
        
        // filter variants prone to switch errors in ONT sequencing.
        if(params.isONT){
            snpFile.filterSNP((*chrIter), *readVariantVec, chr_reference);
        }

        // bam files are partial file or no read support this chromosome's SNP
        if( readVariantVec->size() == 0 ){
            delete readVariantVec;
            continue;
        }

        // create a clip object and prepare to detect Interval
        Clip *clip = new Clip(*chrIter);
        // get the interval of the genomic event
        clip->detectGenomicEventInterval(chrInfo.clipCount, chrInfo.largeGenomicEventInterval, chrInfo.smallGenomicEventRegion);
        // get the region of the LOH
        clip->detectLOHRegion(snpFile, chrInfo.LOHSegments);
        // free memory
        delete clip;
        
        // create a graph object and prepare to phasing.
        VairiantGraph *vGraph = new VairiantGraph(chr_reference, params, (*chrIter));
        chrInfo.vGraph = vGraph;
        // trans read-snp info to edge info
        vGraph->addEdge(readVariantVec);
        // run somatic calling algorithm
        vGraph->somaticCalling(snpFile.getVariants((*chrIter)));
        // run main algorithm
        vGraph->phasingProcess(chrInfo.posPhasingResult, chrInfo.LOHSegments, &chrInfo.ploidyRatioMap);
        std::cerr<< "(" << (*chrIter) << "," << difftime(time(NULL), chrbegin) << "s)";
    }

    std::map<std::string, std::map<double, int>> mergedPloidyRatioMap;
    for (const auto& chrInfo : chrInfoMap) {
        mergedPloidyRatioMap[chrInfo.first] = chrInfo.second.ploidyRatioMap;
    }
    double purity = PurityCalculator::getPurity(mergedPloidyRatioMap, params.resultPrefix);
    std::cerr << std::endl;
    std::cerr << "purity: " << purity << std::endl;
    if(purity > 0.95){
        std::cerr << "second round phasing, ";
    }
    std::cerr << "export phasing result" << std::endl;

    #pragma omp parallel for schedule(dynamic) num_threads(params.numThreads)
    for(std::vector<std::string>::iterator chrIter = chrName.begin(); chrIter != chrName.end() ; chrIter++ ){
        std::time_t chrbegin = time(NULL);
        ChrInfo &chrInfo = chrInfoMap[*chrIter];
        if(chrInfo.vGraph == nullptr)continue;

        VairiantGraph *vGraph = chrInfo.vGraph;
        if(purity > 0.95){
            // convert non-germline variants to somatic variants
            vGraph->convertNonGermlineToSomatic();
            // reset phasing result
            chrInfo.posPhasingResult = PosPhasingResult();
            // run main algorithm
            vGraph->phasingProcess(chrInfo.posPhasingResult, chrInfo.LOHSegments, nullptr);
        }
        // export phasing result
        vGraph->exportPhasingResult(chrInfo.posPhasingResult, chrInfo.LOHSegments);
        // generate dot file
        if(params.generateDot){
            vGraph->writingDotFile((*chrIter));
        }
        
        // release the memory used by the object.
        vGraph->destroy();
        delete vGraph;
        
        std::cerr<< "(" << (*chrIter) << "," << difftime(time(NULL), chrbegin) << "s)";
    }
    hts_tpool_destroy(threadPool.pool);
    
    // Transfer phasing results from chrInfoMap to chrPhasingResult
    ChrPhasingResult chrPhasingResult;
    for (const auto& chrInfo : chrInfoMap) {
        chrPhasingResult[chrInfo.first] = chrInfo.second.posPhasingResult;
    }

    std::cerr<< "\nparsing total:  " << difftime(time(NULL), begin) << "s\n";

    // write result to file
    GenomicWriter genomicWriter(params.resultPrefix, chrName, chrInfoMap);
    genomicWriter.measureTime("LOH", params.outputLOH, [&]() { genomicWriter.writeLOH(); });
    genomicWriter.measureTime("SGE", params.outputSGE, [&]() { genomicWriter.writeSGE(); });
    genomicWriter.measureTime("LGE", params.outputLGE, [&]() { genomicWriter.writeLGE(); });
    genomicWriter.measureTime("GE", params.outputGE, [&]() { genomicWriter.writeAllEvents(); });
    genomicWriter.measureTime("SNP", true, [&]() { snpFile.writeResult(chrPhasingResult); });
    genomicWriter.measureTime("SV", params.svFile != "", [&]() { svFile.writeResult(chrPhasingResult); });
    genomicWriter.measureTime("MOD", params.modFile != "", [&]() { modFile.writeResult(chrPhasingResult); });

    std::cerr<< "\ntotal process: " << difftime(time(NULL), processBegin) << "s\n";

    return;
};

PhasingProcess::~PhasingProcess(){
};

