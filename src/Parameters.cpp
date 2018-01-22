//
//  Parameters.cpp
//  ReadTest
//
//  Created by 臼山 直人 on 1/13/14.
//  Copyright (c) 2014 Univ. of Tokyo. All rights reserved.
//

#include <cstdlib>
#include "Parameters.h"
#include "cmdline.h"
#include <algorithm>


Parameters::Parameters(int argc, const char *argv[]) {
    std::cout << "set all defalt value done" << std::endl;
    getFromCommandLineArguments(argc, argv);
    std::cout << "getCommandLineArg value done" << std::endl;
    std::cout << tumorBam << std::endl;
    std::cout << priorSNP << std::endl;
}


inline void Parameters::parseHyperParameters(std::vector<double> &vec, std::string s) {
    std::vector<std::string> params = Utils::split(s, ",");
    vec.clear();
    for (int i = 0; i < params.size(); i++) {
        vec.push_back(atof(params[i].c_str()));
    }
}

void Parameters::getFromCommandLineArguments(int argc, const char *argv[]) {
    cmdline::parser a;
    a.add<std::string>("ref", 'f', "fasta reference sequence (should be indexed with .fai file)", true, "");
    a.add<std::string>("tumor", 'a', "tumor bam file", true, "");
    a.add<std::string>("normal", 'b', "normal bam file", true, "");
    a.add<std::string>("out", 'o', "file-prefix for output results", true, "");
    a.add<std::string>("windows", 'w', "file with candidate variants to be tested.", false, "");
    
    a.add<std::string>("region", 'R', "pileUp region.", false, "");
    a.add<std::string>("pileup", 'P', "pileUp file", false, "");

    
    a.add("quiet", '\0', "quiet output");
    
    a.add<double>("priorSNP", '\0', "prior probability of a SNP site", false, 1.0 / 100.0);
    a.add<double>("priorIndel", '\0', "prior probability of a detected indel not being a sequencing error", false, 1.0 / 10000.0);
    a.add<int>("maxReads", '\0', "maximum number of reads in a window", false, 20000);
    a.add<int>("minReads", '\0', "minimum number of reads", false, 4);
    // a.add<double>("averageMapQualThreshold", '\0', "lower limit of average read mapping quality for considering haplotypes", false, 0.0);
    // a.add<double>("mapQualThreshold", '\0', "lower limit of read mapping quality", false, 0.9683772233983162);
    a.add<int>("averageMapPhredQualThreshold", '\0', "lower limit of average read mapping quality for considering haplotypes", false, 0);
    a.add<int>("mapPhredQualThreshold", '\0', "lower limit of read mapping quality", false, 15);


    a.add<int>("maxReadLength", '\0', "maximum length of reads", false, 200);
    
    a.add<double>("pError", '\0', "probability of a read indel", false, 5e-4);
    a.add<double>("pMut", '\0', "probability of a mutation in the read", false, 1e-5);
    a.add<int>("maxIndelLength", '\0', "maximum length of a _sequencing error_ indel in read", false, 25);
    
    a.add<int>("maxInsertSize", '\0', "maximum insert size of paired-end reads", false, 800);
    
    a.add<int>("minDistanceSNP", '\0', "minimum distance to the heterozygous SNP", false, 3);
    a.add<int>("minDistanceGermlineIndel", '\0', "minimum distance to the heterozygous indel", false, 25);
    a.add("singleReads", '\0', "Do not use paired-end information");
    a.add<double>("indelPosessionFreqThreshold", '\0', "Threshold for frequency of reads that have indels", false, 1.0);
    a.add<double>("softClipPosessionFreqThreshold", '\0', "Threshold for frequency of reads that have SoftClip", false, 1.0);
    
    a.add("withoutBayesFactor", '\0', "Do not calculate Bayes factor, only shows candidate mutations ");
    
    a.add<std::string>("algorithm", '\0', " algorithm type used in mutation Call, choose one of them : OHVarfinDer2, HapMuC ", false, "OHVarfinDer2");
    a.add<double>("heteroSNPConfidenceInterval", '\0', " check, (ref,obs) is within the ${heteroSNPConfidenceInterval} confidence interval of binom(ref+obs,obs,0.5)", false,0.99);
    a.add<int>("pileUpBufferSize", '\0', "pileUp Parse BufferSize", false, 4000000);
    
    a.add<int>("tumorMinDepth", '\0', "tumor pileUp minDepth ", false, 12);
    a.add<double>("tumorMinObsRate", '\0', "tumor pileUp minObsRate ", false, 0.05);
    a.add<double>("tumorMaxObsRate", '\0', "tumor pileUp maxObsRate ", false, 1.0);
    a.add<int>("tumorMinObsNum", '\0', "tumor pileUp minObsNum ", false, 4);
    a.add<int>("tumorMaxObsNum", '\0', "tumor pileUp maxObsNum ", false, 1000000);
    a.add<double>("tumorMinAvgBaseQuality", '\0', "tumor pileUp min avgBaseQuality for somatic", false, 0);

    a.add<int>("normalMinDepth", '\0', "normal pileUp minDepth ", false, 12);
    a.add<double>("normalMinObsRate", '\0', "normal pileUp minObsRate ", false, 0.0);
    a.add<double>("normalMaxObsRate", '\0', "normal pileUp maxObsRate ", false, 0.1);
    a.add<int>("normalMinObsNum", '\0', "normal pileUp minObsNum ", false, 0);
    a.add<int>("normalMaxObsNum", '\0', "normal pileUp maxObsNum ", false, 1000000);
    a.add<double>("normalMinAvgBaseQuality", '\0', "normal pileUp min avgBaseQuality for somatic", false, 0);
    
    a.add<int>("heteroSNPMinDepth", '\0', "heteroSNP pileUp minDepth ", false, 15);
    a.add<int>("heteroSNPMaxDepth", '\0', "heteroSNP pileUp maxDepth ", false, 1000000);
    a.add<double>("heteroSNPMinObsRate", '\0', "heteroSNP pileUp minObsRate ", false, 0.1);
    a.add<double>("heteroSNPMaxObsRate", '\0', "heteroSNP pileUp maxObsRate ", false, 1.0);
    a.add<int>("heteroSNPMinObsNum", '\0', "heteroSNP pileUp minObsNum ", false, 5);
    a.add<int>("heteroSNPMaxObsNum", '\0', "heteroSNP pileUp maxObsNum ", false, 1000000);
    a.add<double>("heteroSNPMinAvgBaseQuality", '\0', "normal pileUp min avgBaseQuality for SNP", false, 0);

    a.add<double>("triAlleleMinObsRate", '\0', "triAllele pileUp minObsRate ", false, 0.03);
    a.add<int>("triAlleleMinObsNum", '\0', "triAllele pileUp minObsNum ", false, 5);
    
    a.add<int>("minBQ", '\0', " min Base quality used in PileUp", false, 15);
    a.add<int>("fReadFilter", '\0', " -f read filter used in read Collection", false, 2);
    a.add<int>("FReadFilter", '\0', " -F read filter used in read Collection", false, 3840);
    

        

    a.add<std::string>("hapmuc_mutationModelAlpha", '\0', "hyper parameter for tumor haplotype frequencies in mutation model.", false, "1.0,1.0,1.0");
    a.add<std::string>("hapmuc_mutationModelBeta", '\0', "hyper parameter for error rates in mutation model.", false, "0.1,10.0");
    a.add<std::string>("hapmuc_mutationModelGamma", '\0', "hyper parameter for normal haplotype frequencies in mutation model.", false, "1.0,1.0");
    
    a.add<std::string>("hapmuc_errorModelAlpha", '\0', "hyper parameter for tumor haplotype frequencies in error model.", false, "1.0,1.0");
    a.add<std::string>("hapmuc_errorModelBeta", '\0', "hyper parameter for error rates in error model.", false, "1.0,10.0");
    a.add<std::string>("hapmuc_errorModelGamma", '\0', "hyper parameter for normal haplotype frequencies in error model.", false, "1.0,1.0");

    
    a.add<std::string>("ohvar2_mutGammaF", '\0', "ohvar2_mutGammaF", false, "10.0,1.0");
    a.add<std::string>("ohvar2_mutGammaH", '\0', "ohvar2_mutGammaH", false, "5.0,5.0,1.0");
    a.add<std::string>("ohvar2_mutAlphaL", '\0', "ohvar2_mutAlphaL", false, "1.0,100.0");
    a.add<std::string>("ohvar2_mutAlphaH", '\0', "ohvar2_mutAlphaH", false, "1.0,100.0");
    a.add<std::string>("ohvar2_mutGammaEH", '\0', "ohvar2_mutGammaEH",false, "5.0,5.0");
    a.add<std::string>("ohvar2_mutAlphaS",  '\0', "ohvar2_mutAlphaS", false, "1.0,100.0");
    a.add<std::string>("ohvar2_mutAlphaB_E", '\0', "ohvar2_mutAlphaB_E", false, "50.0,50.0");
    a.add<std::string>("ohvar2_mutAlphaB_W", '\0', "ohvar2_mutAlphaB_W", false, "5.0,5.0");
    
    a.add<std::string>("ohvar2_errAlphaL_E",  '\0', "ohvar2_errAlphaL_E", false, "2.0,30.0");
    a.add<std::string>("ohvar2_errAlphaH_E",  '\0', "ohvar2_errAlphaH_E", false, "2.0,30.0");
    a.add<std::string>("ohvar2_errGammaEH_E", '\0', "ohvar2_errGammaEH_E",false, "5.0,5.0");
    a.add<std::string>("ohvar2_errAlphaS_E",  '\0', "ohvar2_errAlphaS_E", false, "2.0,30.0");
    a.add<std::string>("ohvar2_errAlphaB_E",  '\0', "ohvar2_errAlphaB_E", false, "0.05,0.05");

    a.add<std::string>("ohvar2_errAlphaL_W",  '\0', "ohvar2_errAlphaL_W", false, "2.0,30.0");
    a.add<std::string>("ohvar2_errAlphaH_W",  '\0', "ohvar2_errAlphaH_W", false, "2.0,30.0");
    a.add<std::string>("ohvar2_errGammaEH_W", '\0', "ohvar2_errGammaEH_W",false, "5.0,5.0");
    a.add<std::string>("ohvar2_errAlphaS_W",  '\0', "ohvar2_errAlphaS_W", false, "2.0,30.0");
    a.add<std::string>("ohvar2_errAlphaB_W",  '\0', "ohvar2_errAlphaB_W", false, "0.5,0.5");
    
    a.add<std::string>("ohvar2_noCategory_mutGammaF", '\0', "ohvar2_noCategory_mutGammaF", false, "10.0,1.0");
    a.add<std::string>("ohvar2_noCategory_mutGammaH", '\0', "ohvar2_noCategory_mutGammaH", false, "5.0,5.0,1.0");
    a.add<std::string>("ohvar2_noCategory_mutAlphaL", '\0', "ohvar2_noCategory_mutAlphaL", false, "1.0,100.0");
    a.add<std::string>("ohvar2_noCategory_mutAlphaH", '\0', "ohvar2_noCategory_mutAlphaH", false, "1.0,100.0");
    a.add<std::string>("ohvar2_noCategory_mutAlphaB", '\0', "ohvar2_noCategory_mutAlphaB", false, "5.0,5.0");
    a.add<std::string>("ohvar2_noCategory_mutGammaEH", '\0', "ohvar2_noCategory_mutGammaEH",false, "5.0,5.0");
    a.add<std::string>("ohvar2_noCategory_mutAlphaS",  '\0', "ohvar2_noCategory_mutAlphaS", false, "1.0,100.0");

    a.add<std::string>("ohvar2_noCategory_mutAlphaB_E",  '\0', "ohvar2_noCategory_mutAlphaB_N", false, "50.0,50.0");
    a.add<std::string>("ohvar2_noCategory_mutAlphaB_W",  '\0', "ohvar2_noCategory_mutAlphaB_N", false, "5.0,5.0");

    a.add<std::string>("ohvar2_noCategory_errAlphaL_E",  '\0', "ohvar2_noCategory_errAlphaL_E", false, "2.0,30.0");
    a.add<std::string>("ohvar2_noCategory_errAlphaH_E",  '\0', "ohvar2_noCategory_errAlphaH_E", false, "2.0,30.0");
    a.add<std::string>("ohvar2_noCategory_errGammaEH_E", '\0', "ohvar2_noCategory_errGammaEH_E",false, "5.0,5.0");
    a.add<std::string>("ohvar2_noCategory_errAlphaS_E",  '\0', "ohvar2_noCategory_errAlphaS_E", false, "2.0,30.0");
    a.add<std::string>("ohvar2_noCategory_errAlphaB_E",  '\0', "ohvar2_noCategory_errAlphaB_E", false, "0.5,0.5");

    a.add<std::string>("ohvar2_noCategory_errAlphaL_W",  '\0', "ohvar2_noCategory_errAlphaL_W", false, "2.0,30.0");
    a.add<std::string>("ohvar2_noCategory_errAlphaH_W",  '\0', "ohvar2_noCategory_errAlphaH_W", false, "2.0,30.0");
    a.add<std::string>("ohvar2_noCategory_errGammaEH_W", '\0', "ohvar2_noCategory_errGammaEH_W",false, "5.0,5.0");
    a.add<std::string>("ohvar2_noCategory_errAlphaS_W",  '\0', "ohvar2_noCategory_errAlphaS_W", false, "2.0,30.0");
    a.add<std::string>("ohvar2_noCategory_errAlphaB_W",  '\0', "ohvar2_noCategory_errAlphaB_W", false, "0.5,0.5");
    a.parse_check(argc, argv);
    
    maxReads = a.get<int>("maxReads");
    minReads = a.get<int>("minReads");

    averageMapPhredQualThreshold = a.get<int>("averageMapPhredQualThreshold");
    mapPhredQualThreshold = a.get<int>("mapPhredQualThreshold");
    averageMapQualThreshold = 1.0 - pow(10, -1.0 * averageMapPhredQualThreshold/10.0 );
    mapQualThreshold        = 1.0 - pow(10, -1.0 * mapPhredQualThreshold/10.0 );

    // averageMapQualThreshold = a.get<double>("averageMapQualThreshold");
    // mapQualThreshold = a.get<double>("mapQualThreshold");
    maxReadLength = a.get<int>("maxReadLength");
    priorSNP = a.get<double>("priorSNP");
    priorIndel = a.get<double>("priorIndel");
    refFileName = a.get<std::string>("ref");
    outFilePrefix = a.get<std::string>("out");
    tumorBam = a.get<std::string>("tumor");
    normalBam = a.get<std::string>("normal");
    windowFile = a.get<std::string>("windows");
    pError = a.get<double>("pError");
    pMut = a.get<double>("pMut");
    maxIndelLength = a.get<int>("maxIndelLength");
    minDistanceSNP = a.get<int>("minDistanceSNP");
    minDistanceGermlineIndel = a.get<int>("minDistanceGermlineIndel");
    quiet = a.exist("quiet") ? true : false;
    maxInsertSize = a.get<int>("maxInsertSize");
    singleReads = a.exist("singleReads") ? true : false;
    
    softClipPosessionFreqThreshold = a.get<double>("softClipPosessionFreqThreshold");
    indelPosessionFreqThreshold = a.get<double>("indelPosessionFreqThreshold");
    if (mapQualThreshold > 1.0) {
        throw std::string("mapQualThreshold should be 0.0 ~ 1.0. ex: 30 -> 0.999");
    }
    withoutBayesFactor = a.exist("withoutBayesFactor") ? true : false;
    
    pileupFile = a.get<std::string>("pileup");
    region = a.get<std::string>("region");
    
    heteroSNPConfidenceInterval = a.get<double>("heteroSNPConfidenceInterval");

    pileUpBufferSize = a.get<int>("pileUpBufferSize");
    fReadFilter = a.get<int>("fReadFilter");
    FReadFilter = a.get<int>("FReadFilter");
    if(singleReads == true){ fReadFilter &= (~3);}

    pileT = Parameters::PileUpParameters(PileUpParameters::TUMOR);
    pileN = Parameters::PileUpParameters(PileUpParameters::NORMAL);
    pileHetero    = Parameters::PileUpParameters(PileUpParameters::HETEROSNP);
    pileTriAllele = Parameters::PileUpParameters(PileUpParameters::TRIALLELE);
    

    pileT.minDepth = a.get<int>("tumorMinDepth");
    pileT.minObsRate = a.get<double>("tumorMinObsRate");
    pileT.maxObsRate = a.get<double>("tumorMaxObsRate");
    pileT.minObsNum = a.get<int>("tumorMinObsNum");
    pileT.maxObsNum = a.get<int>("tumorMaxObsNum");
    pileT.avgBaseQualityThreshold = a.get<double>("tumorMinAvgBaseQuality");
    
    pileN.minDepth = a.get<int>("normalMinDepth");
    pileN.minObsRate = a.get<double>("normalMinObsRate");
    pileN.maxObsRate = a.get<double>("normalMaxObsRate");
    pileN.minObsNum = a.get<int>("normalMinObsNum");
    pileN.maxObsNum = a.get<int>("normalMaxObsNum");
    pileN.avgBaseQualityThreshold = a.get<double>("normalMinAvgBaseQuality");
    
    pileHetero.minDepth = a.get<int>("heteroSNPMinDepth");
    pileHetero.maxDepth = a.get<int>("heteroSNPMaxDepth");
    pileHetero.minObsRate = a.get<double>("heteroSNPMinObsRate");
    pileHetero.maxObsRate = a.get<double>("heteroSNPMaxObsRate");
    pileHetero.minObsNum = a.get<int>("heteroSNPMinObsNum");
    pileHetero.maxObsNum = a.get<int>("heteroSNPMaxObsNum");
    pileHetero.avgBaseQualityThreshold = a.get<double>("heteroSNPMinAvgBaseQuality");

    pileTriAllele.minDepth   = pileT.minDepth;
    pileTriAllele.minObsRate = std::min(a.get<double>("triAlleleMinObsRate"), pileT.minObsRate);
    pileTriAllele.maxObsRate = 1.0;
    pileTriAllele.minObsNum  = std::min(a.get<int>("triAlleleMinObsNum"), pileT.minObsNum);
    pileTriAllele.maxObsNum  = 1000000;
    
    std::string algorithmType = a.get<std::string>("algorithm");
    if(     algorithmType == "OHVarfinDer2")           method = Parameters::OHVARFINDER2;
    else if(algorithmType == "HapMuC")                 method = Parameters::HAPMUC;
    else if(algorithmType == "HeteroSNPCall")          method = Parameters::HETEROSNPCALL;
    else if(algorithmType == "ReadProbCounter")        method = Parameters::READPROBCOUNTER;
    else if(algorithmType == "OHVarfinDer2NoCategory") method = Parameters::OHVARFINDER2NOCATEGORY;
    else                                               throw std::string(" unknown Algorithm type specified ");
    
    if(method == Parameters::OHVARFINDER2)                getFromCommandLineArgumentsInOHVarfinDer2(a,argc,argv);
    else if(method == Parameters::Parameters::HAPMUC)     getFromCommandLineArgumentsInHapMuC(a,argc,argv);
    else if(method == Parameters::OHVARFINDER2NOCATEGORY) getFromCommandLineArgumentsInOHVarfinDer2NoCategory(a,argc,argv);
}


void Parameters::getFromCommandLineArgumentsInHapMuC(cmdline::parser& a, int argc, const char *argv[]){
    parseHyperParameters(hap3_params.mut_a0, a.get<std::string>("hapmuc_mutationModelAlpha"));
    parseHyperParameters(hap3_params.mut_b0, a.get<std::string>("hapmuc_mutationModelBeta"));
    parseHyperParameters(hap3_params.mut_c0, a.get<std::string>("hapmuc_mutationModelGamma"));
    parseHyperParameters(hap3_params.err_a0, a.get<std::string>("hapmuc_errorModelAlpha"));
    parseHyperParameters(hap3_params.err_b0, a.get<std::string>("hapmuc_errorModelBeta"));
    parseHyperParameters(hap3_params.err_c0, a.get<std::string>("hapmuc_errorModelGamma"));
    hap2_params = hap3_params;
}

void Parameters::getFromCommandLineArgumentsInOHVarfinDer2(cmdline::parser& a, int argc, const char *argv[]){
    OHVarfinDParams2.clear_all();
    parseHyperParameters(OHVarfinDParams2.mut_a0, a.get<std::string>("ohvar2_mutGammaF"));
    parseHyperParameters(OHVarfinDParams2.mut_b0, a.get<std::string>("ohvar2_mutGammaH"));
    parseHyperParameters(OHVarfinDParams2.mut_c0, a.get<std::string>("ohvar2_mutAlphaL"));
    parseHyperParameters(OHVarfinDParams2.mut_d0, a.get<std::string>("ohvar2_mutAlphaH"));
    parseHyperParameters(OHVarfinDParams2.mut_e0, a.get<std::string>("ohvar2_mutGammaEH"));
    parseHyperParameters(OHVarfinDParams2.mut_f0, a.get<std::string>("ohvar2_mutAlphaS"));
    parseHyperParameters(OHVarfinDParams2.mut_g0, a.get<std::string>("ohvar2_mutAlphaB_E"));
    parseHyperParameters(OHVarfinDParams2.mut_h0, a.get<std::string>("ohvar2_mutAlphaB_W"));

    
    parseHyperParameters(OHVarfinDParams2.err_a0, a.get<std::string>("ohvar2_errAlphaL_E"));
    parseHyperParameters(OHVarfinDParams2.err_b0, a.get<std::string>("ohvar2_errAlphaH_E"));
    parseHyperParameters(OHVarfinDParams2.err_c0, a.get<std::string>("ohvar2_errGammaEH_E"));
    parseHyperParameters(OHVarfinDParams2.err_d0, a.get<std::string>("ohvar2_errAlphaS_E"));
    parseHyperParameters(OHVarfinDParams2.err_e0, a.get<std::string>("ohvar2_errAlphaB_E"));

    parseHyperParameters(OHVarfinDParams2.err_f0, a.get<std::string>("ohvar2_errAlphaL_W"));
    parseHyperParameters(OHVarfinDParams2.err_g0, a.get<std::string>("ohvar2_errAlphaH_W"));
    parseHyperParameters(OHVarfinDParams2.err_h0, a.get<std::string>("ohvar2_errGammaEH_W"));
    parseHyperParameters(OHVarfinDParams2.err_i0, a.get<std::string>("ohvar2_errAlphaS_W"));
    parseHyperParameters(OHVarfinDParams2.err_j0, a.get<std::string>("ohvar2_errAlphaB_W"));
}

void Parameters::getFromCommandLineArgumentsInOHVarfinDer2NoCategory(cmdline::parser& a, int argc, const char *argv[]){
    OHVarfinDParams2NoCategory.clear_all();
    parseHyperParameters(OHVarfinDParams2NoCategory.mut_a0, a.get<std::string>("ohvar2_noCategory_mutGammaF"));
    parseHyperParameters(OHVarfinDParams2NoCategory.mut_b0, a.get<std::string>("ohvar2_noCategory_mutGammaH"));
    parseHyperParameters(OHVarfinDParams2NoCategory.mut_c0, a.get<std::string>("ohvar2_noCategory_mutAlphaL"));
    parseHyperParameters(OHVarfinDParams2NoCategory.mut_d0, a.get<std::string>("ohvar2_noCategory_mutAlphaH"));
    parseHyperParameters(OHVarfinDParams2NoCategory.mut_e0, a.get<std::string>("ohvar2_noCategory_mutGammaEH"));
    parseHyperParameters(OHVarfinDParams2NoCategory.mut_f0, a.get<std::string>("ohvar2_noCategory_mutAlphaS"));
    parseHyperParameters(OHVarfinDParams2NoCategory.mut_g0, a.get<std::string>("ohvar2_noCategory_mutAlphaB_E"));
    parseHyperParameters(OHVarfinDParams2NoCategory.mut_h0, a.get<std::string>("ohvar2_noCategory_mutAlphaB_W"));

    parseHyperParameters(OHVarfinDParams2NoCategory.err_a0, a.get<std::string>("ohvar2_noCategory_errAlphaL_E"));
    parseHyperParameters(OHVarfinDParams2NoCategory.err_b0, a.get<std::string>("ohvar2_noCategory_errAlphaH_E"));
    parseHyperParameters(OHVarfinDParams2NoCategory.err_c0, a.get<std::string>("ohvar2_noCategory_errGammaEH_E"));
    parseHyperParameters(OHVarfinDParams2NoCategory.err_d0, a.get<std::string>("ohvar2_noCategory_errAlphaS_E"));
    parseHyperParameters(OHVarfinDParams2NoCategory.err_e0, a.get<std::string>("ohvar2_noCategory_errAlphaB_E"));

    parseHyperParameters(OHVarfinDParams2NoCategory.err_f0, a.get<std::string>("ohvar2_noCategory_errAlphaL_W"));
    parseHyperParameters(OHVarfinDParams2NoCategory.err_g0, a.get<std::string>("ohvar2_noCategory_errAlphaH_W"));
    parseHyperParameters(OHVarfinDParams2NoCategory.err_h0, a.get<std::string>("ohvar2_noCategory_errGammaEH_W"));
    parseHyperParameters(OHVarfinDParams2NoCategory.err_i0, a.get<std::string>("ohvar2_noCategory_errAlphaS_W"));
    parseHyperParameters(OHVarfinDParams2NoCategory.err_j0, a.get<std::string>("ohvar2_noCategory_errAlphaB_W"));


/*
    parseHyperParameters(OHVarfinDParams2NoCategory.mut_a0, a.get<std::string>("ohvar2_noCategory_mutGammaF_T"));
    parseHyperParameters(OHVarfinDParams2NoCategory.mut_b0, a.get<std::string>("ohvar2_noCategory_mutGammaH_T"));
    parseHyperParameters(OHVarfinDParams2NoCategory.mut_c0, a.get<std::string>("ohvar2_noCategory_mutAlphaL_T"));
    parseHyperParameters(OHVarfinDParams2NoCategory.mut_d0, a.get<std::string>("ohvar2_noCategory_mutAlphaH_T"));
    parseHyperParameters(OHVarfinDParams2NoCategory.mut_e0, a.get<std::string>("ohvar2_noCategory_mutAlphaB_T"));

    parseHyperParameters(OHVarfinDParams2NoCategory.mut_f0, a.get<std::string>("ohvar2_noCategory_mutAlphaL_N"));
    parseHyperParameters(OHVarfinDParams2NoCategory.mut_g0, a.get<std::string>("ohvar2_noCategory_mutAlphaH_N"));
    parseHyperParameters(OHVarfinDParams2NoCategory.mut_h0, a.get<std::string>("ohvar2_noCategory_mutGammaEH_N"));
    parseHyperParameters(OHVarfinDParams2NoCategory.mut_i0, a.get<std::string>("ohvar2_noCategory_mutAlphaS_N"));
    parseHyperParameters(OHVarfinDParams2NoCategory.mut_j0, a.get<std::string>("ohvar2_noCategory_mutAlphaB_N"));    
    
    parseHyperParameters(OHVarfinDParams2NoCategory.err_a0, a.get<std::string>("ohvar2_noCategory_errAlphaL_T"));
    parseHyperParameters(OHVarfinDParams2NoCategory.err_b0, a.get<std::string>("ohvar2_noCategory_errAlphaH_T"));
    parseHyperParameters(OHVarfinDParams2NoCategory.err_c0, a.get<std::string>("ohvar2_noCategory_errGammaEH_T"));
    parseHyperParameters(OHVarfinDParams2NoCategory.err_d0, a.get<std::string>("ohvar2_noCategory_errAlphaS_T"));
    parseHyperParameters(OHVarfinDParams2NoCategory.err_e0, a.get<std::string>("ohvar2_noCategory_errAlphaB_T"));

    parseHyperParameters(OHVarfinDParams2NoCategory.err_f0, a.get<std::string>("ohvar2_noCategory_errAlphaL_N"));
    parseHyperParameters(OHVarfinDParams2NoCategory.err_g0, a.get<std::string>("ohvar2_noCategory_errAlphaH_N"));
    parseHyperParameters(OHVarfinDParams2NoCategory.err_h0, a.get<std::string>("ohvar2_noCategory_errGammaEH_N"));
    parseHyperParameters(OHVarfinDParams2NoCategory.err_i0, a.get<std::string>("ohvar2_noCategory_errAlphaS_N"));
    parseHyperParameters(OHVarfinDParams2NoCategory.err_j0, a.get<std::string>("ohvar2_noCategory_errAlphaB_N"));
*/
    
}