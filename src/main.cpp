/*
 * HapMuC
 * http://github.com/usuyama/hapmuc
 *
 * Copyright 2013, Naoto Usuyama
 */

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include "PileUp.h"
#include <boost/lexical_cast.hpp>
#include "Parameters.h"
#include "PileUpChecker.h"
//#include <ext/stdio_filebuf.h>
#include <cstdlib>
#include "PileUpUtils.h"
#include <list>
#include "CandidateVariant.h"
#include "CandidateWindow.h"
#include "VariantUtils.h"
#include "OHVarfinDer.h"
#include "BaseMutationCaller.h"
#include "PileUpManager.h"
#include "MutationCaller.h"

// SAMTOOLS is predefined in the compiling
int main(int argc, const char *argv[]){
    Parameters parameters = Parameters(argc, argv);
    //  
    // Open log file
    //
    FILE* logFile = fopen((parameters.outFilePrefix + ".log").c_str(), "w");
    Output2FILE::Stream() = logFile;
    LOG(logINFO) << "OHVarfinDer ver1.1" << std::endl;
    //
    // Search Candidate Mutation and HeteroSNP
    //
    PileUpManager pileManager = PileUpManager();
    std::vector<CandidateWindow>  allWindows = pileManager.searchCandidateWindow(parameters);
    
    //
    // Open reference genome
    //
    faidx_t *fai = NULL;
    fai = fai_load(parameters.refFileName.c_str());
    if (!fai) {
        LOG(logERROR) << "Cannot open reference sequence file." << std::endl;
        exit(1);
    }
    
    //
    // Prepare bam files
    //
    MyBam tumorBam = MyBam(parameters.tumorBam);
    MyBam normalBam = MyBam(parameters.normalBam);
    BamReader tumorBamReader = BamReader(&tumorBam, fai,   parameters, parameters.maxReads);
    BamReader normalBamReader = BamReader(&normalBam, fai, parameters, parameters.maxReads);
    
    // output Files
    std::ofstream outStream((parameters.outFilePrefix + ".calls.txt").c_str());
    
    if(parameters.method == Parameters::OHVARFINDER){
        outStream << MutationCallResult::getHeader() << std::endl;
        LOG(logINFO) << " start Ohvarfinder2 algorithm" << std::endl;
        OHVarfinDer caller = OHVarfinDer(tumorBamReader, normalBamReader,parameters,fai);
        
        for(int i = 0; i < allWindows.size(); i++){
            try {
                MutationCallResult result = caller.call(allWindows[i]);
                outStream << result.getOutput() << std::endl;
                LOG(logDEBUG) << result.getOutput() << std::endl;
                std::cout << result.getOutput() << std::endl;
            } catch (std::string &s) {
                LOG(logERROR) << "something unexpected happened. exit." << std::endl;
                LOG(logERROR) << s << std::endl;
                exit(1);
            }
        }
    }else if(parameters.method == Parameters::HAPMUC){
        outStream << MutationCallResult::getHeader() << std::endl;
        LOG(logINFO) << " start HapMuC algorithm" << std::endl;
        MutationCaller caller = MutationCaller(tumorBamReader, normalBamReader, parameters, fai);
        for(int i = 0; i < allWindows.size(); i++){
            try {
                MutationCallResult result = caller.call(allWindows[i]);
                outStream << result.getOutput() << std::endl;
                LOG(logDEBUG) << result.getOutput() << std::endl;
                std::cout << result.getOutput() << std::endl;
            } catch (std::string &s) {
                LOG(logERROR) << "something unexpected happened. exit." << std::endl;
                LOG(logERROR) << s << std::endl;
                exit(1);
            }
        }
    }else if(parameters.method == Parameters::HETEROSNPCALL){
        outStream << MutationCallResult::getHeader() << std::endl;
        LOG(logINFO) << " start HETEROSNPCALL algorithm" << std::endl;
        for(int i = 0; i < allWindows.size(); i++){
            try {
                outStream << allWindows[i].toString();
            } catch (std::string &s) {
                LOG(logERROR) << "something unexpected happened. exit." << std::endl;
                LOG(logERROR) << s << std::endl;
                exit(1);
            }
        }
    }else{
        LOG(logERROR) << "unknown algorithm specified." << std::endl;
        exit(1);
    }
    outStream.close();
    fai_destroy(fai);
}


