//
//  PileUpManager.cpp
//  OHVarFinder
//
//  Created by 森山卓也 on 2016/05/24.
//  Copyright © 2016年 森山卓也. All rights reserved.
//

#include "PileUpManager.h"
#include <stdio.h>
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
#include "log.h"

PileUpManager::PileUpManager(){

}


std::vector<CandidateWindow> PileUpManager::searchCandidateWindow(Parameters parameters){

    PileUpChecker pileCheckerHeteroSNP = PileUpChecker(parameters.pileHetero);
    PileUpChecker pileCheckerTriAllele = PileUpChecker(parameters.pileTriAllele);
    PileUpChecker pileCheckerN         = PileUpChecker(parameters.pileN);
    PileUpChecker pileCheckerT         = PileUpChecker(parameters.pileT);
    
    PileUp pileN = PileUp();
    PileUp pileT = PileUp();
    
    std::vector<CandidateVariant> tempCandidates;
    std::vector<CandidateVariant> allVariants;
    std::vector<CandidateWindow>  allWindows;

    FILE *pileUpStdOut;
    bool existFile = true;
    pileUpStdOut = fopen(parameters.pileupFile.c_str(), "r");
    LOG(logINFO) << "pileupFile : " << parameters.pileupFile << std::endl;
    if( parameters.pileupFile == "" || pileUpStdOut == NULL){
        // call subprocess and exec pileUp
        LOG(logINFO) << "samtools : " << SAMTOOLS << std::endl;
        std::string command("");
        command += SAMTOOLS;
        command += " mpileup -q 15 -BQ 0";
        command += " -f ";
        command += parameters.refFileName;
        if( parameters.region != ""){
            command += " -r ";
            command += parameters.region;
        }
        command += " ";
        command += parameters.tumorBam;
        command += " ";
        command += parameters.normalBam;
        LOG(logINFO) << "Pile up command : " << command << std::endl;
        
        existFile = false;
        if(!(pileUpStdOut=popen(command.c_str(),"r"))){
            LOG(logERROR) << " popen failed when doing pileup " << std::endl;
            throw std::string(" failed open process for mpileup ");
        }
        
    }
    
    std::string line;
    std::vector<std::string> cols(9);
    long long temp = 0;
    const int pileUpBufferSize = parameters.pileUpBufferSize;
    char readline[pileUpBufferSize];

    while ( std::fgets(readline, pileUpBufferSize, pileUpStdOut) != NULL ) {
        PileUpUtils::setMpileUpTN(cols, readline, pileUpBufferSize);
        
        int depthT = std::atoi(cols[3].c_str());
        int depthN = std::atoi(cols[6].c_str());
        
        if( depthT >= parameters.pileT.minDepth && depthN >= parameters.pileN.minDepth){            
            // init tempCandidates
            if(tempCandidates.size() > 0) tempCandidates.clear();
            
            if(isHeteroSNP(cols, pileN, tempCandidates, parameters, pileCheckerHeteroSNP)){
                VariantUtils::addToGlobalCandidates(allVariants, tempCandidates);
                continue;
            }
            if(isTriallelic(cols, pileT, tempCandidates, parameters, pileCheckerTriAllele) ){
                continue;
            }
            if(isSomaticCandidate(cols, pileN, pileT, tempCandidates, parameters, pileCheckerN, pileCheckerT)){
                VariantUtils::addToGlobalCandidates(allVariants, tempCandidates);
                continue;
            }
        }
    }
    LOG(logINFO) << " allVariants.size :" << allVariants.size() << std::endl;
    // make Candidate Windows
    int maxWindow = (parameters.maxInsertSize + parameters.maxReadLength);
    VariantUtils::makeWindowsFromSortedCandidates(allVariants, allWindows, parameters.heteroSNPConfidenceInterval, parameters.minDistanceGermlineIndel, maxWindow);
    LOG(logINFO) << " size of all windows  :" << allWindows.size() << std::endl;

    
    for(int i = 0 ; i < allWindows.size(); i++){
        LOG(logINFO) << allWindows[i].toString() << std::endl;
    }
    
    if(existFile){
        std::cout << "pileup File is closed" << std::endl;
        fclose(pileUpStdOut);
    }else {
        std::cout<< "pileup process is closed" << std::endl;
        pclose(pileUpStdOut);
    }

    return allWindows;
}

bool PileUpManager::isHeteroSNP(const std::vector<std::string> &cols, PileUp &pileN,
                                std::vector<CandidateVariant> &tempCandidates,
                                const Parameters &param, PileUpChecker &pileChecker){
    tempCandidates.clear();
    pileChecker.setCandidateMutations(tempCandidates,pileN,cols[0],cols[1],cols[2], cols[6], cols[7], cols[8]);
    if(tempCandidates.size() > 0){
        pileN.setPileUpSummaryInDetail(cols[0],cols[1],cols[2],cols[6],cols[7],cols[8], param.pileN.minBQ);
        pileChecker.filterCandidates(tempCandidates, pileN, false);
        VariantUtils::removeFalseHeteroSNP(tempCandidates, param.heteroSNPConfidenceInterval);
    }
    return (tempCandidates.size() > 0);
}


bool PileUpManager::isTriallelic(const std::vector<std::string> &cols, PileUp &pileT,
                                 std::vector<CandidateVariant> &tempCandidates, 
                                 const Parameters &param, PileUpChecker &pileChecker){
    tempCandidates.clear();
    pileChecker.setCandidateMutations(tempCandidates,pileT,cols[0],cols[1],cols[2],cols[3],cols[4],cols[5]);
    if(tempCandidates.size() == 0 ) return false; 
    pileT.setPileUpSummaryInDetail(cols[0],cols[1],cols[2],cols[3],cols[4],cols[5], param.pileT.minBQ);
    pileChecker.filterCandidates(tempCandidates, pileT, true);

    if(tempCandidates.size() > 1 ) return true;
    else                           return false;
}

bool PileUpManager::isSomaticCandidate(const std::vector<std::string> &cols, PileUp &pileN, PileUp &pileT,
                                       std::vector<CandidateVariant> &tempCandidates, 
                                       const Parameters &param, PileUpChecker &pileCheckerN, PileUpChecker &pileCheckerT){
    tempCandidates.clear();
    pileT.setPileUpSummaryInDetail(cols[0],cols[1],cols[2],cols[3],cols[4],cols[5], param.pileT.minBQ);
    pileCheckerT.setCandidateMutations(tempCandidates,pileT,cols[0],cols[1],cols[2],cols[3],cols[4],cols[5]);
    pileCheckerT.filterCandidates(tempCandidates, pileT, true);
    if(tempCandidates.size() == 0 ) return false;

    pileN.setPileUpSummary(cols[0],cols[1],cols[2], cols[6], cols[7], cols[8]);
    pileN.setPileUpSummaryInDetail(cols[0],cols[1],cols[2], cols[6], cols[7], cols[8], param.pileN.minBQ);
    pileCheckerN.filterCandidates(tempCandidates, pileN, false);
    if(tempCandidates.size() == 0 ) return false;
    else                            return true;
}
