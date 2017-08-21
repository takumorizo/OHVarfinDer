//
//  PileUpChecker.cpp
//  OHVarFinder
//
//  Created by 森山卓也 on 2016/03/14.
//  Copyright (c) 2016年 森山卓也. All rights reserved.
//

#include "PileUpChecker.h"
#include <stdio.h>
#include <vector>
#include <string>
#include <map>
#include <ctype.h>
#include "Utils.h"
#include <iostream>
#include "Parameters.h"
#include "PileUpUnit.h"
#include "PileUp.h"
#include "Variant.h"

PileUpChecker::PileUpChecker(Parameters::PileUpParameters param){
    this->param = param;
}


void PileUpChecker::setCandidateMutations(std::vector<CandidateVariant> &candidates, PileUp &pileUp,
                                          const std::string &Chr,    const std::string &position,
                                          const std::string &ref,    const std::string &depth,
                                          const std::string &bases,  const std::string &qualities){
    
    //    candidates.clear();
    //    std::vector<PileUpUnit>().swap(candidates);
    pileUp.setPileUpSummary(Chr, position, ref, depth, bases, qualities);
    
    //   === pre-checking candidate mutations ===
    int d = pileUp.depth;
    if ( !(param.minDepth <= d )){
        return ;
    }
    // LOG(logINFO) << " passed minDepth condition @setCandidateMutations" << std::endl;
    int refAll = pileUp.refNumP + pileUp.refNumM;
    int insAll = pileUp.insNum;
    int delAll = pileUp.delNum;
    if ( !(param.minObsNum <= (d - refAll + insAll + delAll)) ){
        return ;
    }
    // LOG(logINFO) << " passed minObs condition @setCandidateMutations" << std::endl;
    
    
    bool snvOK = (param.minObsNum <= (d-refAll));
    bool insOK = (param.minObsNum <= pileUp.insNum);
    bool delOK = (param.minObsNum <= pileUp.delNum);
    if( !(snvOK || insOK || delOK)) return;
    // LOG(logINFO) << " passed minObs condition in detail  @setCandidateMutations" << std::endl;

    pileUp.setPileUpSummaryInDetail(Chr, position, ref,depth, bases, qualities, param.minBQ);
    PileUpSummary summary = pileUp.seeDeteiledSummary();
    
    // === check SNV base ===
    char cs[] = {'A','T','G','C'};
    if( snvOK ){
        for(int i = 0 ; i < 4 ;i++)if(cs[i] != pileUp.ref){
            candidates.push_back( CandidateVariant( PileUpUnit(PileUpUnit::SNP,PileUpUnit::PLUS, cs[i], 1000000),Chr,std::atoi(position.c_str()),pileUp.ref ));
        }
    }
    // === check InDel ===
    if(insOK){
        std::map<std::string, int> allIns;
        for(std::map<std::string, int>::iterator it = summary.insPs.begin(); it != summary.insPs.end(); it++) {
            allIns[it->first] = 1;
        }
        for(std::map<std::string, int>::iterator it = summary.insMs.begin(); it != summary.insMs.end(); it++) {
            allIns[it->first] = 1;
        }
        for(std::map<std::string, int>::iterator it = allIns.begin(); it != allIns.end(); it++) {
            candidates.push_back( CandidateVariant( PileUpUnit(PileUpUnit::INS,PileUpUnit::PLUS, it->first, 1000000),Chr,std::atoi(position.c_str()),pileUp.ref));
        }
    }
    if(delOK){
        std::map<std::string, int> allDel;
        for(std::map<std::string, int>::iterator it = summary.delPs.begin(); it != summary.delPs.end(); it++) {
            allDel[it->first] = 1;
        }
        for(std::map<std::string, int>::iterator it = summary.delMs.begin(); it != summary.delMs.end(); it++) {
            allDel[it->first] = 1;
        }
        for(std::map<std::string, int>::iterator it = allDel.begin(); it != allDel.end(); it++) {
            candidates.push_back( CandidateVariant( PileUpUnit(PileUpUnit::DEL,PileUpUnit::PLUS, it->first, 1000000),Chr,std::atoi(position.c_str()),pileUp.ref));
        }
    }
}


void PileUpChecker::filterCandidates(std::vector<CandidateVariant> &candidates, const PileUp& pileUp, bool isTumor){

    // LOG(logINFO) << "param.minObsRate" << param.minObsRate << std::endl; 
    // LOG(logINFO) << "param.minObsRatePlus" << param.minObsRatePlus << std::endl; 
    // LOG(logINFO) << "param.minObsRateMinus" << param.minObsRateMinus << std::endl; 

    // LOG(logINFO) << "param.maxObsRate" << param.maxObsRate << std::endl; 
    // LOG(logINFO) << "param.maxObsRatePlus" << param.maxObsRatePlus << std::endl; 
    // LOG(logINFO) << "param.maxObsRateMinus" << param.maxObsRateMinus << std::endl; 

    // LOG(logINFO) << "param.minObsNum" << param.minObsNum << std::endl;       
    // LOG(logINFO) << "param.minObsNumPlus" << param.minObsNumPlus << std::endl;   
    // LOG(logINFO) << "param.minObsNumMinus" << param.minObsNumMinus << std::endl;  

    // LOG(logINFO) << "param.maxObsNum" << param.maxObsNum << std::endl;       
    // LOG(logINFO) << "param.maxObsNumPlus" << param.maxObsNumPlus << std::endl;   
    // LOG(logINFO) << "param.maxObsNumMinus" << param.maxObsNumMinus << std::endl;  

    // LOG(logINFO) << "param.minDepth : " << param.minDepth << std::endl;
    // LOG(logINFO) << "param.minDepthPlus : " << param.minDepthPlus << std::endl;
    // LOG(logINFO) << "param.minDepthMinus : " << param.minDepthMinus << std::endl;
    // LOG(logINFO) << "param.minRefNum : " << param.minRefNum << std::endl;
    // LOG(logINFO) << "param.minRefNumPlus : " << param.minRefNumPlus << std::endl;
    // LOG(logINFO) << "param.minRefNumMinus : " << param.minRefNumMinus << std::endl;

    // LOG(logINFO) << "param.maxDepth : " << param.maxDepth << std::endl;
    // LOG(logINFO) << "param.maxDepthPlus : " << param.maxDepthPlus << std::endl;
    // LOG(logINFO) << "param.maxDepthMinus : " << param.maxDepthMinus << std::endl;
    // LOG(logINFO) << "param.maxRefNum : " << param.maxRefNum << std::endl;
    // LOG(logINFO) << "param.maxRefNumPlus : " << param.maxRefNumPlus << std::endl;
    // LOG(logINFO) << "param.maxRefNumMinus : " << param.maxRefNumMinus << std::endl;


    std::vector<CandidateVariant> ansVec;
    PileUpSummary summary = pileUp.seeDeteiledSummary();
    
    if( !(summary.averageBaseQuality >= param.avgBaseQualityThreshold) ){
        // fail aveBaseQuality filter
        swap(candidates,ansVec);    return;
    }

    for(int i = 0 ; i < candidates.size() ; i++){
        // LOG(logINFO) << "start checking candidate @filterCandidates : " << candidates[i].getSymbol() << std::endl;
        // checking for depth and refBase numbers for SNV and indels
        int depth = 0;
        int depthPlus = 0;
        int depthMinus = 0;
        
        char cs[] = {'A','T','G','C','N'};
        for(int b = 0 ; b < 5; b++){
            std::map<char, int>::iterator itp = summary.basePs.find(cs[b]); if(itp != summary.basePs.end() ) depthPlus  += itp->second;
            std::map<char, int>::iterator itm = summary.baseMs.find(cs[b]); if(itm != summary.baseMs.end() ) depthMinus += itm->second;
        }
        for(std::map<std::string, int>::iterator it = summary.delPs.begin(); it != summary.delPs.end(); it++) depthPlus  += it->second;
        for(std::map<std::string, int>::iterator it = summary.delMs.begin(); it != summary.delMs.end(); it++) depthMinus += it->second;
        for(std::map<std::string, int>::iterator it = summary.insPs.begin(); it != summary.insPs.end(); it++) depthPlus  += it->second;
        for(std::map<std::string, int>::iterator it = summary.insMs.begin(); it != summary.insMs.end(); it++) depthMinus += it->second;
        
        depth = (depthPlus + depthMinus);
        
        int refNum = 0;
        int refNumPlus = 0;
        int refNumMinus = 0;
        std::map<char, int>::iterator itp = summary.basePs.find(pileUp.ref);   if(itp != summary.basePs.end())  refNumPlus  += itp->second;
        std::map<char, int>::iterator itm = summary.baseMs.find(pileUp.ref);   if(itm != summary.baseMs.end())  refNumMinus += itm->second;
        refNum = refNumPlus + refNumMinus;
        
        
        if (!(param.minDepth         <= depth        && depth        <= param.maxDepth        &&
              param.minDepthPlus     <= depthPlus    && depthPlus    <= param.maxDepthPlus    &&
              param.minDepthMinus    <= depthMinus   && depthMinus   <= param.maxDepthMinus   &&
              param.minRefNum        <= refNum       && refNum       <= param.maxRefNum       &&
              param.minRefNumPlus    <= refNumPlus   && refNumPlus   <= param.maxRefNumPlus   &&
              param.minRefNumMinus   <= refNumMinus  && refNumMinus  <= param.maxRefNumMinus ) ){
            // LOG(logINFO) << "filter candidate, fail depth, ref condition @filterCandidates " << std::endl;
            continue;
        }
        // for specific SNV, indels
        // see obs Num
        int    obsNum       = 0;
        int    obsNumPlus   = 0;
        int    obsNumMinus  = 0;
        // see obs Rate
        double obsRate      = 0.0;
        double obsRatePlus  = 0.0;
        double obsRateMinus = 0.0;
        
        if(candidates[i].type == Variant::SNP){
            itp = summary.basePs.find(candidates[i].obs[0]); if(itp != summary.basePs.end()) obsNumPlus  += itp->second;
            itm = summary.baseMs.find(candidates[i].obs[0]); if(itm != summary.baseMs.end()) obsNumMinus += itm->second;
            obsNum += obsNumPlus + obsNumMinus;
        
        }else if(candidates[i].type == Variant::INS){
            std::map<std::string, int>::iterator itIp = summary.insPs.find(candidates[i].obs); if(itIp != summary.insPs.end()) obsNumPlus  += itIp->second;
            std::map<std::string, int>::iterator itIm = summary.insMs.find(candidates[i].obs); if(itIm != summary.insMs.end()) obsNumMinus += itIm->second;
            obsNum += obsNumPlus + obsNumMinus;
        
        }else if(candidates[i].type == Variant::DEL){
            std::map<std::string, int>::iterator itDp = summary.delPs.find(candidates[i].ref); if(itDp != summary.delPs.end()) obsNumPlus  += itDp->second;
            std::map<std::string, int>::iterator itDm = summary.delMs.find(candidates[i].ref); if(itDm != summary.delMs.end()) obsNumMinus += itDm->second;
            obsNum += obsNumPlus + obsNumMinus;

        }
        
        // obsNum condition
        if(!(param.minObsNum      <= obsNum      && obsNum      <= param.maxObsNum     &&
             param.minObsNumPlus  <= obsNumPlus  && obsNumPlus  <= param.maxObsNumPlus &&
             param.minObsNumMinus <= obsNumMinus && obsNumMinus <= param.maxObsNumMinus)){
            // LOG(logINFO) << "filter candidate, fail obsNum condition @filterCandidates " << std::endl;
            continue;
        }
        // obs Rate condition
        if(depth > 0){
            obsRate = (obsNum * 1.0 ) / depth;
        }
        if(depthPlus > 0 ){
            obsRatePlus = (obsNumPlus * 1.0) /depthPlus;
        }
        if(depthMinus > 0 ){
            obsRateMinus = (obsNumMinus * 1.0) /depthMinus;
        }
        if (!(param.minObsRate      <= obsRate      && obsRate      <= param.maxObsRate     &&
              param.minObsRatePlus  <= obsRatePlus  && obsRatePlus  <= param.maxObsRatePlus &&
              param.minObsRateMinus <= obsRateMinus && obsRateMinus <= param.maxObsRateMinus)){
            // LOG(logINFO) << "filter candidate, fail obsRate condition @filterCandidates " << std::endl;
            // LOG(logINFO) << "minObsRate : " << param.minObsRate << std::endl;
            // LOG(logINFO) << "obsRate : " << obsRate << std::endl;
            // LOG(logINFO) << "maxObsRate : " << param.maxObsRate << std::endl;

            // LOG(logINFO) << "minObsRate(+) : " << param.minObsRatePlus << std::endl;
            // LOG(logINFO) << "obsRate(+) : " << obsRatePlus << std::endl;
            // LOG(logINFO) << "maxObsRate(+) : " << param.maxObsRatePlus << std::endl;

            // LOG(logINFO) << "minObsRate(-) : " << param.minObsRateMinus<< std::endl;
            // LOG(logINFO) << "obsRate(-) : " << obsRateMinus << std::endl;
            // LOG(logINFO) << "maxObsRate(-) : " << param.maxObsRateMinus << std::endl;

            continue;
        }
        // passed all filter
        candidates[i].detail.setNum(refNum, obsNum, isTumor);
        ansVec.push_back(candidates[i]);
    }
    swap(candidates,ansVec);
}


