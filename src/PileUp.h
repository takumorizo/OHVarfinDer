//
//  PileUp.h
//  OHVarFinder
//
//  Created by 森山卓也 on 2016/03/11.
//  Copyright (c) 2016年 森山卓也. All rights reserved.
//

#ifndef __OHVarFinder__PileUp__
#define __OHVarFinder__PileUp__


#include <vector>
#include <string>
#include <map>
#include "PileUpUnit.h"
#include <cstdio>
class PileUpSummary{
public:
    PileUpSummary(){
        insNum=0;delNum=0;depth=0;
    }
    PileUpSummary(const PileUpSummary& p){
        this->basePs = p.basePs; this->baseMs = p.baseMs;
        this->insNum = p.insNum; this->delNum = p.delNum; this->depth = p.depth;
        this->insPs = p.insPs;   this->delPs = p.delPs;   this->insMs = p.insMs; this->delMs = p.delMs;
        this->averageBaseQuality = p.averageBaseQuality;
    }
    PileUpSummary& operator=(const PileUpSummary& p){
        this->basePs = p.basePs; this->baseMs = p.baseMs;
        this->insNum = p.insNum; this->delNum = p.delNum; this->depth = p.depth;
        this->insPs = p.insPs;   this->delPs = p.delPs;   this->insMs = p.insMs; this->delMs = p.delMs;
        this->averageBaseQuality = p.averageBaseQuality;
        return *this;
    }
    void init(){
        basePs.clear(); baseMs.clear();
        insNum=0;delNum=0;
        insPs.clear(); delPs.clear();
        insMs.clear(); delMs.clear();
    }
    std::map<char,int> basePs;
    std::map<char,int> baseMs;
    
    int insNum;
    int delNum;
    int depth;
    double averageBaseQuality;
    std::map<std::string,int> insPs;
    std::map<std::string,int> delPs;
    std::map<std::string,int> insMs;
    std::map<std::string,int> delMs;
};

class PileUp {
public:
    PileUp(){
        refNumP=0; refNumM=0;
        insNum=0;delNum=0;
        pos=-1;ref='-';depth=0;
    };
    // simple Summary
    int refNumP;
    int refNumM;
    int insNum;
    int delNum;
    // all Info
    std::vector<PileUpUnit> pileUpUnits;
    
    // pileUpLine basic info
    std::string Chr;
    long long pos;
    char ref;
    int depth;
    std::string bases;
    std::vector<int> qualities;

    void setPileUpSummary(const  std::string& Chr,    const std::string& position,
                          const  std::string& ref,    const std::string& depth,
                          const  std::string& bases,  const std::string& qualities);
    void setPileUpSummaryInDetail(const std::string& Chr,    const std::string& position,
                                  const std::string& ref,    const std::string& depth,
                                  const std::string& bases,  const std::string& qualities,int minBQ);
    PileUpSummary seeDeteiledSummary()const;
//    void printPileUpUnits();

private:
    void correspondBaseAndBQ( const std::string& Chr,    const std::string& position,
                              const std::string& ref,    const std::string& depth,
                              const std::string& bases,  const std::string& qualities);
    PileUpSummary summary; // detailed summary
};


#endif /* defined(__OHVarFinder__PileUp__) */

