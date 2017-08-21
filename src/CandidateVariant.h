//
//  CandidateVariant.h
//  OHVarFinder
//
//  Created by 森山卓也 on 2016/03/31.
//  Copyright (c) 2016年 森山卓也. All rights reserved.
//

#ifndef __OHVarFinder__CandidateVariant__
#define __OHVarFinder__CandidateVariant__

#include <stdio.h>
#include <string>
#include <iostream>
#include "PileUpUnit.h"
#include "VariantFromSAM.h"
#include "Variant.h"
#include <string.h>

//    int length, startInGenome, endInGenome;
//    std::string seq, originalString;
//    std::string ref, obs;
//    typedef enum { INS, DEL, SNP } Type;
//    Type type;
//    bool isIndel() const;
//    bool isSNP() const;

class CandidateVariant : public Variant {
public:
    class VariantDetail{
    public:
        VariantDetail(){
            isTriallilic = false;
            indelCoverCheck = false;
            isHeteroSNP = false;
            refNumT=-1;obsNumT=-1;refNumN=-1;obsNumN=-1;
            rateT=-1.0;rateN=-1.0;
        }
        void setNum(int ref, int obs, bool isTumor){
            if(isTumor){
                refNumT = ref; obsNumT = obs;
                rateT=0.0;
                if(refNumT+obsNumT) rateT = 1.0*obsNumT/(refNumT+obsNumT);
            }else{
                refNumN = ref; obsNumN = obs;
                rateN=0.0;
                if(refNumN+obsNumN) rateN = 1.0*obsNumN/(refNumN+obsNumN);
            }
        }
        int refNumT,obsNumT,refNumN,obsNumN;
        double rateT,rateN;
        bool isTriallilic;    // true if there is another obs base
        bool indelCoverCheck; // true if there is a indel nearby
        bool isHeteroSNP;     // true if this variant is heteroSNP
    };
    VariantDetail detail;
    std::string Chr;
    CandidateVariant();
    CandidateVariant(PileUpUnit pile,std::string _Chr, int _startInGenome, char _ref);
//    bool operator<(const CandidateVariant& rhs) const;
    int toPileUpStart(Variant::Type type, int start){
        if(type == Variant::SNP){
            return start+1;
        }else if(type == Variant::INS){
            return start;
        }else if(type == Variant::DEL){
            return start;
        }else{
            return -1;
        }
    };
 };


#endif /* defined(__Variant__) */


