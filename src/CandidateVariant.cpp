//
//  CandidateVariant.cpp
//  OHVarFinder
//
//  Created by 森山卓也 on 2016/03/31.
//  Copyright (c) 2016年 森山卓也. All rights reserved.
//

#include "CandidateVariant.h"


CandidateVariant::CandidateVariant(){
    Chr = "NA";
}

// _startInGenome is 1-based index, and cnvert to 0-based bed format
CandidateVariant::CandidateVariant(PileUpUnit pile,std::string _Chr, int _startInGenome, char _ref){
    Chr = _Chr;
    if(pile.type == PileUpUnit::SNP){
        type = Variant::SNP;
        ref = "";
        ref += _ref;
        obs = pile.obs;
        length = 1;
        
        startInGenome = _startInGenome - 1; // convert to 0-based index
        endInGenome   = startInGenome +1;
        seq = "";
        seq += ref; seq += "=>"; seq += obs;
        originalString = seq;
    }else if(pile.type == PileUpUnit::INS){
        type = Variant::INS;
        ref = "-";
        obs = pile.obs;
        length = (int)obs.length();

        startInGenome = _startInGenome; // convert to 0-based index for insertion start
        endInGenome   = startInGenome;
        seq = obs;
        originalString = "+";
        originalString += obs;
    }else if(pile.type == PileUpUnit::DEL){
        type = Variant::DEL;
        ref = pile.obs;
        obs = "-";
        length = (int)ref.length();
        
        startInGenome = _startInGenome; // convert to 0-based index for deletion start
        endInGenome   = startInGenome + length;
        seq = ref;
        originalString = "-";
        originalString += ref;
    }
}

//bool CandidateVariant::operator<(const CandidateVariant& rhs) const{
//    if(this->Chr != rhs.Chr){
//        return strcmp(this->Chr.c_str(), rhs.Chr.c_str());
//    }else if(this->startInGenome != rhs.startInGenome){
//        return this->startInGenome < rhs.startInGenome;
//    }else if(this->endInGenome != rhs.endInGenome) {
//        return this->endInGenome < rhs.endInGenome;
//    }else{
//        return false;
//    }
//}
