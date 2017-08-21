//
//  PileUpUnit.h
//  OHVarFinder
//
//  Created by 森山卓也 on 2016/03/11.
//  Copyright (c) 2016年 森山卓也. All rights reserved.
//

#ifndef __OHVarFinder__PileUpUnit__
#define __OHVarFinder__PileUpUnit__

#include <vector>
#include <string>
#include <map>
#include <ctype.h>
#include "Utils.h"
#include <iostream>
#define INSDELBQ 100000

class PileUpUnit{
public:
    typedef enum { INS, DEL, SNP } Type;
    typedef enum { PLUS, MINUS } Strand;
    Type type;
    Strand strand;
    std::string obs;
    int bq;
    
    // If this unit is defined as INS or DEL, bq is INFINITE
    PileUpUnit(Type t, Strand s,char obs,int bq){
        this->type = t;   this->strand = s;
        this->obs  = "";  this->obs  += toupper(obs);
        this-> bq    = bq;
        if(t == INS || t == DEL) this->bq = INSDELBQ;
    }
    PileUpUnit(Type t, Strand s,char obs, char bq){
        this->type = t;   this->strand = s;
        this->obs  = "";  this->obs  += toupper(obs);
        this-> bq    = ((int)bq)-33;
        if(t == INS || t == DEL) this->bq = INSDELBQ;
    }
    PileUpUnit(Type t, Strand s,const std::string& obs, int bq){
        this->type = t;   this->strand = s;
        this->obs = obs;
        Utils::toUpperString(this->obs);
        this-> bq    = bq;
        if(t == INS || t == DEL) this->bq = INSDELBQ;
    }
    PileUpUnit(Type t, Strand s,const std::string& obs, char bq){
        this->type = t;   this->strand = s;
        this->obs = obs;
        Utils::toUpperString(this->obs);
        this-> bq    = ((int)bq)-33;
        if(t == INS || t == DEL) this->bq = INSDELBQ;
    }
    
    PileUpUnit(const PileUpUnit& p){
        this->type = p.type; this->strand = p.strand; this->obs = p.obs; this->bq = p.bq;
    }
    PileUpUnit & operator=(const PileUpUnit& p){
        this->type = p.type; this->strand = p.strand; this->obs = p.obs; this->bq = p.bq;
        return *this;
    }
    ~PileUpUnit(){};
    
    void printMySelf(){
        if(type == SNP) std::cout << "SNP: ";
        if(type == INS) std::cout << "INS: ";
        if(type == DEL) std::cout << "DEL: ";

        std::cout << obs << std::endl;
        
        if(strand == PLUS)  std::cout << "(+): ";
        if(strand == MINUS) std::cout << "(-): ";
        
    }
    
//    void printMyself(){
//        std::cerr << "(" << type << "," << strand << "," << obs << "," << bq << ")" ;
//    }

};
#endif /* defined(__OHVarFinder__PileUpUnit__) */
