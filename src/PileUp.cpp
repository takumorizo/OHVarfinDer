//
//  PileUp.cpp
//  OHVarFinder
//
//  Created by 森山卓也 on 2016/03/11.
//  Copyright (c) 2016年 森山卓也. All rights reserved.
//

#include "PileUp.h"
#include <vector>
#include <string>
#include <map>
#include "Utils.h"
#include "PileUpUnit.h"
#include "log.h"
#include <boost/lexical_cast.hpp>
#include <iostream>
void PileUp::setPileUpSummary( const std::string& Chr,    const std::string& position,
                               const std::string& ref,    const std::string& depth,
                               const std::string& bases,  const std::string& qualities){
    this->Chr = Chr;
    this->pos = boost::lexical_cast<long long>(position);
    this->ref = toupper(ref[0]);
    this->depth = boost::lexical_cast<int>(depth);
    this->bases = bases;

    this->qualities.clear();
    std::vector<int>().swap(this->qualities);
    for(int i = 0 ; i < qualities.length(); i++){
        this->qualities.push_back((int)qualities[i]-33);
    }
    this->refNumP=Utils::countChar(this->bases, '.');
    this->refNumM=Utils::countChar(this->bases, ',');
    
    this->insNum=Utils::countChar(this->bases, '+');
    this->delNum=Utils::countChar(this->bases, '-');

}


void PileUp::setPileUpSummaryInDetail(const std::string& Chr,   const std::string& position,
                                      const std::string& ref,   const std::string& depth,
                                      const std::string& bases, const std::string& qualities,
                                      int minBQ){
    // make detailed summary counts
    this->correspondBaseAndBQ(Chr, position, ref, depth, bases, qualities);
    
//     count SNV,INS,DEL, considering minBQ
    this->summary.init();
    double totalBaseQuality=0.0;
    for(std::vector<PileUpUnit>::iterator it = this->pileUpUnits.begin(); it != this->pileUpUnits.end(); it++){
        // check Base quality
        totalBaseQuality += it->bq;

        // count Base/In/Del : sutisfy the base quality condition.
        if(it->bq < minBQ) continue;
        if(it->type == PileUpUnit::SNP){
            if(it->strand == PileUpUnit::PLUS ){
                if(this->summary.basePs.count(it->obs[0]) == 0) this->summary.basePs[it->obs[0]] = 1;
                else                                            this->summary.basePs[it->obs[0]]++;
            }
            if(it->strand == PileUpUnit::MINUS){
                if(this->summary.baseMs.count(it->obs[0]) == 0) this->summary.baseMs[it->obs[0]] = 1;
                else                                            this->summary.baseMs[it->obs[0]]++;
            }
            
        }else if(it->type == PileUpUnit::INS ) {
            if(it->strand == PileUpUnit::PLUS ){
                if(this->summary.insPs.count(it->obs) == 0) this->summary.insPs[it->obs] = 1;
                else                                        this->summary.insPs[it->obs]++;
            }
            if(it->strand == PileUpUnit::MINUS){
                if(this->summary.insMs.count(it->obs) == 0) this->summary.insMs[it->obs] = 1;
                else                                        this->summary.insPs[it->obs]++;
            }
            
        }else if(it->type == PileUpUnit::DEL ) {
            if(it->strand == PileUpUnit::PLUS ){
                if(this->summary.delPs.count(it->obs) == 0) this->summary.delPs[it->obs] = 1;
                else                                        this->summary.delPs[it->obs]++;
            }
            if(it->strand == PileUpUnit::MINUS){
                if(this->summary.delMs.count(it->obs) == 0) this->summary.delMs[it->obs] = 1;
                else                                        this->summary.delMs[it->obs]++;
            }
        }
    }
    this->summary.averageBaseQuality = (totalBaseQuality*1.0) / (1.0 * pileUpUnits.size());
}


void PileUp::correspondBaseAndBQ(const std::string& Chr,   const std::string& position,
                                 const std::string& ref,   const std::string& depth,
                                 const std::string& bases, const std::string& qualities){
    this->pileUpUnits.clear();
    std::vector<PileUpUnit>().swap(this->pileUpUnits);

    int b = 0; int q = 0 ;  // index for the base and basequality
    this->setPileUpSummary(Chr,position,ref,depth,bases,qualities);
    std::string tempStr("");
    while( b < this->bases.length() ){
        if (this->bases[b] == '.' ) { // match +
            this->pileUpUnits.push_back( PileUpUnit(PileUpUnit::SNP, PileUpUnit::PLUS, this->ref, this->qualities[q]) );
            b++;
            q++;
        }else if(this->bases[b] == ','){ //match -
            this->pileUpUnits.push_back( PileUpUnit(PileUpUnit::SNP, PileUpUnit::MINUS, this->ref, this->qualities[q]) );
            b++;
            q++;
        }else if(this->bases[b] == 'A' || this->bases[b] == 'C' || this->bases[b] == 'G' || this->bases[b] == 'T'|| this->bases[b] == 'N' ){
            // SNP +
            this->pileUpUnits.push_back( PileUpUnit(PileUpUnit::SNP, PileUpUnit::PLUS, this->bases[b], this->qualities[q]) );
            b++;
            q++;
        }else if(this->bases[b] == 'a' || this->bases[b] == 'c' || this->bases[b] == 'g' || this->bases[b] == 't' || this->bases[b] == 'n' ){
            // SNP -
            this->pileUpUnits.push_back( PileUpUnit(PileUpUnit::SNP, PileUpUnit::MINUS, toupper(this->bases[b]), this->qualities[q]) );
            b++;
            q++;
        }else if(this->bases[b]=='$'){ // read Start
            b++;
        }else if(this->bases[b]=='^'){ // skip MAPQ score
            b += 2;
        }else if(this->bases[b]=='*'){ // deletion
            b++;
            q++;
        }else if(this->bases[b] == '>' || this->bases[b] == '<'){ // reference skip
            b++;
            q++;
        }else if(this->bases[b] == '+' || this->bases[b] == '-'){ // indel starting
//            LOG(logERROR) << "with indel bases line : " << this->bases << std::endl;
            tempStr = "";
            int tmp_b = b+1;
            while('0' <= this->bases[tmp_b] && this->bases[tmp_b] <= '9'){
                tempStr += this->bases[tmp_b];
                tmp_b++;
            }
            int insDelLen = boost::lexical_cast<int>(tempStr);
            std::string indel = bases.substr(tmp_b,insDelLen);
//            LOG(logERROR) << "indel : " << indel << std::endl;
            PileUpUnit::Strand direction;
            if( toupper(indel[0]) == indel[0] ){
                direction = PileUpUnit::PLUS;
            }else{
                direction = PileUpUnit::MINUS;
            }
            
            Utils::toUpperString(indel);
            if(this->bases[b]=='+'){
                this->pileUpUnits.push_back( PileUpUnit(PileUpUnit::INS, direction, indel, INSDELBQ));
            }else{
                this->pileUpUnits.push_back( PileUpUnit(PileUpUnit::DEL, direction, indel, INSDELBQ));
            }
            b = tmp_b + insDelLen;
        }else{
            LOG(logERROR) << "pileup format error !!!" << std::endl;
            throw std::string("pileup_error");
        }
    }
}
// return copy of the detailed summary information
PileUpSummary  PileUp::seeDeteiledSummary()const{
    return this->summary;
}
