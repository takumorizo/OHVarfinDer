//
//  PileUpUtils.h
//  OHVarFinder
//
//  Created by 森山卓也 on 2016/03/30.
//  Copyright (c) 2016年 森山卓也. All rights reserved.
//

#ifndef PileUpUtils_h
#define PileUpUtils_h



#include <stdio.h>
#include <vector>
#include <string>
#include <cmath>
#include <ctype.h>
#include <cmath>
class PileUpUtils {
public:
    // chr pos ref depthT baseT bqT depthN baseN bqN
    static inline void setMpileUpTN(std::vector<std::string> &cols, const char* pileup, const int buffSize, int colSize = 9){
        char buffer[buffSize];
        if(cols.size() != colSize) cols.resize(colSize);
        int colIdx = 0;
        int before = 0;
        int next = 0;
        
        int i = 0;
        while(pileup[i] != '\0'){
            if(pileup[i] == '\t' || pileup[i] == '\n'){
                before = next;
                next = i;
                if(colIdx > 0 && next-before-1 <= buffSize){
                    strncpy(buffer, pileup+before+1, next-before-1 );
                    buffer[next-before-1] = '\0';
                    cols[colIdx] = buffer;
                }else if(next-before <= buffSize){
                    strncpy(buffer, pileup+before, next-before);
                    buffer[next-before] = '\0';
                    cols[colIdx] = buffer;
                }else{
                    throw std::string("unexpected buffSize or colSize");
                }
                colIdx++;
            }
            i++;
        }
        if(colIdx != colSize){
            throw std::string("unexpected number of cols");
        }
    }
    // chr pos ref depthT baseT bqT depthN baseN bqN
    static inline void setMpileUpTN(std::vector<std::string> &cols,const std::string &pileup, int colSize = 9 ){
        if(cols.size() != colSize) cols.resize(colSize);
        int colIdx = 0;
        int before = 0;
        int next = 0;
        for(int i = 0 ; i < pileup.length(); i++){
            if(pileup[i] == '\t' || pileup[i] == '\n' ){
                before = next;
                next = i;
                if(colIdx > 0)  cols[colIdx] = pileup.substr(before+1,next-before-1); // ignore tab
                else            cols[colIdx] = pileup.substr(before,next-before);
                colIdx++;
            }
        }
        if(colIdx != colSize){
            throw std::string("unexpected number of cols");
        }
    }
};


#endif