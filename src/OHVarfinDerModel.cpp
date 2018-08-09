//
//  OHVarfinder2Model.cpp
//  OHVarFinder
//
//  Created by 森山卓也 on 2016/07/06.
//  Copyright © 2016年 森山卓也. All rights reserved.
//

#include "OHVarfinDerModel.h"
#include <vector>
#include <cstdlib>
#include "ModelUtils.h"
#include "log.h"
#include <math.h>

// TODO : remove comments
OHVarfinDerModel::OHVarfinDerModel(const std::vector<std::vector<double> > &tumorLiks_, const std::vector<std::vector<double> > &normalLiks_,
                                   const std::vector<std::vector<int> > &tumorIndicator_, const std::vector<std::vector<int> > &normalIndicator_,
                                   const Parameters::BayesEMParameters &prior_,
                                   int updateCount_ , double minLqConvergence_, int thresModelChangeDepth_):
tumorLiks(tumorLiks_), normalLiks(normalLiks_), tumorIndicator(tumorIndicator_), normalIndicator(normalIndicator_),
prior(prior_),updateCount(updateCount_), minLqConvergence(minLqConvergence_), thresModelChangeDepth(thresModelChangeDepth_),
ln2(log(2.0)) {
    if (tumorLiks.size() == 0) {
        LOG(logINFO) << "tumorLiks.size() == 0 " << std::endl;
        throw std::string("tumorLiks.size() == 0");
    }
    if (normalLiks.size() == 0) {
        LOG(logINFO) << "normalLiks.size() == 0 " << std::endl;
        throw std::string("normalLiks.size() == 0");
    }
    if (tumorLiks[0].size()  != 22 ){
        LOG(logINFO) << "tumorLiks[0].size() != 22" << std::endl;
        throw std::string("tumorLiks[0].size() != 22");
    }
    if (normalLiks[0].size() != 22 ){
        LOG(logINFO) << "tumorLiks[0].size() != 22" << std::endl;
        throw std::string("normalLiks[0].size() != 22");
    }
}

void OHVarfinDerModel::updateParams(std::vector<double>& param, const std::vector<std::vector<int> > &rule,
                                     const std::vector<double> &EZ, const std::vector<int> &allIndicatorInLine,
                                     int readNum){
    for(int j = 0; j < rule.size(); j++ ){
        for(int i = 0; i < readNum ; i++) for(int k = 0 ; k < rule[j].size(); k++){
            int idx = rule[j][k];
            if(allIndicatorInLine[at(i,idx)] != 0 ){
                param[j] += EZ[at(i,idx)];
            }
        }
    }
}

std::vector<double> OHVarfinDerModel::evaluateEqContentsTumor(
                                                int readsNum, const std::vector<double> &allLiksInLine, const std::vector<int> &allIndicatorInLine,
                                                const std::vector<double> &digPIF,  const std::vector<double> &digPIH,
                                                const std::vector<double> &digEPSL, const std::vector<double> &digEPSH, const std::vector<double> &digEPSB ){
    std::vector<double> EqContents(22*readsNum,0.0);
    for(int i = 0; i < readsNum ; i++){
        if( allIndicatorInLine[at(i,0)]  != 0 ) EqContents[at(i,0)]  =  (  digPIF[0] + digEPSL[1]             + allLiksInLine[at(i,0)] ) ;
        if( allIndicatorInLine[at(i,1)]  != 0 ) EqContents[at(i,1)]  =  (  digPIF[1]                          + allLiksInLine[at(i,1)] ) ;
        if( allIndicatorInLine[at(i,2)]  != 0 ) EqContents[at(i,2)]  =  (  digEPSB[0]+ digPIF[0] + digEPSL[0] + allLiksInLine[at(i,2)] ) ;
        if( allIndicatorInLine[at(i,3)]  != 0 ) EqContents[at(i,3)]  =  (  digEPSB[1]+ digPIF[0] + digEPSL[0] + allLiksInLine[at(i,3)] ) ;

        if( allIndicatorInLine[at(i,4)]  != 0 ) EqContents[at(i,4)]  =  (  digPIH[0]              + allLiksInLine[at(i,4)] ) ;
        if( allIndicatorInLine[at(i,5)]  != 0 ) EqContents[at(i,5)]  =  (  digPIH[1] + digEPSH[1] + allLiksInLine[at(i,5)] ) ;
        if( allIndicatorInLine[at(i,6)]  != 0 ) EqContents[at(i,6)]  =  (  digPIH[2]              + allLiksInLine[at(i,6)] ) ;
        if( allIndicatorInLine[at(i,7)]  != 0 ) EqContents[at(i,7)]  =  (  digPIH[1] + digEPSH[0] + allLiksInLine[at(i,7)] ) ;

        if( allIndicatorInLine[at(i,8)]  != 0 ) EqContents[at(i,8)]  =  (  digPIH[0]  +              digEPSL[1]              + allLiksInLine[at(i,8)] ) ;
        if( allIndicatorInLine[at(i,9)]  != 0 ) EqContents[at(i,9)]  =  (  digPIH[1]  + digEPSL[1] + digEPSL[1]              + allLiksInLine[at(i,9)] ) ;
        if( allIndicatorInLine[at(i,10)] != 0 ) EqContents[at(i,10)] =  (  digPIH[2]                                         + allLiksInLine[at(i,10)] ) ;
        if( allIndicatorInLine[at(i,11)] != 0 ) EqContents[at(i,11)] =  (  digPIH[1]  + digEPSL[0] + digEPSL[0]              + allLiksInLine[at(i,11)] ) ;
        if( allIndicatorInLine[at(i,12)] != 0 ) EqContents[at(i,12)] =  (  digEPSB[0] + digPIH[0]  + digEPSL[0]              + allLiksInLine[at(i,12)] ) ;
        if( allIndicatorInLine[at(i,13)] != 0 ) EqContents[at(i,13)] =  (  ln2 + digEPSB[0] + digPIH[1]  + digEPSL[1] + digEPSL[0] + allLiksInLine[at(i,13)] ) ;
        if( allIndicatorInLine[at(i,14)] != 0 ) EqContents[at(i,14)] =  (  digEPSB[1] + digPIH[0]  + digEPSL[0]              + allLiksInLine[at(i,14)] ) ;
        if( allIndicatorInLine[at(i,15)] != 0 ) EqContents[at(i,15)] =  (  ln2 + digEPSB[1] + digPIH[1]  + digEPSL[0] + digEPSL[1] + allLiksInLine[at(i,15)] ) ;

        if( allIndicatorInLine[at(i,16)] != 0 ) EqContents[at(i,16)] =  (  digPIF[0]              + allLiksInLine[at(i,16)] ) ;
        if( allIndicatorInLine[at(i,17)] != 0 ) EqContents[at(i,17)] =  (  digEPSB[1] + digPIF[1] + allLiksInLine[at(i,17)] ) ;
        if( allIndicatorInLine[at(i,18)] != 0 ) EqContents[at(i,18)] =  (  digEPSB[0] + digPIF[1] + allLiksInLine[at(i,18)] ) ;

        if( allIndicatorInLine[at(i,19)] != 0 ) EqContents[at(i,19)] =  (  digPIF[0]              + allLiksInLine[at(i,19)] ) ;
        if( allIndicatorInLine[at(i,20)] != 0 ) EqContents[at(i,20)] =  (  digEPSB[0] + digPIF[1] + allLiksInLine[at(i,20)] ) ;
        if( allIndicatorInLine[at(i,21)] != 0 ) EqContents[at(i,21)] =  (  digEPSB[1] + digPIF[1] + allLiksInLine[at(i,21)] ) ;
    }
    return EqContents;
}

std::vector<double> OHVarfinDerModel::evaluateEqContentsError(
                                                int readsNum, const std::vector<double> &allLiksInLine, const std::vector<int> &allIndicatorInLine,
                                                const std::vector<double> &digEPSS, const std::vector<double> &digPIEH,
                                                const std::vector<double> &digEPSL, const std::vector<double> &digEPSH, const std::vector<double> &digEPSB){
    std::vector<double> EqContents(22*readsNum,0.0);
    for(int i = 0; i < readsNum ; i++){
        if( allIndicatorInLine[at(i,0)] != 0 )  EqContents[at(i,0)]  =  (  2 * digEPSL[1]                                    + allLiksInLine[at(i,0)] ) ;
        if( allIndicatorInLine[at(i,1)] != 0 )  EqContents[at(i,1)]  =  (  2 * digEPSL[0]                                    + allLiksInLine[at(i,1)] ) ;
        if( allIndicatorInLine[at(i,2)] != 0 )  EqContents[at(i,2)]  =  (  ln2 + digEPSB[0] + digEPSL[0] + digEPSL[1]              + allLiksInLine[at(i,2)] ) ;
        if( allIndicatorInLine[at(i,3)] != 0 )  EqContents[at(i,3)]  =  (  ln2 + digEPSB[1] + digEPSL[0] + digEPSL[1]              + allLiksInLine[at(i,3)] ) ;

        if( allIndicatorInLine[at(i,4)] != 0 )  EqContents[at(i,4)]  =  (  digPIEH[0] + digEPSH[1]                           + allLiksInLine[at(i,4)] ) ;
        if( allIndicatorInLine[at(i,5)] != 0 )  EqContents[at(i,5)]  =  (  digPIEH[1] + digEPSH[1]                           + allLiksInLine[at(i,5)] ) ;
        if( allIndicatorInLine[at(i,6)] != 0 )  EqContents[at(i,6)]  =  (  digPIEH[0] + digEPSH[0]                           + allLiksInLine[at(i,6)] ) ;
        if( allIndicatorInLine[at(i,7)] != 0 )  EqContents[at(i,7)]  =  (  digPIEH[1] + digEPSH[0]                           + allLiksInLine[at(i,7)] ) ;

        if( allIndicatorInLine[at(i,8)] != 0 )  EqContents[at(i,8)]  =  (  digPIEH[0] + digEPSL[1] + digEPSL[1]              + allLiksInLine[at(i,8)] ) ;
        if( allIndicatorInLine[at(i,9)] != 0 )  EqContents[at(i,9)]  =  (  digPIEH[1] + digEPSL[1] + digEPSL[1]              + allLiksInLine[at(i,9)] ) ;
        if( allIndicatorInLine[at(i,10)] != 0 ) EqContents[at(i,10)] =  (  digPIEH[0] + digEPSL[0] + digEPSL[0]              + allLiksInLine[at(i,10)] ) ;
        if( allIndicatorInLine[at(i,11)] != 0 ) EqContents[at(i,11)] =  (  digPIEH[1] + digEPSL[0] + digEPSL[0]              + allLiksInLine[at(i,11)] ) ;
        if( allIndicatorInLine[at(i,12)] != 0 ) EqContents[at(i,12)] =  (  ln2 + digPIEH[0] + digEPSL[0] + digEPSL[1] + digEPSB[0] + allLiksInLine[at(i,12)] ) ;
        if( allIndicatorInLine[at(i,13)] != 0 ) EqContents[at(i,13)] =  (  ln2 + digPIEH[1] + digEPSL[0] + digEPSL[1] + digEPSB[0] + allLiksInLine[at(i,13)] ) ;
        if( allIndicatorInLine[at(i,14)] != 0 ) EqContents[at(i,14)] =  (  ln2 + digPIEH[0] + digEPSL[0] + digEPSL[1] + digEPSB[1] + allLiksInLine[at(i,14)] ) ;
        if( allIndicatorInLine[at(i,15)] != 0 ) EqContents[at(i,15)] =  (  ln2 + digPIEH[1] + digEPSL[0] + digEPSL[1] + digEPSB[1] + allLiksInLine[at(i,15)] ) ;

        if( allIndicatorInLine[at(i,16)] != 0 ) EqContents[at(i,16)] =  (  digEPSS[1]                                        + allLiksInLine[at(i,16)] ) ;
        if( allIndicatorInLine[at(i,17)] != 0 ) EqContents[at(i,17)] =  (  digEPSS[0] + digEPSB[1]                           + allLiksInLine[at(i,17)] ) ;
        if( allIndicatorInLine[at(i,18)] != 0 ) EqContents[at(i,18)] =  (  digEPSS[0] + digEPSB[0]                           + allLiksInLine[at(i,18)] ) ;

        if( allIndicatorInLine[at(i,19)] != 0 ) EqContents[at(i,19)] =  (  digEPSS[1]                                        + allLiksInLine[at(i,19)] ) ;
        if( allIndicatorInLine[at(i,20)] != 0 ) EqContents[at(i,20)] =  (  digEPSS[0] + digEPSB[0]                           + allLiksInLine[at(i,20)] ) ;
        if( allIndicatorInLine[at(i,21)] != 0 ) EqContents[at(i,21)] =  (  digEPSS[0] + digEPSB[1]                           + allLiksInLine[at(i,21)] ) ;
    }
    return EqContents;
}






OHVarfinDerModel::Result OHVarfinDerModel::calcTrueLowerBounds(const std::vector<std::vector<double> > &allLiksT,
                                               const std::vector<std::vector<double> > &allLiksN,
                                               const std::vector<std::vector<int> > &allIndicatorT,
                                               const std::vector<std::vector<int> > &allIndicatorN,
                                               const Parameters::BayesEMParameters &prior,
                                               bool isExome){
    std::vector<double > gammaF;    std::vector<double > prior_gammaF;
    std::vector<double > gammaH;    std::vector<double > prior_gammaH;
    std::vector<double > alphaL;    std::vector<double > prior_alphaL;
    std::vector<double > alphaH;    std::vector<double > prior_alphaH;
    std::vector<double > alphaB;    std::vector<double > prior_alphaB;
    std::vector<double > gammaEH;   std::vector<double > prior_gammaEH;
    std::vector<double > alphaS;    std::vector<double > prior_alphaS;
    {
        gammaF = prior_gammaF = prior.mut_a0;
        gammaH = prior_gammaH = prior.mut_b0;
        alphaL = prior_alphaL = prior.mut_c0;
        alphaH = prior_alphaH = prior.mut_d0;
        gammaEH = prior_gammaEH = prior.mut_e0;
        alphaS  = prior_alphaS  = prior.mut_f0;
        if(     isExome){   alphaB = prior_alphaB = prior.mut_g0;}
        else if(!isExome){  alphaB = prior_alphaB = prior.mut_h0;}
    }

    LOG(logINFO) << "=== OHVarfinDerModel::calcTumorLowerBounds Both Parameters ===" << std::endl;
    for(int i = 0 ; i < prior_gammaF.size();i++){LOG(logINFO) << "prior_gammaF [" << i  << "] :" << prior_gammaF[i] << std::endl;}
    for(int i = 0 ; i < prior_gammaH.size();i++){LOG(logINFO) << "prior_gammaH [" << i  << "] :" << prior_gammaH[i] << std::endl;}
    for(int i = 0 ; i < prior_alphaL.size();i++){LOG(logINFO) << "prior_alphaL [" << i  << "] :" << prior_alphaL[i] << std::endl;}
    for(int i = 0 ; i < prior_alphaH.size();i++){LOG(logINFO) << "prior_alphaH [" << i  << "] :" << prior_alphaH[i] << std::endl;}
    for(int i = 0 ; i < prior_gammaEH.size();i++){LOG(logINFO) << "prior_gammaEH [" << i  << "] :" << prior_gammaEH[i] << std::endl;}
    for(int i = 0 ; i < prior_alphaS.size();i++){LOG(logINFO) << "prior_alphaS [" << i  << "] :" << prior_alphaS[i] << std::endl;}
    LOG(logINFO) << "==========================================================" << std::endl;

    int readsNumT = allLiksT.size();
    int readsNumN = allLiksN.size();

    std::vector<double > EZT(allLiksT.size()*22, 0.0);
    std::vector<int> allIndicatorInLineT(readsNumT * 22);
    std::vector<double> allLiksInLineT(readsNumT * 22);
    {
        make2DArrayInLine<int>(allIndicatorInLineT, allIndicatorT, readsNumT, 22);
        make2DArrayInLine<double>(allLiksInLineT, allLiksT, readsNumT, 22);
        initEZ(EZT, allIndicatorInLineT, readsNumT);
    }

    std::vector<double > EZN(allLiksN.size()*22, 0.0);
    std::vector<int> allIndicatorInLineN(readsNumN * 22);
    std::vector<double> allLiksInLineN(readsNumN * 22);
    {
        make2DArrayInLine<int>(allIndicatorInLineN, allIndicatorN, readsNumN, 22);
        make2DArrayInLine<double>(allLiksInLineN, allLiksN, readsNumN, 22);
        initEZ(EZN, allIndicatorInLineN, readsNumN);
    }

    double beforeLq = -1.0 * 1e50;
    double nextLq = -1.0 * 1e50;
    bool first = true;
    int update = 0;

    //
    // Tumor param rules
    //

    std::vector<int> gamF0; std::vector<int> gamF1;
    {
        gamF0.push_back(0); gamF0.push_back(2); gamF0.push_back(3);  gamF0.push_back(16); gamF0.push_back(19);
        gamF1.push_back(1); gamF1.push_back(17);gamF1.push_back(18); gamF1.push_back(20); gamF1.push_back(21);
    }
    std::vector<int> gamH0; std::vector<int> gamH1; std::vector<int> gamH2;
    {
        gamH0.push_back(4); gamH0.push_back(8); gamH0.push_back(12); gamH0.push_back(14);
        gamH1.push_back(5); gamH1.push_back(7); gamH1.push_back(9); gamH1.push_back(11);
        gamH1.push_back(13); gamH1.push_back(15);
        gamH2.push_back(6); gamH2.push_back(10);
    }
    std::vector<int> alL0;  std::vector<int> alL1;
    {
        alL0.push_back(2);  alL0.push_back(3);
        alL0.push_back(11); alL0.push_back(11); alL0.push_back(12); alL0.push_back(13); alL0.push_back(14); alL0.push_back(15);
        alL1.push_back(0);  alL1.push_back(8);  alL1.push_back(9);  alL1.push_back(9);  alL1.push_back(13); alL1.push_back(15);
    }
    std::vector<int> alH0;  std::vector<int> alH1;
    {
        alH0.push_back(7); alH1.push_back(5);
    }
    std::vector<int> alB0; std::vector<int> alB1;
    {
        alB0.push_back(2); alB0.push_back(12); alB0.push_back(13); alB0.push_back(18); alB0.push_back(20);
        alB1.push_back(3); alB1.push_back(14); alB1.push_back(15); alB1.push_back(17); alB1.push_back(21);
    }
    //
    // Normal param rules
    //
    std::vector<int >  alLN0;  std::vector<int >  alLN1;
    {
        alLN0.push_back(1);  alLN0.push_back(1);  alLN0.push_back(2);  alLN0.push_back(3);
        alLN0.push_back(10); alLN0.push_back(10); alLN0.push_back(11); alLN0.push_back(11); alLN0.push_back(12); alLN0.push_back(13); alLN0.push_back(14); alLN0.push_back(15);
        alLN1.push_back(0);  alLN1.push_back(0);  alLN1.push_back(2);  alLN1.push_back(3);
        alLN1.push_back(8);  alLN1.push_back(8);  alLN1.push_back(9);  alLN1.push_back(9);  alLN1.push_back(12); alLN1.push_back(13); alLN1.push_back(14); alLN1.push_back(15);
    }
    std::vector<int >  alHN0; std::vector<int >  alHN1;
    {
        alHN0.push_back(6); alHN0.push_back(7);
        alHN1.push_back(4); alHN1.push_back(5);
    }
    std::vector<int >  gamEHN0; std::vector<int >  gamEHN1;
    {
        gamEHN0.push_back(4); gamEHN0.push_back(6); gamEHN0.push_back(8); gamEHN0.push_back(10); gamEHN0.push_back(12); gamEHN0.push_back(14);
        gamEHN1.push_back(5); gamEHN1.push_back(7); gamEHN1.push_back(9); gamEHN1.push_back(11); gamEHN1.push_back(13); gamEHN1.push_back(15);
    }
    std::vector<int >  alSN0; std::vector<int >  alSN1;
    {
        alSN0.push_back(17); alSN0.push_back(18); alSN0.push_back(20); alSN0.push_back(21);
        alSN1.push_back(16); alSN1.push_back(19);
    }
    std::vector<int >  alBN0; std::vector<int >  alBN1;
    {
        alBN0.push_back(2); alBN0.push_back(12); alBN0.push_back(13); alBN0.push_back(18); alBN0.push_back(20);
        alBN1.push_back(3); alBN1.push_back(14); alBN1.push_back(15); alBN1.push_back(17); alBN1.push_back(21);
    }

    std::vector<std::vector<int> > gamF; { gamF.push_back(gamF0);       gamF.push_back(gamF1);}
    std::vector<std::vector<int> > gamH; { gamH.push_back(gamH0);       gamH.push_back(gamH1); gamH.push_back(gamH2);}
    std::vector<std::vector<int> > alL;  { alL.push_back(alL0);         alL.push_back(alL1);}
    std::vector<std::vector<int> > alH;  { alH.push_back(alH0);         alH.push_back(alH1);}
    std::vector<std::vector<int> > alB;  { alB.push_back(alB0);         alB.push_back(alB1);}

    std::vector<std::vector<int> > alLN;  { alLN.push_back(alLN0);       alLN.push_back(alLN1);}
    std::vector<std::vector<int> > alHN;  { alHN.push_back(alHN0);       alHN.push_back(alHN1);}
    std::vector<std::vector<int> > gamEHN;{ gamEHN.push_back(gamEHN0);   gamEHN.push_back(gamEHN1);}
    std::vector<std::vector<int> > alSN;  { alSN.push_back(alSN0);       alSN.push_back(alSN1);}
    std::vector<std::vector<int> > alBN;  { alBN.push_back(alBN0);       alBN.push_back(alB1);}

    while( first || std::abs(beforeLq - nextLq) >=  minLqConvergence  ){
        // LOG(logINFO) << " ------ update count : " << update << " @ calcTumorLowerBounds in OHVarfinDerModel2" << std::endl;
        if(first) first = false;
        if( beforeLq > nextLq + minLqConvergence ){
            LOG(logERROR) << "decrease lower bound @ calcTrueLowerBounds Both in OHVarfinDerModel2" << std::endl;
        }

        beforeLq = nextLq;
        nextLq = 0.0;

        //M-step computation
        gammaF = prior_gammaF;
        gammaH = prior_gammaH;
        alphaL = prior_alphaL;
        alphaH = prior_alphaH;
        alphaB = prior_alphaB;
        gammaEH = prior_gammaEH;
        alphaS  = prior_alphaS;
        {
            updateParams(gammaF, gamF, EZT, allIndicatorInLineT, readsNumT);
            updateParams(gammaH, gamH, EZT, allIndicatorInLineT, readsNumT);
            updateParams(alphaL, alL,  EZT, allIndicatorInLineT, readsNumT);
            updateParams(alphaH, alH,  EZT, allIndicatorInLineT, readsNumT);
            updateParams(alphaB, alB,  EZT, allIndicatorInLineT, readsNumT);

            updateParams(gammaEH, gamEHN, EZN, allIndicatorInLineN, readsNumN);
            updateParams(alphaS,  alSN,   EZN, allIndicatorInLineN, readsNumN);
            updateParams(alphaL,  alLN,   EZN, allIndicatorInLineN, readsNumN);
            updateParams(alphaH,  alHN,   EZN, allIndicatorInLineN, readsNumN);
            updateParams(alphaB,  alBN,   EZN, allIndicatorInLineN, readsNumN);
        }

        std::vector<double > digPIF  = ModelUtils::calDirExp(gammaF);
        std::vector<double > digPIH  = ModelUtils::calDirExp(gammaH);
        std::vector<double > digEPSL = ModelUtils::calDirExp(alphaL);
        std::vector<double > digEPSH = ModelUtils::calDirExp(alphaH);
        std::vector<double > digEPSB  = ModelUtils::calDirExp(alphaB);
        std::vector<double > digPIEH = ModelUtils::calDirExp(gammaEH);
        std::vector<double > digEPSS = ModelUtils::calDirExp(alphaS);

        std::vector<double> EqContentsT = evaluateEqContentsTumor(readsNumT, allLiksInLineT, allIndicatorInLineT,
                                                                  digPIF, digPIH, digEPSL, digEPSH, digEPSB);
        std::vector<double> EqContentsN = evaluateEqContentsError(readsNumN, allLiksInLineN, allIndicatorInLineN,
                                                                  digEPSS, digPIEH, digEPSL, digEPSH, digEPSB);

        logSumExpAs2DArray(EZT, EqContentsT, allIndicatorInLineT, readsNumT, 22);
        logSumExpAs2DArray(EZN, EqContentsN, allIndicatorInLineN, readsNumN, 22);

        nextLq = 0.0;
        {
            nextLq += ModelUtils::calELogDir(digPIF, prior_gammaF);
            nextLq += ModelUtils::calELogDir(digPIH, prior_gammaH);
            nextLq += ModelUtils::calELogDir(digEPSL, prior_alphaL);
            nextLq += ModelUtils::calELogDir(digEPSH, prior_alphaH);
            nextLq += ModelUtils::calELogDir(digEPSB, prior_alphaB);
            nextLq += ModelUtils::calELogDir(digPIEH, prior_gammaEH);
            nextLq += ModelUtils::calELogDir(digEPSS, prior_alphaS);

            for(int i = 0; i < readsNumT; i++) for(int j = 0 ; j < 22; j++){
                if( allIndicatorInLineT[at(i,j)] != 0 ) nextLq += EZT[at(i,j)]  *   (  EqContentsT[at(i,j)] );
            }
            for(int i = 0; i < readsNumN; i++) for(int j = 0 ; j < 22; j++){
                if( allIndicatorInLineN[at(i,j)] != 0 ) nextLq += EZN[at(i,j)]  *   (  EqContentsN[at(i,j)] );
            }

            nextLq -= ModelUtils::calELogDir(digPIF, gammaF);
            nextLq -= ModelUtils::calELogDir(digPIH, gammaH);
            nextLq -= ModelUtils::calELogDir(digEPSL, alphaL);
            nextLq -= ModelUtils::calELogDir(digEPSH, alphaH);
            nextLq -= ModelUtils::calELogDir(digEPSB, alphaB);
            nextLq -= ModelUtils::calELogDir(digPIEH, gammaEH);
            nextLq -= ModelUtils::calELogDir(digEPSS, alphaS);

            for(int i = 0 ; i < readsNumT; i++) for (int j = 0 ; j < 22; j++){
                if( allIndicatorInLineT[at(i,j)] != 0 && EZT[at(i,j)] > 0.0 ) nextLq -= EZT[at(i,j)]*( log(EZT[at(i,j)]) );
            }

            for(int i = 0 ; i < readsNumN; i++) for (int j = 0 ; j < 22; j++){
                if( allIndicatorInLineN[at(i,j)] != 0 && EZN[at(i,j)] > 0.0 ) nextLq -= EZN[at(i,j)]*( log(EZN[at(i,j)]) );
            }
        }

        LOG(logINFO) << " end Lq comp" << std::endl;
        LOG(logINFO) << " before Lq "<< beforeLq << std::endl;
        LOG(logINFO) << " next Lq" << nextLq <<std::endl;
        update++;
        if(update > updateCount) break;
    }


    Result ans;
    Parameters::BayesEMParameters post;
{
    post.mut_a0 = gammaF;   post.mut_b0 = gammaH;
    post.mut_c0 = alphaL;   post.mut_d0 = alphaH;
    post.mut_e0 = alphaB;
}
    ans.posterior = post; ans.value = nextLq;

    std::vector<double> eLnHap = ModelUtils::calDirExp(gammaEH);
    if (!(eLnHap[0] > log(0.15) && eLnHap[0] < log(0.85))) {
        ans.st = UMBALANCED_NORMAL_HAPLOTYPE;
    }else{
        ans.st = SUCCESS;
    }

    LOG(logINFO) << "=== OHVarfinDerModel::calcTumorLowerBounds Both Parameters Posterior ===" << std::endl;
    for(int i = 0 ; i < gammaF.size();i++){LOG(logINFO) << "gammaF [" << i  << "] :" << gammaF[i] << std::endl;}
    for(int i = 0 ; i < gammaH.size();i++){LOG(logINFO) << "gammaH [" << i  << "] :" << gammaH[i] << std::endl;}
    for(int i = 0 ; i < alphaL.size();i++){LOG(logINFO) << "alphaL [" << i  << "] :" << alphaL[i] << std::endl;}
    for(int i = 0 ; i < alphaH.size();i++){LOG(logINFO) << "alphaH [" << i  << "] :" << alphaH[i] << std::endl;}
    for(int i = 0 ; i < gammaEH.size();i++){LOG(logINFO) << "gammaEH [" << i  << "] :" << gammaEH[i] << std::endl;}
    for(int i = 0 ; i < alphaS.size();i++){LOG(logINFO) << "alphaS [" << i  << "] :" << alphaS[i] << std::endl;}
    for(int i = 0 ; i < alphaB.size();i++){LOG(logINFO) << "alphaB [" << i  << "] :" << alphaB[i] << std::endl;}
    LOG(logINFO) << "==========================================================" << std::endl;
    return ans;
}




OHVarfinDerModel::Result OHVarfinDerModel::calcErrorLowerBounds(const std::vector<std::vector<double> > &allLiksT,
                                               const std::vector<std::vector<double> > &allLiksN,
                                               const std::vector<std::vector<int> > &allIndicatorT,
                                               const std::vector<std::vector<int> > &allIndicatorN,
                                               const Parameters::BayesEMParameters &prior,
                                               bool isExome){
    std::vector<double > alphaL;    std::vector<double > prior_alphaL;
    std::vector<double > alphaH;    std::vector<double > prior_alphaH;
    std::vector<double > gammaEH;   std::vector<double > prior_gammaEH;
    std::vector<double > alphaS;    std::vector<double > prior_alphaS;
    std::vector<double > alphaB;    std::vector<double > prior_alphaB;
    {
        if(isExome){
            alphaL  = prior_alphaL  = prior.err_a0;
            alphaH  = prior_alphaH  = prior.err_b0;
            gammaEH = prior_gammaEH = prior.err_c0;
            alphaS  = prior_alphaS  = prior.err_d0;
            alphaB  = prior_alphaB  = prior.err_e0;
        }else if(!isExome){
            alphaL  = prior_alphaL  = prior.err_f0;
            alphaH  = prior_alphaH  = prior.err_g0;
            gammaEH = prior_gammaEH = prior.err_h0;
            alphaS  = prior_alphaS  = prior.err_i0;
            alphaB  = prior_alphaB  = prior.err_j0;
        }
    }

    LOG(logINFO) << "=== OHVarfinDerModel::calcErrorLowerBounds Both Parameters ===" << std::endl;
    for(int i = 0 ; i < prior_alphaL.size();i++){LOG(logINFO) << "prior_alphaL [" << i  << "] :" << prior_alphaL[i] << std::endl;}
    for(int i = 0 ; i < prior_alphaH.size();i++){LOG(logINFO) << "prior_alphaH [" << i  << "] :" << prior_alphaH[i] << std::endl;}
    for(int i = 0 ; i < prior_gammaEH.size();i++){LOG(logINFO) << "prior_gammaEH [" << i  << "] :" << prior_gammaEH[i] << std::endl;}
    for(int i = 0 ; i < prior_alphaS.size();i++){LOG(logINFO) << "prior_alphaS [" << i  << "] :" << prior_alphaS[i] << std::endl;}
    LOG(logINFO) << "==========================================================" << std::endl;

    int readsNumT = allLiksT.size();
    int readsNumN = allLiksN.size();


    std::vector<double > EZT(allLiksT.size()*22, 0.0);
    std::vector<int> allIndicatorInLineT(readsNumT * 22);
    std::vector<double> allLiksInLineT(readsNumT   * 22);
    {
        make2DArrayInLine<int>(allIndicatorInLineT, allIndicatorT, readsNumT, 22);
        make2DArrayInLine<double>(allLiksInLineT, allLiksT, readsNumT, 22);
        initEZ(EZT, allIndicatorInLineT, readsNumT);
    }

    std::vector<double > EZN(allLiksN.size()*22, 0.0);
    std::vector<int> allIndicatorInLineN(readsNumN * 22);
    std::vector<double> allLiksInLineN(readsNumN   * 22);
    {
        make2DArrayInLine<int>(allIndicatorInLineN, allIndicatorN, readsNumN, 22);
        make2DArrayInLine<double>(allLiksInLineN, allLiksN, readsNumN, 22);
        initEZ(EZN, allIndicatorInLineN, readsNumN);
    }


    double beforeLq = -1.0 * 1e50;
    double nextLq   = -1.0 * 1e50;
    bool first = true;
    int update = 0;

    std::vector<int >  alL0;  std::vector<int >  alL1;
    {
        alL0.push_back(1);  alL0.push_back(1);  alL0.push_back(2);  alL0.push_back(3);
        alL0.push_back(10); alL0.push_back(10); alL0.push_back(11); alL0.push_back(11); alL0.push_back(12); alL0.push_back(13); alL0.push_back(14); alL0.push_back(15);
        alL1.push_back(0);  alL1.push_back(0);  alL1.push_back(2);  alL1.push_back(3);
        alL1.push_back(8);  alL1.push_back(8);  alL1.push_back(9);  alL1.push_back(9);  alL1.push_back(12); alL1.push_back(13); alL1.push_back(14); alL1.push_back(15);
    }
    std::vector<int >  alH0; std::vector<int >  alH1;
    {
        alH0.push_back(6); alH0.push_back(7);
        alH1.push_back(4); alH1.push_back(5);
    }
    std::vector<int >  gamEH0; std::vector<int >  gamEH1;
    {
        gamEH0.push_back(4); gamEH0.push_back(6); gamEH0.push_back(8); gamEH0.push_back(10); gamEH0.push_back(12); gamEH0.push_back(14);
        gamEH1.push_back(5); gamEH1.push_back(7); gamEH1.push_back(9); gamEH1.push_back(11); gamEH1.push_back(13); gamEH1.push_back(15);
    }
    std::vector<int >  alS0; std::vector<int >  alS1;
    {
        alS0.push_back(17); alS0.push_back(18); alS0.push_back(20); alS0.push_back(21);
        alS1.push_back(16); alS1.push_back(19);
    }
    std::vector<int >  alB0; std::vector<int >  alB1;
    {
        alB0.push_back(2); alB0.push_back(12); alB0.push_back(13); alB0.push_back(18); alB0.push_back(20);
        alB1.push_back(3); alB1.push_back(14); alB1.push_back(15); alB1.push_back(17); alB1.push_back(21);
    }
    std::vector<std::vector<int> > alL;{   alL.push_back(alL0);      alL.push_back(alL1);}
    std::vector<std::vector<int> > alH;{   alH.push_back(alH0);      alH.push_back(alH1);}
    std::vector<std::vector<int> > gamEH;{ gamEH.push_back(gamEH0);  gamEH.push_back(gamEH1);}
    std::vector<std::vector<int> > alS;{   alS.push_back(alS0);      alS.push_back(alS1);}
    std::vector<std::vector<int> > alB;{   alB.push_back(alB0);      alB.push_back(alB1);}

    while( first || std::abs(beforeLq - nextLq) >=  minLqConvergence  ){
        // LOG(logINFO) << " ------ update count : " << update << " @ calcErrorLowerBounds in OHVarfinDerModel2" << std::endl;
        if(first) first = false;
        if( beforeLq > nextLq + minLqConvergence ){
            LOG(logERROR) << "decrease lower bound @ calcErrorLowerBounds Both in OHVarfinDerModel2" << std::endl;
        }
        beforeLq = nextLq;
        nextLq = 0.0;

        //M-step computation
        alphaL = prior_alphaL;
        alphaH = prior_alphaH;
        gammaEH = prior_gammaEH;
        alphaS = prior_alphaS;
        alphaB = prior_alphaB;
        {
            updateParams(alphaL, alL,  EZT, allIndicatorInLineT, readsNumT);
            updateParams(alphaH, alH,  EZT, allIndicatorInLineT, readsNumT);
            updateParams(gammaEH, gamEH, EZT, allIndicatorInLineT, readsNumT);
            updateParams(alphaS, alS,  EZT, allIndicatorInLineT, readsNumT);
            updateParams(alphaB, alB,  EZT, allIndicatorInLineT, readsNumT);

            updateParams(alphaL, alL,  EZN, allIndicatorInLineN, readsNumN);
            updateParams(alphaH, alH,  EZN, allIndicatorInLineN, readsNumN);
            updateParams(gammaEH, gamEH, EZN, allIndicatorInLineN, readsNumN);
            updateParams(alphaS, alS,  EZN, allIndicatorInLineN, readsNumN);
            updateParams(alphaB, alB,  EZN, allIndicatorInLineN, readsNumN);
        }

        std::vector<double > digEPSL = ModelUtils::calDirExp(alphaL);
        std::vector<double > digEPSH = ModelUtils::calDirExp(alphaH);
        std::vector<double > digPIEH = ModelUtils::calDirExp(gammaEH);
        std::vector<double > digEPSS = ModelUtils::calDirExp(alphaS);
        std::vector<double > digEPSB  = ModelUtils::calDirExp(alphaB);

        //E-step computation
        std::vector<double> EqContentsT = evaluateEqContentsError(readsNumT, allLiksInLineT, allIndicatorInLineT,
                                                                  digEPSS, digPIEH, digEPSL, digEPSH, digEPSB);
        std::vector<double> EqContentsN = evaluateEqContentsError(readsNumN, allLiksInLineN, allIndicatorInLineN,
                                                                  digEPSS, digPIEH, digEPSL, digEPSH, digEPSB);

        logSumExpAs2DArray(EZT, EqContentsT, allIndicatorInLineT, readsNumT, 22);
        logSumExpAs2DArray(EZN, EqContentsN, allIndicatorInLineN, readsNumN, 22);

        // compute Lowerbound in each step
        nextLq = 0.0;
        {
            nextLq += ModelUtils::calELogDir(digEPSL, prior_alphaL);
            nextLq += ModelUtils::calELogDir(digEPSH, prior_alphaH);
            nextLq += ModelUtils::calELogDir(digPIEH, prior_gammaEH);
            nextLq += ModelUtils::calELogDir(digEPSS, prior_alphaS);
            nextLq += ModelUtils::calELogDir(digEPSB, prior_alphaB);

            for(int i = 0; i < readsNumT; i++)for(int j = 0 ; j < 22; j++){
                if( allIndicatorInLineT[at(i,j)] != 0 ) nextLq += EZT[at(i,j)]  *   (  EqContentsT[at(i,j)] );
            }

            for(int i = 0; i < readsNumN; i++)for(int j = 0 ; j < 22; j++){
                if( allIndicatorInLineN[at(i,j)] != 0 ) nextLq += EZN[at(i,j)]  *   (  EqContentsN[at(i,j)] );
            }

            nextLq -= ModelUtils::calELogDir(digEPSL, alphaL);
            nextLq -= ModelUtils::calELogDir(digEPSH, alphaH);
            nextLq -= ModelUtils::calELogDir(digPIEH, gammaEH);
            nextLq -= ModelUtils::calELogDir(digEPSS, alphaS);
            nextLq -= ModelUtils::calELogDir(digEPSB, alphaB);

            for(int i = 0 ; i < readsNumT; i++)for (int j = 0 ; j < 22; j++){
                if( allIndicatorInLineT[at(i,j)] != 0 && EZT[at(i,j)] > 0.0 ) nextLq -= EZT[at(i,j)]*( log(EZT[at(i,j)]) );
            }

            for(int i = 0 ; i < readsNumN; i++)for (int j = 0 ; j < 22; j++){
                if( allIndicatorInLineN[at(i,j)] != 0 && EZN[at(i,j)] > 0.0 ) nextLq -= EZN[at(i,j)]*( log(EZN[at(i,j)]) );
            }
        }
        LOG(logINFO) << " end Lq comp" << std::endl;
        LOG(logINFO) << " before Lq "<< beforeLq << std::endl;
        LOG(logINFO) << " next Lq" << nextLq <<std::endl;
        update++;
        if(update > updateCount) break;
    }

    LOG(logINFO) << "=== OHVarfinDerModel::calcErrorLowerBounds Both PosteriorParameters ===" << std::endl;
    for(int i = 0 ; i < alphaL.size();i++){LOG(logINFO) << "alphaL [" << i  << "] :" << alphaL[i] << std::endl;}
    for(int i = 0 ; i < alphaH.size();i++){LOG(logINFO) << "alphaH [" << i  << "] :" << alphaH[i] << std::endl;}
    for(int i = 0 ; i < gammaEH.size();i++){LOG(logINFO) << "gammaEH [" << i  << "] :" << gammaEH[i] << std::endl;}
    for(int i = 0 ; i < alphaS.size();i++){LOG(logINFO) << "alphaS [" << i  << "] :" << alphaS[i] << std::endl;}
    LOG(logINFO) << "==========================================================" << std::endl;

    Result ans;
    Parameters::BayesEMParameters post;
{
     // (!isTumor) && (!isTrueModel)  ){
    post.err_f0 = alphaL  ;        post.err_g0 = alphaH  ;
    post.err_h0 = gammaEH ;        post.err_i0 = alphaS  ;
    post.err_j0 = alphaB  ;

}
    ans.posterior = post; ans.value = nextLq;
    ans.st = SUCCESS;
    return ans;
}


OHVarfinDerModel::Result OHVarfinDerModel::calcErrorLowerBounds(const std::vector<std::vector<double> > &allLiks,
                                              const std::vector<std::vector<int> > &allIndicator,
                                              const Parameters::BayesEMParameters &prior, bool isExome){
    std::vector<double > alphaL;    std::vector<double > prior_alphaL;
    std::vector<double > alphaH;    std::vector<double > prior_alphaH;
    std::vector<double > gammaEH;   std::vector<double > prior_gammaEH;
    std::vector<double > alphaS;    std::vector<double > prior_alphaS;
    std::vector<double > alphaB;    std::vector<double > prior_alphaB;
    {
        if( isExome ){
            alphaL  = prior_alphaL  = prior.err_a0;
            alphaH  = prior_alphaH  = prior.err_b0;
            gammaEH = prior_gammaEH = prior.err_c0;
            alphaS  = prior_alphaS  = prior.err_d0;
            alphaB  = prior_alphaB  = prior.err_e0;
        }else if( !isExome ) {
            alphaL  = prior_alphaL  = prior.err_f0;
            alphaH  = prior_alphaH  = prior.err_g0;
            gammaEH = prior_gammaEH = prior.err_h0;
            alphaS  = prior_alphaS  = prior.err_i0;
            alphaB  = prior_alphaB  = prior.err_j0;
        }
    }

    LOG(logINFO) <<  "alphaB.size() : " << alphaB.size() << std::endl;
    LOG(logINFO) << "=== OHVarfinDerModel::calcErrorLowerBounds Parameters ===" << std::endl;
    if(isExome) LOG(logINFO) << "Exome Error Model" << std::endl;
    else if(!isExome) LOG(logINFO) << "Whole Error Model" << std::endl;

    for(int i = 0 ; i < prior_alphaL.size();i++){ LOG(logINFO) << "prior_alphaL [" << i  << "] :" << prior_alphaL[i] << std::endl;}
    for(int i = 0 ; i < prior_alphaH.size();i++){ LOG(logINFO) << "prior_alphaH [" << i  << "] :" << prior_alphaH[i] << std::endl;}
    for(int i = 0 ; i < prior_gammaEH.size();i++){LOG(logINFO) << "prior_gammaEH [" << i  << "] :" << prior_gammaEH[i] << std::endl;}
    for(int i = 0 ; i < prior_alphaS.size();i++){ LOG(logINFO) << "prior_alphaS [" << i  << "] :" << prior_alphaS[i] << std::endl;}
    for(int i = 0 ; i < prior_alphaB.size();i++){LOG(logINFO) << "prior_alphaB [" << i  << "] :" << prior_alphaB[i] << std::endl;}
    LOG(logINFO) << "==========================================================" << std::endl;

    int readsNum = allLiks.size();
    std::vector<double > EZ(allLiks.size()*22, 0.0);
    std::vector<int> allIndicatorInLine(readsNum * 22);
    std::vector<double> allLiksInLine(readsNum   * 22);
    {
        make2DArrayInLine<int>(allIndicatorInLine, allIndicator, readsNum, 22);
        make2DArrayInLine<double>(allLiksInLine, allLiks, readsNum, 22);
        initEZ(EZ, allIndicatorInLine, readsNum);
    }

    double beforeLq = -1.0 * 1e50;
    double nextLq   = -1.0 * 1e50;
    bool first = true;
    int update = 0;

    std::vector<int >  alL0;  std::vector<int >  alL1;
    {
        alL0.push_back(1);  alL0.push_back(1);  alL0.push_back(2);  alL0.push_back(3);
        alL0.push_back(10); alL0.push_back(10); alL0.push_back(11); alL0.push_back(11); alL0.push_back(12); alL0.push_back(13); alL0.push_back(14); alL0.push_back(15);
        alL1.push_back(0);  alL1.push_back(0);  alL1.push_back(2);  alL1.push_back(3);
        alL1.push_back(8);  alL1.push_back(8);  alL1.push_back(9);  alL1.push_back(9);  alL1.push_back(12); alL1.push_back(13); alL1.push_back(14); alL1.push_back(15);
    }
    std::vector<int >  alH0; std::vector<int >  alH1;
    {
        alH0.push_back(6); alH0.push_back(7);
        alH1.push_back(4); alH1.push_back(5);
    }
    std::vector<int >  gaEH0; std::vector<int >  gaEH1;
    {
        gaEH0.push_back(4); gaEH0.push_back(6); gaEH0.push_back(8); gaEH0.push_back(10); gaEH0.push_back(12); gaEH0.push_back(14);
        gaEH1.push_back(5); gaEH1.push_back(7); gaEH1.push_back(9); gaEH1.push_back(11); gaEH1.push_back(13); gaEH1.push_back(15);
    }
    std::vector<int >  alS0; std::vector<int >  alS1;
    {
        alS0.push_back(17); alS0.push_back(18); alS0.push_back(20); alS0.push_back(21);
        alS1.push_back(16); alS1.push_back(19);
    }
    std::vector<int >  alB0; std::vector<int >  alB1;
    {
        alB0.push_back(2); alB0.push_back(12); alB0.push_back(13); alB0.push_back(18); alB0.push_back(20);
        alB1.push_back(3); alB1.push_back(14); alB1.push_back(15); alB1.push_back(17); alB1.push_back(21);
    }
    std::vector<std::vector<int> > alL, alH, gaEH, alS, alB;
    {
        alL.push_back(alL0);   alL.push_back(alL1);
        alH.push_back(alH0);   alH.push_back(alH1);
        gaEH.push_back(gaEH0); gaEH.push_back(gaEH1);
        alS.push_back(alS0);   alS.push_back(alS1);
        alB.push_back(alB0);   alB.push_back(alB1);
    }

    while( first || std::abs(beforeLq - nextLq) >=  minLqConvergence  ){
        LOG(logINFO) << " ------ update count : " << update << " @ calcErrorLowerBounds in OHVarfinDerModel2" << std::endl;
        if(first) first = false;
        if( beforeLq > nextLq + minLqConvergence ){
            LOG(logERROR) << "decrease lower bound @ calcErrorLowerBounds in OHVarfinDerModel2" << std::endl;
        }

        beforeLq = nextLq;
        nextLq = 0.0;

        // M-step computation
        alphaL = prior_alphaL;
        alphaH = prior_alphaH;
        gammaEH = prior_gammaEH;
        alphaS = prior_alphaS;
        alphaB = prior_alphaB;
        {
            updateParams(alphaL,  alL,  EZ, allIndicatorInLine, readsNum);
            updateParams(alphaH,  alH,  EZ, allIndicatorInLine, readsNum);
            updateParams(gammaEH, gaEH, EZ, allIndicatorInLine, readsNum);
            updateParams(alphaS,  alS,  EZ, allIndicatorInLine, readsNum);
            updateParams(alphaB,  alB,  EZ, allIndicatorInLine, readsNum);
        }

        std::vector<double > digEPSL = ModelUtils::calDirExp(alphaL);
        std::vector<double > digEPSH = ModelUtils::calDirExp(alphaH);
        std::vector<double > digPIEH = ModelUtils::calDirExp(gammaEH);
        std::vector<double > digEPSS = ModelUtils::calDirExp(alphaS);
        std::vector<double > digEPSB = ModelUtils::calDirExp(alphaB);

        //E-step computation
        std::vector<double> EqContents = evaluateEqContentsError(readsNum, allLiksInLine, allIndicatorInLine,
                                                                  digEPSS, digPIEH, digEPSL, digEPSH, digEPSB);

        logSumExpAs2DArray(EZ, EqContents, allIndicatorInLine, readsNum, 22);
        // compute Lowerbound in each step
        nextLq = 0.0;
        {
            nextLq += ModelUtils::calELogDir(digEPSL, prior_alphaL);
            nextLq += ModelUtils::calELogDir(digEPSH, prior_alphaH);
            nextLq += ModelUtils::calELogDir(digPIEH, prior_gammaEH);
            nextLq += ModelUtils::calELogDir(digEPSS, prior_alphaS);
            nextLq += ModelUtils::calELogDir(digEPSB, prior_alphaB);

            for(int i = 0; i < readsNum; i++)for(int j = 0 ; j < 22; j++){
                if( allIndicatorInLine[at(i,j)] != 0 ) nextLq += EZ[at(i,j)]  *   (  EqContents[at(i,j)] );
            }
            nextLq -= ModelUtils::calELogDir(digEPSL, alphaL);
            nextLq -= ModelUtils::calELogDir(digEPSH, alphaH);
            nextLq -= ModelUtils::calELogDir(digPIEH, gammaEH);
            nextLq -= ModelUtils::calELogDir(digEPSS, alphaS);
            nextLq -= ModelUtils::calELogDir(digEPSB, alphaB);

            for(int i = 0 ; i < readsNum; i++)for (int j = 0 ; j < 22; j++){
                if( allIndicatorInLine[at(i,j)] != 0 && EZ[at(i,j)] > 0.0 ) nextLq -= EZ[at(i,j)]*( log(EZ[at(i,j)]) );
            }
        }
        LOG(logINFO) << " end Lq comp" << std::endl;
        LOG(logINFO) << " before Lq "  << beforeLq << std::endl;
        LOG(logINFO) << " next Lq"     << nextLq   << std::endl;
        update++;
        if(update > updateCount) break;
    }

    Result ans;
    Parameters::BayesEMParameters post;
    {
        if( isExome ){
            post.err_a0 = alphaL  ;        post.err_b0 = alphaH  ;
            post.err_c0 = gammaEH ;        post.err_d0 = alphaS  ;
            post.err_e0 = alphaB  ;
        }else if( !isExome ) {
            post.err_f0 = alphaL  ;        post.err_g0 = alphaH  ;
            post.err_h0 = gammaEH ;        post.err_i0 = alphaS  ;
            post.err_j0 = alphaB  ;
        }
    }
    ans.posterior = post; ans.value = nextLq;
    ans.st = SUCCESS;
    return ans;
}

OHVarfinDerModel::Result OHVarfinDerModel::estimate(){
	Result ans;
    Parameters::BayesEMParameters priorError = prior;
    Parameters::BayesEMParameters priorTumor = prior;
    if(thresModelChangeDepth < (tumorLiks.size() + normalLiks.size()) * 0.5 ){
        LOG(logINFO) << " exome model " << std::endl;
        Result rd1   =  calcTrueLowerBounds(tumorLiks, normalLiks, tumorIndicator, normalIndicator, priorTumor, true);
        Result ri3_e = calcErrorLowerBounds(normalLiks, normalIndicator, priorError, true);
        Result ri4_e = calcErrorLowerBounds(tumorLiks,  tumorIndicator,  priorError, true);
        ans.value = rd1.value - (ri3_e.value + ri4_e.value);
        ans.posterior = rd1.posterior;
        ans.st = getState3(rd1, ri3_e, ri4_e);
    }else{
        LOG(logINFO) << " whole model " << std::endl;
        Result rd1   =  calcTrueLowerBounds(tumorLiks, normalLiks, tumorIndicator, normalIndicator, priorTumor, false);
        Result rd2_w = calcErrorLowerBounds(tumorLiks, normalLiks, tumorIndicator, normalIndicator, priorError, false);
        ans.value = rd1.value - rd2_w.value;
        ans.posterior = rd1.posterior;
        ans.st = getState2(rd1, rd2_w);
    }
    return ans;
}

OHVarfinDerModel::State OHVarfinDerModel::getState2(Result r1, Result r2){
    if(     r1.st != SUCCESS){ return r1.st; }
    else if(r2.st != SUCCESS){ return r2.st; }

    return SUCCESS;
}
OHVarfinDerModel::State OHVarfinDerModel::getState3(Result r1, Result r2, Result r3){
    State tmp;
    if((tmp = getState2(r1,r2)) != SUCCESS){ return tmp;  }
    else if(r3.st != SUCCESS){               return r3.st;}

    return SUCCESS;
}

void OHVarfinDerModel::initEZ(std::vector<double> &EZ, const std::vector<int> &allIndicatorInLine, int readsNum){
    for(int i = 0; i < readsNum; i++){
        if(allIndicatorInLine[at(i,0)] != 0){
            EZ[at(i,0)]  = 1.0;
        }else if(allIndicatorInLine[at(i,4)] != 0){
            EZ[at(i,4)]  = 1.0;
        }else if(allIndicatorInLine[at(i,8)] != 0){
            EZ[at(i,8)]  = 1.0;
        }else if(allIndicatorInLine[at(i,16)]!= 0){
            EZ[at(i,16)] = 0.5; EZ[at(i,17)] = 0.5;
        }else if(allIndicatorInLine[at(i,19)]!= 0){
            EZ[at(i,19)] = 0.5; EZ[at(i,20)] = 0.5;
        }
    }
}

template <typename T>
void OHVarfinDerModel::make2DArrayInLine(std::vector<T> &ans, const std::vector< std::vector<T> > &from, int colSize, int rowSize){
    if(ans.size() != colSize*rowSize){ ans = std::vector<T>(colSize*rowSize);}

    for(int i = 0; i < colSize; i++)for(int j = 0 ; j < rowSize ; j++){ans[at(i,j)] = from[i][j];}
}


std::vector<double> OHVarfinDerModel::logSumExpAs2DArray(const std::vector<double> &Eq, const std::vector<int> &mask, int colSize, int rowSize){
    std::vector< double > ans(colSize * rowSize, 0.0);
    if(Eq.size() != colSize*rowSize) {
        LOG(logINFO) << "Eq.size() != colSize*rowSize @ OHVarfinDerModel::logSumExpAs2DArray" << std::endl;
        throw std::string("Eq.size() != colSize*rowSize @ OHVarfinDerModel::logSumExpAs2DArray");
        return ans;
    }
    for(int i = 0; i < colSize; i++){
        double sum   = 0.0;
        double maxEq = -1e100;
        for(int j = 0; j < rowSize; j++) ans[at(i,j)] = 0.0;

        for(int j = 0; j < rowSize; j++)if( mask[at(i,j)] != 0 ){
            maxEq  = std::max( maxEq, Eq[at(i,j)] );
        }
        for(int j = 0; j < rowSize; j++)if( mask[at(i,j)] != 0 ){
            ans[at(i,j)]  += exp( Eq[at(i,j)] - maxEq );
        }
        for(int j = 0; j < rowSize; j++) sum += ans[at(i,j)];
        for(int j = 0; j < rowSize; j++) ans[at(i,j)] /= sum;
    }
    return ans;
}

void OHVarfinDerModel::logSumExpAs2DArray(std::vector<double>& ans, const std::vector<double> &Eq, const std::vector<int> &mask, int colSize, int rowSize){
    if(ans.size() != colSize * rowSize){
        ans = std::vector< double >(colSize * rowSize, 0.0);
    }
    if(Eq.size() != colSize*rowSize) {
        LOG(logINFO) << "Eq.size() != colSize*rowSize @ OHVarfinDerModel::logSumExpAs2DArray" << std::endl;
        throw std::string("Eq.size() != colSize*rowSize @ OHVarfinDerModel::logSumExpAs2DArray");
    }
    for(int i = 0; i < colSize; i++){
        double sum   = 0.0;
        double maxEq = -1e100;
        for(int j = 0; j < rowSize; j++) ans[at(i,j)] = 0.0;

        for(int j = 0; j < rowSize; j++)if( mask[at(i,j)] != 0 ){
            maxEq  = std::max( maxEq, Eq[at(i,j)] );
        }
        for(int j = 0; j < rowSize; j++)if( mask[at(i,j)] != 0 ){
            ans[at(i,j)]  += exp( Eq[at(i,j)] - maxEq );
        }
        for(int j = 0; j < rowSize; j++) sum += ans[at(i,j)];
        for(int j = 0; j < rowSize; j++) ans[at(i,j)] /= sum;
    }
}
