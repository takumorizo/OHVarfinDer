//
//  OHVarfinder2Model.hpp
//  OHVarFinder
//
//  Created by 森山卓也 on 2016/07/06.
//  Copyright © 2016年 森山卓也. All rights reserved.
//

#ifndef OHVarfinder2Model_hpp
#define OHVarfinder2Model_hpp


#include <stdio.h>
#include "math.h"
#include <stdio.h>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <fstream>
#include <utility>
#include "Variant.h"
#include "Parameters.h"
#include "Haplotype.hpp"

// TODO : remove comments
class OHVarfinDer2Model  {
public:
    const Parameters::BayesEMParameters &prior;
    const std::vector<std::vector<double> > &tumorLiks, &normalLiks;
    const std::vector<std::vector<int> > &tumorIndicator, &normalIndicator;
    bool isErrorModel;
    const int  updateCount;
    const double minLqConvergence;
    const int thresModelChangeDepth;
    typedef enum { SUCCESS, UMBALANCED_NORMAL_HAPLOTYPE } State;
    struct Result {
        Parameters::BayesEMParameters posterior;
        State  st;
        double value;
    };
    OHVarfinDer2Model(const std::vector<std::vector<double> > &tumorLiks, const std::vector<std::vector<double> > &normalLiks,
                     const std::vector<std::vector<int> > &tumorIndicator, const std::vector<std::vector<int> > &normalIndicator,
                     const Parameters::BayesEMParameters &prior,
                    int updateCount = 500, double minLqConvergence = 0.0001, int thresModelChangeDepth = 100);  //est_type=mutation or non-mutation
    // double estimate(bool isTrueModel);
    Result estimate();
private:
    // double calcTumorLowerBounds(const std::vector<std::vector<double> > &allLiks,
    //                                     const std::vector<std::vector<int> > &allIndicator );
    // double calcErrorLowerBounds(const std::vector<std::vector<double> > &allLiks,
    //                                     const std::vector<std::vector<int> > &allIndicator, bool isTumor, bool isTureModel );
    Result calcTrueLowerBounds(const std::vector<std::vector<double> > &allLiksT,
                               const std::vector<std::vector<double> > &allLiksN,
                               const std::vector<std::vector<int> > &allIndicatorT,
                               const std::vector<std::vector<int> > &allIndicatorN,
                               const Parameters::BayesEMParameters &prior,
                               bool isExome);

    Result calcErrorLowerBounds(const std::vector<std::vector<double> > &allLiks,
                                const std::vector<std::vector<int> > &allIndicator,
                                const Parameters::BayesEMParameters &prior,
                                bool isExome);

    Result calcErrorLowerBounds(const std::vector<std::vector<double> > &allLiksT,
                                const std::vector<std::vector<double> > &allLiksN,
                                const std::vector<std::vector<int> > &allIndicatorT,
                                const std::vector<std::vector<int> > &allIndicatorN,
                                const Parameters::BayesEMParameters &prior,
                                bool isExome);
    void updateParams(std::vector<double>& param,    const std::vector<std::vector<int> > &rule,
                      const std::vector<double> &EZ, const std::vector<int> &allIndicatorInLine,
                      int readNum);

    std::vector<double> evaluateEqContentsTumor(
                                int readsNum, const std::vector<double> &allLiksInLine, const std::vector<int> &allIndicatorInLine,
                                const std::vector<double> &digPIF,  const std::vector<double> &digPIH, 
                                const std::vector<double> &digEPSL, const std::vector<double> &digEPSH, const std::vector<double> &digEPSB );

    std::vector<double> evaluateEqContentsError(
                                int readsNum, const std::vector<double> &allLiksInLine, const std::vector<int> &allIndicatorInLine,
                                const std::vector<double> &digEPSS, const std::vector<double> &digPIEH,
                                const std::vector<double> &digEPSL, const std::vector<double> &digEPSH, const std::vector<double> &digEPSB);


    int at(int x, int y){
        return 22*x+y;
    }

    OHVarfinDer2Model::State getState2(Result r1, Result r2);
    OHVarfinDer2Model::State getState3(Result r1, Result r2, Result r3);

    void initEZ(std::vector<double> &EZ, const std::vector<int> &allIndicatorInLine, int readsNum);
    template <typename T> void make2DArrayInLine(std::vector<T> &ans, const std::vector< std::vector<T> > &from, int colSize, int rowSize = 22);
    std::vector<double> logSumExpAs2DArray(const std::vector<double> &Eq, const std::vector<int> &mask, int colSize,  int rowSize =22);
    void logSumExpAs2DArray(std::vector<double> &ans, const std::vector<double> &Eq, const std::vector<int> &mask, int colSize,  int rowSize =22);
};


#endif /* OHVarfinder2Model_hpp */
