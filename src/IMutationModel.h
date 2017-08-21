//
//  IMutationModel.h
//  OHVarFinder
//
//  Created by 森山卓也 on 2016/04/04.
//  Copyright (c) 2016年 森山卓也. All rights reserved.
//

#ifndef __OHVarFinder__IMutationModel__
#define __OHVarFinder__IMutationModel__

#include <stdio.h>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <fstream>
#include "Variant.h"
#include "Parameters.h"
#include "Haplotype.hpp"

class IMutationModel {
public:
    virtual double estimate(bool useTumor = true, bool setPosterior = true) = 0;
private:
    virtual double calcTumorLowerBounds(const std::vector<std::vector<double> > &allLiks,
                                        const std::vector<std::vector<int> > &allIndicator, bool setPosterior = true) = 0;
    virtual double calcErrorLowerBounds(const std::vector<std::vector<double> > &allLiks,
                                        const std::vector<std::vector<int> > &allIndicator, bool setPosterior = true) = 0;
    virtual int at(int x, int y) = 0;
};




#endif /* defined(__OHVarFinder__IMutationModel__) */
