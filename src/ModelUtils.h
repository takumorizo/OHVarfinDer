//
//  ModelUtils.h
//  OHVarFinder
//
//  Created by 森山卓也 on 2016/04/05.
//  Copyright (c) 2016年 森山卓也. All rights reserved.
//

#ifndef ModelUtils_h
#define ModelUtils_h

#include <math.h>
#include <boost/math/special_functions/digamma.hpp>
#include <boost/foreach.hpp>
#include "MutationModel.hpp"
#include "log.h"
#include "Haplotype.hpp"

class ModelUtils{
public:
    static inline double sumVec(std::vector<double> ar) {
        double sum = 0.0;
        for (int i = 0; i < ar.size(); i++) {
            sum += ar[i];
        }
        return sum;
    }
    
    static inline double sumVec(std::vector<double> ar,std::vector<int> indicator) {
        double sum = 0.0;
        for (int i = 0; i < ar.size(); i++)if(indicator[i]!=0) {
            sum += ar[i];
        }
        return sum;
    }
    
    //l_e_pit <- digamma(new_al) - digammma(sum(new_al))
    static inline std::vector<double> calDirExp(std::vector<double> ar) {
        double sum = sumVec(ar);
        std::vector<double> ans(ar.size(), 0.0);
        try {
            LOG(logINFO) << "start cal dig_sum for sum : " << sum << std::endl;
            double dig_sum = boost::math::digamma(sum);
            LOG(logINFO) << "start cal dig iteratively " << std::endl;
            for (int i = 0; i < ar.size(); i++) {
                ans[i] = boost::math::digamma(ar[i]) - dig_sum;
            }
        } catch (std::exception &e) {
            std::string message = std::string("error_exception_").append(e.what());
            LOG(logDEBUG) << message << std::endl;
            LOG(logDEBUG) << "calDirExp" << std::endl;
            BOOST_FOREACH(double x, ar) {
                LOG(logDEBUG) << x << " ";
            }
            LOG(logDEBUG) << std::endl;
            LOG(logDEBUG) << "sum: " << sum << std::endl;
            LOG(logDEBUG).flush();
            throw;
        }
        return ans;
    }
    
    
    static inline std::vector<double> normVec(std::vector<double> ar){
        double sum = sumVec(ar);
        std::vector<double> ans(ar.size(), 0.0);
        for (int i = 0; i < ar.size(); i++) {
            ans[i] = ar[i] / sum;
        }
        return ans;
    }
    
    static inline std::vector<double> normVec(std::vector<double> ar,std::vector<int> indicator){
        double sum = sumVec(ar,indicator);
        std::vector<double> ans(ar.size(), 0.0);
        for (int i = 0; i < ar.size(); i++) {
            ans[i] = ar[i] / sum;
        }
        return ans;
    }
    
    static inline std::vector<double> expVec(std::vector<double> ar) {
        std::vector<double> ans(ar.size(), 0.0);
        for (int i = 0; i < ar.size(); i++) {
            ans[i] = exp(ar[i]);
        }
        return ans;
    }
    
    static inline std::vector<double> expVec(std::vector<double> ar,std::vector<int> indicator) {
        std::vector<double> ans(ar.size(), 0.0);
        for (int i = 0; i < ar.size(); i++)if( indicator[i] != 0 ) {
            ans[i] = exp(ar[i]);
        }
        return ans;
    }
    // log_theta[0] =  Eq[ log_Var[0] ]
    static inline double calELogDir(std::vector<double> log_theta, std::vector<double> param) {
        LOG(logDEBUG) << " start calELogDir " << std::endl;
        //E[log p(theta|param)] : Dir
        //- sum(lgamma(param)) + lgamma(sum(param)) + sum((param-1) * log_theta)
        double sum_lgam_param = 0.0;
        for (int i = 0; i < param.size(); i++) {
            sum_lgam_param -= lgamma(param[i]);
            LOG(logDEBUG) << " done lgamma " << std::endl;
        }
        double lgam_sum_param = lgamma(sumVec(param));
        double sum_par_lt = 0.0;
        for (int i = 0; i < param.size(); i++) {
            sum_par_lt += (param[i] - 1) * log_theta[i];
        }
        if (log_theta.size() != param.size()) {
            LOG(logDEBUG) << "--- calELogDir ---" << std::endl;
            printVec(log_theta);
            LOG(logDEBUG) << std::endl;
            printVec(param);
            LOG(logDEBUG) << std::endl;
            LOG(logDEBUG) << sum_lgam_param << " " << lgam_sum_param << " " << sum_par_lt << std::endl;
            throw std::string("calELogDir");
        }
        return sum_lgam_param + lgam_sum_param + sum_par_lt;
    }
    template <class X> static void printVec(std::vector<X> vec) {
        BOOST_FOREACH(X x, vec) {
            LOGP(logDEBUG) << x << " ";
        }
    }

};


#endif
