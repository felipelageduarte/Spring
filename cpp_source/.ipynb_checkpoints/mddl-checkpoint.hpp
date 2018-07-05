#ifndef RCPP_MDDL_H
#define RCPP_MDDL_H

#include "metric.hpp"
#include "matrix.hpp"

class mddl: public metric{
    
    public: 
        double evaluate(vector<double> x, vector<double> y);
    
    private:
        int whichMin( double x, double y, double z);
        matrix cost_matrix(vector<double> x, vector<double> y);
        matrix align(matrix costMat);

};

#endif