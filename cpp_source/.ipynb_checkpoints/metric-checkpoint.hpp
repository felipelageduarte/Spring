#ifndef metric_h
#define metric_h

#include<vector>
using namespace std;

class metric{
    public:
        virtual double evaluate(vector<double> x, vector<double> y);
};

#endif