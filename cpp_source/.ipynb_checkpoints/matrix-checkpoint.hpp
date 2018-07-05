#ifndef MATRIX_H
#define MATRIX_H

#include <iostream>
#include <iomanip>
#include <vector>
#include <cassert>
using namespace std; 

class matrix {
    protected:
        vector< vector< double > > data; 
        int m, n; 
    public:
        explicit matrix(int m_, int n_);
        vector<double>& operator[](int i);
        double get(int i, int j);
        int rows();
        int cols();
        matrix operator+(matrix& B);
        matrix operator*(double c);
        matrix operator*(matrix& B);
        double trace();
        matrix tp();
        void write(ostream& out);
        void read(istream& in);
}; 

#endif