#include "matrix.hpp"

#include <iostream>
#include <iomanip>
#include <vector>
#include <cassert>
using namespace std; 

matrix::matrix(int m_, int n_) {
    m = m_; n = n_; 
    data.resize(m+1); 
    for(int i=0; i<=m; i++)
        data[i].resize(n+1); 
}

vector<double>& matrix::operator[](int i) { return data[i]; }
double matrix::get(int i, int j) { return data[i][j]; }
int matrix::rows() { return m; }
int matrix::cols() { return n; }

matrix matrix::operator+(matrix& B) { // Return self+B
    assert(m==B.m && n==B.n); 
    matrix C(m, n); 
    for(int i=1; i<=m; i++)
        for(int j=1; j<=n; j++)
            C.data[i][j] = data[i][j] + B.data[i][j]; 
    return C; 
}

matrix matrix::operator*(double c) { // Scalar multiplication
    matrix C(m, n);
    for(int i=1; i<=m; i++)
        for(int j=1; j<=n; j++)
            C.data[i][j] = c * data[i][j]; 
    return C; 
}

matrix matrix::operator*(matrix& B) { // Return self * B
    int p = B.n; 
    assert(n == B.m); 
    matrix C(m, p); 
    for(int i=1; i<=m; i++)
        for(int j=1; j<=p; j++)
            for(int k=1; k<=n; k++)
                C.data[i][j] += data[i][k] * B.data[k][j]; 
    return C; 
}

double matrix::trace() { // Return the sum of the elements on the main diagonal
    assert(m == n); 
    double s = 0.0; 
    for(int i=1; i<=n; i++) s += data[i][i]; 
    return s; 
}

matrix matrix::tp() { // Return the transpose of A
    matrix At(n, m); 
    for(int i=1; i<=n; i++)
        for(int j=1; j<=m; j++)
            At.data[i][j] = data[j][i]; 
    return At; 
}

void matrix::write(ostream& out) { // Write to output stream
    out << rows() << " " << cols() << "\n"; 
    out << fixed << setprecision(2); 
    for(int i=1; i<=rows(); i++) {
        for(int j=1; j<=cols(); j++) {
            out << data[i][j]; 
            if(j < cols()) out << "\t"; 
        }
        if(i < rows()) out << "\n"; 
    }
    out << endl; 
}

void matrix::read(istream& in) { // Read from input stream
    for(int i=1; i<=rows(); i++)
        for(int j=1; j<=cols(); j++)
            in >> data[i][j]; 
    in >> ws; 
}
