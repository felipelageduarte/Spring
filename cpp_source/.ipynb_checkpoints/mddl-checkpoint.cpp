#include "mddl.hpp"
#include <math.h>

int mddl::whichMin( double x, double y, double z){
    if(x < y && x < z) return 1;
    if(y < x && y < z) return 2;
    return 3;
}

matrix mddl::cost_matrix(vector<double> x, vector<double> y){
	matrix costMat(x.size(), y.size());
    double cost;
	for (int i = 0; i < x.size(); ++i){
		for (int j = 0; j < y.size(); ++j){
            cost = sqrt(pow(x[i] - y[j], 2));
            if(i == 0 && j == 0){
                costMat[i][j] = cost;
            } else if(i == 0){
                costMat[i][j] = cost + costMat.get(i, j-1);
            } else if(j == 0){
                costMat[i][j] = cost + costMat.get(i-1, j);
            } else {
                costMat[i][j] = min(min(costMat.get(i-1, j) + cost, costMat.get(i  , j-1) + cost), costMat.get(i-1, j-1) + 2 * cost);
            }
		}
	}

	return costMat;
}

matrix mddl::align(matrix costMat){
    matrix index(max(costMat.rows(), costMat.cols()), 2);
    
    int i = 0;
    int j = 0;
    int c = 1;
    
    index[c][0] = i;
    index[c][1] = j;
    c++;

    while(1){
        if(i == costMat.rows() && j == costMat.cols()){
            break;
        } else if(i == costMat.rows()) {
            j++;
        } else if(j == costMat.cols()) {
            i++;
        } else {
            int idx = whichMin(costMat.get(i,j+1), costMat.get(i+1,j+1), costMat.get(i+1,j));
            if(idx == 1){
                j++;
            } else if(idx == 2){
                j++;
                i++;
            } else if(idx == 3){
                i++;
            } 
        }

        index[c][0] = i;
        index[c][1] = j;
        c++;
    }
    
    return index;
}

//mddl deve receber os indices 1(x) e 2(y) como parametro da funcao]
//sizeTS e sizeMDDL representam o tamanho das series comparadas (supondo 
//que sao do mesmo tamanho) e o tamanho dos warping paths
double
mddl (double *x, double *y, int sizeTS, int sizeMDDL){
	double result = 0;
	double step = ((double)sizeTS-1)/((double)sizeMDDL-1);
	double interval = 1;

	for(int i = 0; i < sizeMDDL; i++){
		result += sqrt(pow(x[i] - interval, 2) + pow(y[i] - interval,2));
		interval += step;
	}

	return result/sizeMDDL;
}

//versao normalizada recebe o valor da mddl e o tamanho das duas series
// -- supondo que as series sao do mesmo tamanho
double
normalizedMDDL (double mddl, int sizeTS){
	double result = 0;
	int j = 0, k = 0;
	double step = ((double)sizeTS-1)/((double)sizeTS*2-1);
	double interval = 1;

	for(int i = 0; i < (2*sizeTS); i++){
		if(i < sizeTS){
			result += sqrt(pow(1-interval,2) + pow(j+1-interval,2));
			j++;
		}else{
			result += sqrt(pow(k+1-interval,2) + pow(j-interval,2));
			k++;
		}
		interval += step;
	}

	return (mddl*sizeTS*2/result);
}


double mddl::evaluate(vector<double> x, vector<double> y){
    matrix costMat = cost_matrix(x,y);
    matrix index = align(costMat);
    return 0.0;
}