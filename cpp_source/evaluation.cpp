/*
-------------------------------------------------------------------------
  Time Series Decomposition Using Spring System Applied on Phase Spaces  
                                                                         
 Copyright (c) 2016-2019 Felipe S. L. G. Duarte, Ricardo A. Rios,        
 Eduardo R. Hruschka and Rodrigo F. de Mello, Sao Carlos/SP, Brazil.     
 All Rights Reserved.                                                    
                                                                         
 you can redistribute it and/or modify it under the terms of the GNU     
 General Public License as published by the Free Software Foundation,    
 either version 3 of the License, or (at your option) any later version. 
                                                                         
 Spring is distributed in the hope that it will be useful, but WITHOUT   
 ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or   
 FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License    
 for more details.                                                       
                                                                         
 Contributor(s):                                                         
 * Felipe S. L. G. Duarte - felipe.duarte@itau-unibanco.com.br           
                            fgduarte@icmc.usp.br                         
 * Ricardo A. Rios - ricardoar@ufba.br                                   
 * Eduardo R. Hruschka - eduardo.hruschka@itau-unibanco.com.br           
                         edu@poli                                        
 * Rodrigo F. de Mello - mello@icmc.usp.br                               
                                                                         
 You should have received a copy of the GNU General Public License along 
 with Spring. If not, see <http://www.gnu.org/licenses/>.                
                                                                         
 based on the publication:                                               
                                                                         
  @Article{spring2018,                                                   
   Title    = {Time Series Decomposition Using Spring System Applied on  
               Phase Spaces},                                            
   Author   = {Felipe S. L. G. Duarte, Ricardo A. Rios, Eduardo R.       
               Hruschka and Rodrigo F. de Mello},                        
   Journal  = {},                                                        
   Year     = {2018},                                                    
   Month    = {10},                                                      
   Number   = {},                                                        
   Pages    = {},                                                        
   Volume   = {},                                                        
   ISSN     = {},                                                        
   Doi      = {},                                                        
 }                                                                       
                                                                         
 The software is provided "As is", without warranty of any kind, express 
 or implied, including but not limited to the warranties of              
 merchantability, fitness for a particular purpose and noninfringement.  
 In no event shall the authors or copyright holders be liable for any    
 claim, damages or other liability, whether in an action of contract,    
 tort or otherwise, arising from, out of or in connection with the       
 software or the use or other dealings in the software.                  
-------------------------------------------------------------------------
*/

#include <algorithm>
#include <Rcpp.h>
#include <math.h>

using namespace Rcpp;

bool find_string(std::string str, Rcpp::List l){
    for(int i = 0; i < l.size(); ++i){
        if(str.compare(Rcpp::as<std::string>(l[i])) == 0){ 
            return true;
        }
    }
    return false;
}

double var(Rcpp::NumericVector x){
    return (sum(pow(x - mean(x),2)))/x.size();
}

double stdev(Rcpp::NumericVector x){
    return sqrt(var(x));
}

double md_dist(Rcpp::NumericMatrix nm1, Rcpp::NumericMatrix nm2){
    if(nm1.size() != nm2.size()) {
        Rcerr << "Could not execute Multidimension Distance. Matrix's has different sizes.\n";
    }

    double dist = 0.0;
    for(int i=0; i<nm1.nrow(); ++i){
        dist += sqrt(sum(pow(nm1(i,_) - nm2(i,_), 2)));
    }

    return dist/nm1.nrow();
}

Rcpp::NumericMatrix embedd(Rcpp::NumericVector ts, int m, int d){
    int nrow = ts.size() - ((m-1)*d);
    int ncol = m;
    NumericMatrix at(nrow, ncol);
    for(int i=0; i<ncol; i++){
        for(int j=0; j<nrow; j++){
            at(j,i) = ts[j+(i*d)];
        }
    }
    return at;
}

int whichMin( double x, double y, double z){
    if(x < y && x < z) return 1;
    if(y < x && y < z) return 2;
    return 3;
}


int whichMin(Rcpp::NumericVector dfo){
    int min_idx = 0;
    int min = dfo[0];
    for(int i = 0; i < dfo.size(); ++i){
        if(dfo[i] < min){
            min = dfo[i];
            min_idx = i;
        }
    }
    return min_idx;
}


double min(double x, double y, double z){
    if(x < y && x < z) return x;
    if(y < x && y < z) return y;
    return z;
}

// [[Rcpp::export]]
Rcpp::List mddl(Rcpp::NumericVector ts1, Rcpp::NumericVector ts2, Rcpp::Function dtw_f){

    Rcpp::List dtw = dtw_f(ts1, ts2);
    Rcpp::NumericVector x = dtw["index1"];
    Rcpp::NumericVector y = dtw["index2"];
    int sizeTS = ts1.size();
    int sizeMDDL = x.size();

    //double mddl (double *x, double *y, int sizeTS, int sizeMDDL){

    double result = 0;
    double step = ((double)sizeTS-1)/((double)sizeMDDL-1);
    double interval = 1;

    for(int i = 0; i < sizeMDDL; i++){
        result += sqrt(pow(x[i] - interval, 2) + pow(y[i] - interval,2));
        interval += step;
    }
    double val = result/sizeMDDL;
    return List::create(_["name"]="MDDL", _["best"]="min", _["val"]=val);
}

// [[Rcpp::export]]
Rcpp::List mda(Rcpp::NumericVector ts1, Rcpp::NumericVector ts2, int m, int d){
    NumericMatrix at1 = embedd(ts1, m, d);
    NumericMatrix at2 = embedd(ts2, m, d);
    //int nrow = at1.nrow();
    //int delta = at1.nrow() - at2.nrow();
    //if(delta == 0) 
    double val = md_dist(at1, at2);
    return List::create(_["name"]="MDA", _["best"]="min", _["val"]=val);
}

// [[Rcpp::export]]
Rcpp::NumericVector whichBest(Rcpp::NumericMatrix tbl){
    double mddl_max = max(tbl(_,1));
    double mda_max  = max(tbl(_,2));
    Rcpp::NumericVector mddl = tbl(_,1)/( (mddl_max < 0.0000001)? 1.0 : mddl_max );
    Rcpp::NumericVector mda  = tbl(_,2)/( (mda_max  < 0.0000001)? 1.0 : mda_max  );
    Rcpp::NumericVector dfo  = sqrt(pow(mddl,2)+pow(mda,2));
    int idx = whichMin(dfo);
    Rcpp::NumericVector result = tbl(idx,_);
    result.push_back(dfo[idx]);

    return NumericVector::create(_["idx"]    = idx + 1,
                                 _["mddl"]   = tbl(idx,1),
                                 _["mda"]    = tbl(idx,2),
                                 _["n_mddl"] = mddl[idx],
                                 _["n_mda"]  = mda[idx],
                                 _["dfo"]    = dfo[idx]
                                );
}

// [[Rcpp::export]]
Rcpp::DataFrame evaluation_table_cpp(Rcpp::NumericVector ts1, Rcpp::NumericVector ts2, int m, int d, Rcpp::Function dtw_f){
    //ts1 always the biggest time series
    if(ts1.size() < ts2.size()){
        NumericVector ts_aux = ts1;
        ts1 = ts2;
        ts2 = ts_aux;
    }
    List eval_result(0);
    eval_result.push_back(mddl(ts1, ts2, dtw_f));
    eval_result.push_back(mda(ts1, ts2, m, d));
    //cria DataFrame to return
    StringVector  name(eval_result.size());
    StringVector  best(eval_result.size());
    NumericVector val(eval_result.size());
    for(int i=0; i < eval_result.size(); ++i){
        List er = eval_result[i];
        name[i] = Rcpp::as<std::string>(er["name"]);
        best[i] = Rcpp::as<std::string>(er["best"]);
        val[i]  = er["val"];
    }

    return DataFrame::create(_["name"]=name, _["best"]=best, _["val"]=val, _["stringsAsFactors"] = false );
}

// [[Rcpp::export]]
Rcpp::List evaluate_cpp(Rcpp::List result, Rcpp::Function dtw_f){
    Rcpp::List gs = result["gs"];
    Rcpp::List fgs = gs[0];
    Rcpp::List fgsp = fgs["params"];
    Rcpp::List series_obj = result["series_obj"];
    int m = series_obj["det.embDim"];
    int d = series_obj["det.sepDim"];
    Rcpp::NumericVector det = series_obj["det.series"];
    Rcpp::NumericMatrix tbl(gs.size(), 3);
    Rcpp::CharacterVector scale(2);
    for(int i = 0; i < gs.size(); ++i){
        List gs_i = gs[i];
        Rcpp::DataFrame df = evaluation_table_cpp(gs_i["det"], det, m, d, dtw_f);
        scale = df["best"];
        CharacterVector name = df["name"];
        name.push_front("test_id");
        colnames(tbl) = name;
        NumericVector laux = df["val"];
        tbl(i, 0) = i; //add test id
        for (int j = 1; j <= laux.size(); ++j){
            tbl(i, j) = laux[j-1];
        }
    }

    return Rcpp::List::create(_["tbl"] = tbl, _["scale"] = scale);
}
