#include <Rcpp.h>
#include <time.h>

Rcpp::List getByRow(Rcpp::DataFrame df, int idx){
    Rcpp::List row(df.ncol());
    for(int i = 0; i < df.ncol(); ++i){
        Rcpp::List col = df[i];
        row[i] = col[idx];
    }
    return row;
}

Rcpp::NumericVector rep(double val, int size){
    Rcpp::NumericVector vec(size);
    for(int i=0; i < size; ++i) vec[i] = val;
    return vec;
}

// [[Rcpp::export]]
Rcpp::List gridSearch_cpp(Rcpp::List series_obj, Rcpp::Function f, Rcpp::Nullable<Rcpp::DataFrame> params){
    Rcpp::NumericVector series = series_obj["series"];
    Rcpp::NumericVector obsDet = series_obj["det.series"];
    Rcpp::NumericVector obsSto = series_obj["sto.series"];
    Rcpp::DataFrame p = Rcpp::as<Rcpp::DataFrame>(params);

    Rcpp::List gsResult;
    gsResult["series_obj"] = series_obj;
    gsResult["gs"] = Rcpp::List(0);

    Rcpp::List gs((p.nrow() > 1)?p.nrow():1);
    
    for (int i = 0; i < p.nrow() || i == 0; ++i){
        clock_t t1 = clock();
        Rcpp::NumericVector det = f(series, (params.isNull()) ? 0 : getByRow(p, i));
        clock_t t2 = clock();
        
        gs[i] = Rcpp::List::create(Rcpp::_["series"] = series,
                                   Rcpp::_["params"] = (params.isNull()) ? 0 : getByRow(p, i),
                                   Rcpp::_["f"]      = f,
                                   Rcpp::_["det"]    = det,
                                   Rcpp::_["sto"]    = series - det,
                                   Rcpp::_["rt"]     = ((float)t2 - (float)t1) / CLOCKS_PER_SEC
                                  );        
    }
    
    gsResult["gs"] = gs;

    return gsResult;
}
