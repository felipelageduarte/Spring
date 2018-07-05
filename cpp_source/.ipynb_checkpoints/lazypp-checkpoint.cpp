#include "utils.hpp"
#include <Rcpp.h>

// [[Rcpp::export]]
Rcpp::NumericVector lazypp(Rcpp::NumericVector series, double frac, int m, int d, int numIt){
    
    Rcpp::NumericMatrix s_emb  = embedd(series, m, d);
    Rcpp::NumericMatrix dist_m = dist(s_emb);
    Rcpp::NumericVector row;
    Rcpp::NumericVector det(s_emb.nrow());
    int idx_frac = static_cast<int>(s_emb.nrow() * frac);
        
    for(int it = 0; it < numIt; ++it){
        for(int i = 0; i < s_emb.nrow(); ++i){
            row = dist_m(i,Rcpp::_);
            Rcpp::NumericVector o = orderVec(row);
            det(i) = o(idx_frac);
        }
    }
    
    return det;
/*   
    s.emb  = t(mapply(function(i, y, s.emb, frac){
      r  = as.numeric(quantile(y[-i], probs = frac))
      nn = s.emb[which(y <= r),]
      dp = apply(nn, 2, sd)
      np = colMeans(nn)
    }, 1:nrow(s.emb), as.data.frame(dist.mat), MoreArgs=list(frac=frac, s.emb=s.emb)))
    }

    return(toTimeSpace(s.emb,m,d))
    
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
        
        Rcpp::List l;
        l["series"] = series;
        l["params"] = (params.isNull()) ? 0 : getByRow(p, i);
        l["f"]      = f;
        l["det"]    = det;
        l["sto"]    = series - det;
        l["rt"]     = ((float)t2 - (float)t1) / CLOCKS_PER_SEC;

        gs[i] = l;
        
    }
    
    gsResult["gs"] = gs;

    return gsResult;
*/
}
