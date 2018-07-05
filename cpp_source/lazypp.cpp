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
