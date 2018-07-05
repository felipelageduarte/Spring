#-------------------------------------------------------------------------#
#  Time Series Decomposition Using Spring System Applied on Phase Spaces  #
#                                                                         #
# Copyright (c) 2016-2019 Felipe S. L. G. Duarte, Ricardo A. Rios,        #
# Eduardo R. Hruschka and Rodrigo F. de Mello, Sao Carlos/SP, Brazil.     #
# All Rights Reserved.                                                    #
#                                                                         #
# you can redistribute it and/or modify it under the terms of the GNU     #
# General Public License as published by the Free Software Foundation,    #
# either version 3 of the License, or (at your option) any later version. #
#                                                                         #
# Spring is distributed in the hope that it will be useful, but WITHOUT   #
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or   #
# FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License    #
# for more details.                                                       #
#                                                                         #
# Contributor(s):                                                         #
# * Felipe S. L. G. Duarte - felipe.duarte@itau-unibanco.com.br           #
#                            fgduarte@icmc.usp.br                         #
# * Ricardo A. Rios - ricardoar@ufba.br                                   #
# * Eduardo R. Hruschka - eduardo.hruschka@itau-unibanco.com.br           #
#                         edu@poli                                        #
# * Rodrigo F. de Mello - mello@icmc.usp.br                               #
#                                                                         #                
# You should have received a copy of the GNU General Public License along #
# with Spring. If not, see <http://www.gnu.org/licenses/>.                #
#                                                                         #
# based on the publication:                                               #
#                                                                         #
#  @Article{spring2018,                                                   #
#   Title    = {Time Series Decomposition Using Spring System Applied on  #
#               Phase Spaces},                                            #
#   Author   = {Felipe S. L. G. Duarte, Ricardo A. Rios, Eduardo R.       #
#               Hruschka and Rodrigo F. de Mello},                        #
#   Journal  = {},                                                        # 
#   Year     = {2018},                                                    #
#   Month    = {10},                                                      #
#   Number   = {},                                                        #
#   Pages    = {},                                                        #
#   Volume   = {},                                                        # 
#   ISSN     = {},                                                        #  
#   Doi      = {},                                                        #  
# }                                                                       #
#                                                                         #
# The software is provided "As is", without warranty of any kind, express #
# or implied, including but not limited to the warranties of              #
# merchantability, fitness for a particular purpose and noninfringement.  #
# In no event shall the authors or copyright holders be liable for any    #
# claim, damages or other liability, whether in an action of contract,    #
# tort or otherwise, arising from, out of or in connection with the       #
# software or the use or other dealings in the software.                  #
#-------------------------------------------------------------------------#

library(tictoc)
library(Rcpp)

loadSeriesFile <- function(seriesFolder){
  #select all series file in the data folder
  seriesFile = list.files(path = seriesFolder, full.names = TRUE)
  seriesList = list()
  for(i in 1:length(seriesFile))
    seriesList[[i]] = get(load(seriesFile[i]))
  return(seriesList)
}

toTimeSpace = function(s.emb, m, d){
    if(m <= 1 | ncol(s.emb) != m){ stop('Erro') }
    n  = nrow(s.emb)
    k  = ceiling(m/2)
    ts = c()
    if(m == 2){
        ts = c(s.emb[,1], s.emb[(n-d+1):n,2])
    } else {
        for(i in 1:max(c(1,k-1))) ts = c(ts, s.emb[1:d, i])
        ts = c(ts, s.emb[, k])
        for(i in (k+1):m) ts = c(ts, s.emb[(n-d+1):n,i])
    }
    return(ts)
}

df_to_string <- function(id, df){
    sprintf('%s) %s', id, paste(df, collapse=", "))
}
     
save_results <- function(result, evaluation, tech_name, folder){
    tic('saving')
    det.name = result[[1]]$series_obj$det.model
    save(result,     file=paste(folder,'/',det.name,'_',tech_name,'_result.RData', sep=''))
    save(evaluation, file=paste(folder,'/',det.name,'_',tech_name,'_evaluation.RData', sep=''))
    toc()
}


