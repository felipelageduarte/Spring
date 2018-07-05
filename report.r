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

library(data.table)

report <- function(result, evaluation, report.idx = 1){
    
    series_obj   = result$series_obj
    params       = rbindlist(lapply(result$gs, function(x) x$params))
    eval.metrics = c()
    metrics      = colnames(evaluation$tbl)
    eval.metrics = evaluation$tbl[,-1]
    metrics      = colnames(eval.metrics)
    
    tbl = data.frame(name   = character(),
                     best   = character(),
                     it     = integer(),
                     params = character(),
                     val    = double(), 
                     stringsAsFactors=FALSE)

    for(i in 1:length(metrics)){
    tbl[i,'name']   = metrics[i]
    tbl[i,'best']   = evaluation$scale[i]
    tbl[i,'it']     = ifelse(tbl[i,'best']=='min', which.min(evaluation$tbl[,i]), which.max(evaluation$tbl[,i]))
    tbl[i,'params'] = paste(params[tbl[i,'it'],], collapse=", ")
    tbl[i,'val']    = evaluation$tbl[ tbl[i,'it'] ,i]
}
    show(tbl)
    
    if(nrow(params) != 1){
        par(mfrow = c(1,3))
        for(i in 1:length(metrics)) hist(eval.metrics[,i], main=metrics[i], xlab = metrics[i], breaks = 100)
        par(mfrow = c(1,1))
    
        rt = sapply(result$gs, function(x) x$rt)
        ts.plot(rt, main='Running Time (s)', xlab='Running Time (s)')
        abline(h=mean(rt), col=2)

        data = cbind(params, eval.metrics)
        for(i in names(data)){
            if(is.character(data[[i]])){
                data[[i]] = as.numeric(as.factor(data[[i]]))
            }
        }

        for(i in 1:length(metrics)){
            if((ncol(params)) == 1){
                plot(x=unlist(params[,1]), y=unlist(eval.metrics[,i]), 
                     ylab = metrics[i], main = metrics[i], xlab=colnames(params)[1], type='b')
            } else {
                plot(wireframe(
                    as.formula(paste(c(metrics[i], paste(colnames(params)[1:2] ,collapse ='*') ), collapse='~')),
                    data = data,
                    drape = TRUE,
                    colorkey = TRUE
                ))
            } 
        }  
    }  

    m = result$series_obj$det.embDim
    d = result$series_obj$det.sepDim
    for(i in 1:nrow(tbl)){
        idx = as.numeric(tbl$it[i])
        series.det = result$gs[[idx]]$det
        plot(tseriesChaos::embedd(result$gs[[idx]]$series, m, d), main=paste('Det: best', tbl$name[i]))
        points(tseriesChaos::embedd(series.det, m, d), col='blue')
        ts.plot(result$gs[[idx]]$series, main=paste('Det: best',tbl$name[i]))
        lines(series.det, col='blue')
        acf(result$gs[[idx]]$sto, main=paste('Sto: best',tbl$name[i]))
    }   
}
