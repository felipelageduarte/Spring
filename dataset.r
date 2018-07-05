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

require(fNonlinear)
require(tseriesChaos)
require(nonlinearTseries)
source('utils.r')

seriesSize = 1000
dataFolder = 'data'

getDeterministicSeries <- function(comp, size = 1000){
  if(missing(comp)){
    stop("You must select one of those options for deterministic component: sine, lorenz, rossler, logistic or henon")
  }

  if(tolower(comp) == 'sine'){
    return(list( det.series = normalize(sin(2*pi*seq(0,9,len=size))),
                 det.sepDim = 20,
                 det.embDim = 2,
                 det.model  = "sine" ,
                 det.length = size)
           )
  } else if(tolower(comp) == 'lorenz') {
    return(
      list(det.series = normalize(lorentzSim(doplot = FALSE,
                            times  = seq(0, 50, by = (50/size)),
                            parms  = c(sigma = 16, r = 45.92, b = 4),
                            start  = c(-14, -13, 47)
      )[,2][1:size]),
      det.sepDim = 1,
      det.embDim = 3,
      det.model  = "lorenz",
      det.length = size)
    )
  } else if(tolower(comp) == 'rossler') {
      return(
        list(det.series = normalize(sim.cont(rossler.syst,
                            start.time=0,
                            end.time=650,
                            dt=650/size,
                            start.x=c(0,0,0),
                            parms=c(0.15, 0.2, 10))[1:size]),
             det.sepDim = 1,
             det.embDim = 3,
             det.model  = "rossler",
             det.length = size)
      )
  } else if(tolower(comp) == 'logistic') {
    return(
      list(det.series = normalize(logisticSim(n=size, parms = c(r = 3.8), start = 0.5, doplot = FALSE)),
           det.sepDim = 1,
           det.embDim = 2,
           det.model  = "logistic",
           det.length = size)
    )
  } else if(tolower(comp) == 'henon') {
    return(
      list(det.series = normalize(nonlinearTseries::henon(
                                start=c(-0.006423277,-0.473545134),
                                n.sample = size,
                                n.transient=10,
                                do.plot=FALSE)$x[1:size]),
           det.sepDim = 1,
           det.embDim = 2,
           det.model  = "henon",
           det.length = size)
    )
  } else {
    stop("You must select one of those options for deterministic component: sine, lorenz, rossler, logistic or henon")
  }
}

getStochasticSeries <- function(comp, params=list(), size = 1000){
  if(missing(comp)){
    stop("You must select one of those options for stochastic component: zero, uniforme or normal")
  }
    
  if(tolower(comp) == 'zero'){
    return(list( sto.series = rep(0, size),
                 sto.model  = "zero" ,
                 sto.params = list(),
                 sto.length = size)
    )
  } else if(tolower(comp) == 'uniforme') {
    min = ifelse(is.null(params$min), -1, params$min)
    max = ifelse(is.null(params$max),  1, params$max)
    return(list( sto.series = runif(size, min, max),
                 sto.model  = paste("uniforme_",min,"_",max, sep=""),
                 sto.params = list(min=min, max=max),
                 sto.length = size)
    )
  } else if(tolower(comp) == 'normal') {
    mean = ifelse(is.null(params$mean), 0, params$mean)
    sd   = ifelse(is.null(params$sd), 1, params$sd)
    return(list( sto.series = rnorm(size, mean=mean, sd=sd),
                 sto.model  = paste("normal_",mean,"_",sd, sep=""),
                 sto.params = list(mean=mean, sd=sd),
                 sto.length = size)
    )
  } else {
    stop("You must select one of those options for stochastic component: zero, uniforme or normal")
  }
}

timeSeriesFactor <- function(det.comp, sto.comp, sto.params=list(), size = 1000){
  det = getDeterministicSeries(det.comp, size)
  sto = getStochasticSeries(sto.comp, sto.params, size)
  tsObj = list(series = det$det.series + sto$sto.series, size = size)
  tsObj = c(tsObj, det, sto)
  return(tsObj)
}

det = list("sine", "lorenz", "rossler", "logistic", "henon")
sto = list(list(comp="zero"),
           list(comp="uniforme",params=list(min=-0.01, max=0.01)),
           list(comp="uniforme",params=list(min=-0.05, max=0.05)),
           list(comp="uniforme",params=list(min=-0.10, max=0.10)),
           list(comp="uniforme",params=list(min=-0.15, max=0.15)),
           list(comp="uniforme",params=list(min=-0.20, max=0.20)),
           list(comp="normal",params=list(sd=0.01)),
           list(comp="normal",params=list(sd=0.05)),
           list(comp="normal",params=list(sd=0.10)),
           list(comp="normal",params=list(sd=0.15)),
           list(comp="normal",params=list(sd=0.20))
          )

set.seed(42)
for(it in 1:30){
    idx = 1
    for(i in 1:length(det)){
        for(j in 1:length(sto)){
          seriesObj = timeSeriesFactor(det[[i]], sto[[j]]$comp, sto[[j]]$params, seriesSize)
          filename  = paste(dataFolder, sprintf('/series_%02d_%02d', it, idx),'.RData',sep='')
          save(seriesObj, file=filename)
          idx = idx + 1
        }
    }
}
