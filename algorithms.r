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


library(tseriesChaos)
library(TSDecomposition)
library(Rssa)
library(OpenImageR)
library(forecast)
source('utils.r')

fourierDec <- function(series, params){
  freq.cutoff = unlist(params[1])

  coeffs = fft(series)
  mags   = coeffs[1:(length(coeffs)/2)]
  mags   = 1+sqrt(Re(mags)^2+Im(mags)^2)
  o.idx  = order(mags, decreasing = T)
  idx    = (1:length(mags))[-o.idx[1:freq.cutoff]]

  coeffs[idx] = complex(real=0, imaginary=0)
  coeffs[length(coeffs) - idx + 1] = complex(real=0, imaginary=0)
  det = Re(fft(coeffs, inverse=T)) / length(series)

  return(det)
}

waveletDec <- function(series, params){
    filter   = unlist(params[1])
    n.levels = as.numeric(unlist(params[2]))
    r.wavelet = wavelets::dwt(series,
                            filter = filter,
                            n.levels = n.levels,
                            fast=TRUE)
    for (i in 1:length(r.wavelet@W)) {
        r.wavelet@W[[i]] = cbind(rep(0, length(r.wavelet@W[[i]])))
    }
    det = wavelets::idwt(r.wavelet)
    return(det)
}

ssaDec <- function(series, par){

  L        = as.numeric(unlist(par[1]))
  neig     = as.numeric(unlist(par[2]))

  #execute
  s = ssa(series, L=L, neig=neig)#, kind=kind)
  r = reconstruct(s, groups = seq(1:(L/2)))

  #mutual information to separate deterministic components
  mi = c()
  for(i in 1:(length(r)-1))
    if(!all(r[[i]] == 0))
        mi = c(mi, FNN::mutinfo(r[[i]], r[[i+1]]))

  det.idx = 1
  if(length(mi) > 1) det.idx = which.max(abs(diff(mi))) + 1

  #sum det. comp.
  detComp = rep(0, length(r[[1]]))
  for(i in 1:det.idx)
    detComp = detComp + r[[i]]

  return(detComp)
}
      
emdrpDec   <- function(series, par){
    det = tryCatch({
          detlevel = unlist(par[1])
          thresh   = unlist(par[2])
          delay    = unlist(par[3])
          embedded = unlist(par[4])
          emdrp    = rpemdDecomposition(series, detlevel, thresh, delay, embedded)
          emdrp@deterministic
    }, error = function(e) {
        emdrpDec(series, par)
    })
    return(det)
}
      
emdmiDec <- function(series, par){
    source("emd-mi.r")
    return(DVEMD(series)$deterministic)
}

lazyDec <- function(series, par){
    l = length(series)
    x = 0
    c = 1
    m = unlist(par[3])
    d = unlist(par[4])
    i = unlist(par[2])
    r = unlist(par[1])
    v = 0

    ifile = tempfile()
    ofile = tempfile()
    write.table(series, file = ifile, row.names = F)
    nrlazy.exec = '/usr/local/Tisean_3.0.1/bin/nrlazy'
    shell.command = sprintf("'%s' %s -i%d -x%d -c%d -m1,%d -d%d -r%f -v0 -l%d -V0 > %s", 
                            nrlazy.exec, ifile, i, x, c, m, d, r, l, ofile)
    system(shell.command)
    det = unlist(read.table(ofile)[,1])
    unlink(ifile)
    unlink(ofile)
    return(det)
}
   
springDec <- function(series, params){
  frac   = unlist(params[1])
  num.it = unlist(params[2])
  m      = unlist(params[3])
  d      = unlist(params[4])

  s.emb = tseriesChaos::embedd(series, m=m, d=d)
  dist.mat = as.matrix(dist(s.emb))

  for(j in 1:num.it){
    s.emb  = t(mapply(function(i, y, s.emb, frac){
      r  = as.numeric(quantile(y[-i], probs = frac))
      nn = s.emb[which(y <= r),]
      dp = apply(nn, 2, sd)
      np = colMeans(nn)
    }, 1:nrow(s.emb), as.data.frame(dist.mat), MoreArgs=list(frac=frac, s.emb=s.emb)))
  }

  return(toTimeSpace(s.emb,m,d))
}