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

source('utils.r')
source('gridSearch.r')
source('algorithms.r')

dataFolder     = 'data'
pipelineFolder = 'pipeline'
seriesList     = loadSeriesFile(dataFolder)

params = expand.grid(
  L = 130:180,
  neig = 1:2
)
runPipeline(seriesList, ssaDec, params, pipelineFolder, 'ssa', run_r = TRUE)
