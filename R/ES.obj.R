# Copyright 2009 John R. Stevens
# Distributed as part of the metahdep package, under the terms of the GNU General Public License (see DESCRIPTION file)

setClass(Class="ES.obj",
         representation=representation(  
                    gn="character", 
                    ES.mat="matrix", 
                    Cov.mat="matrix", 
                    chip="character", 
                    covariates="data.frame", 
                    dep.grp="integer"
                         ),
          package="metahdep"
         )

