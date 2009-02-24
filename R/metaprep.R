# Copyright 2009 John R. Stevens
# Distributed as part of the metahdep package, under the terms of the GNU General Public License (see DESCRIPTION file)

setClass(Class="metaprep",
         representation=representation(  
                  theta="vector", 
                  V="matrix", 
                  X="matrix", 
                  M="matrix",
                  max.k="integer",
                  row.indices="matrix",
                  gene="character"
                         ),
         package="metahdep"
         )

