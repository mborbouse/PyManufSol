import sympy as sp

#* 1/ Create manuf solution for each phase
#   Inputs:
#       - independent variables
#       - possibly symbolic unknown parameters
#   Outputs:
#       - an object containing:
#           > symbolic solutions vectors depending on indep variables 
#             + possible unkwown parameters to be determined 
#             => either defined by the user directly or through the use of pre-determined
#             manuf sol shape classes
#           > the independent variables
#           > the unknown parameters


#* 2/ Create physical model for each phase
#   Inputs:
#       -  object for symbolic manuf solution for this phase
#   Outputs:
#       -  object containing:
#           > 
#             