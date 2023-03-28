from enum import Enum

#* --- TAGS FOR SETS OF VARIABLES --- *#
# - Abstract level
class VarSetTags(Enum):
    pass

# - By-default variables set
class DefaultVarSetTags(VarSetTags):
    MYVARSETTAG = 1

# - Tags for Navier-Stokes fluid systems
class CompressibleFlowVarSetTags(VarSetTags):
    CONSERVATIVE = 1
    PRIMITIVE_PVT = 2
    PRIMITIVE_PVRHO = 3


#* --- TAGS FOR SOLUTION VARIABLES --- *#
# - Abstract level
class SolutionTags(Enum):
    pass

# - By-default solution tags
class DefaultSolutionTags(SolutionTags):
    VAR_1 = 1
    VAR_2 = 2
    VAR_3 = 3
    VAR_4 = 4
    VAR_5 = 5
    VAR_6 = 6
    VAR_7 = 7
    VAR_8 = 8
    VAR_9 = 9
    VAR_10 = 10

# - Tags for Navier-Stokes fluid systems
class FluidSolutionTags(SolutionTags):
    DENSITY = 1
    PRESSURE = 2
    TEMPERATURE = 3
    VELOCITY_X = 4
    VELOCITY_Y = 5
    VELOCITY_Z = 6
    VELOCITY_VEC = 7
    MOMENTUM_X = 8
    MOMENTUM_Y = 9
    MOMENTUM_Z = 10
    MOMENTUM_VEC = 11
    INTERNALENERGY = 12
    KINETICENERGY = 13
    TOTALENERGY = 14
    RHOTOTALENERGY = 15
    ENTHALPY = 16
    TOTALENTHALPY = 17

#* --- TAGS FOR GRADIENT VARIABLES --- *#
# - Abstract level
class SolutionGradientTags(Enum):
    pass

# # - By-default solution tags
# class DefaultGradientSolutionTags(SolutionGradientTags):
#     GRADVAR_1 = 1
#     GRADVAR_2 = 2
#     GRADVAR_3 = 3
#     GRADVAR_4 = 4
#     GRADVAR_5 = 5
#     GRADVAR_6 = 6
#     GRADVAR_7 = 7
#     GRADVAR_8 = 8
#     GRADVAR_9 = 9
#     GRADVAR_10 = 10

# - Tags for Navier-Stokes fluid systems
class FluidSolutionGradientTags(SolutionGradientTags):
    # DUDX = 1
    # DUDY = 2
    # DUDZ = 3
    # DVDX = 4
    # DVDY = 5
    # DVDZ = 6
    # DWDX = 7
    # DWDY = 8
    # DWDZ = 9
    GRADVEL = 10
    VORTICITY = 11
    # DTDX = 12
    # DTDY = 13
    # DTDZ = 14
    GRADT = 15
    SHEARSTRESS = 16
    HEATFLUX = 17
    DIVERGENCEVEL = 18
    VISCOUSDISSIPATION = 19
    HEATFLUX_MINUS_VISCOUSDISSIPATION = 20
    SHEARSTRESS_XX = 21
    SHEARSTRESS_XY = 22
    SHEARSTRESS_YY = 23
    HEATFLUX_X = 24
    HEATFLUX_Y = 25

#* --- TAGS FOR COORDINATES SYSTEMS TYPES --- *#
class CoordinatesSystemType(Enum):
    CARTESIAN = 1
    CYLINDRICAL = 2

#* --- TAGS FOR PROJECTION TYPES --- *#
class ProjectionType(Enum):
    NOPROJECTION = 1
    NORMAL = 2
    TANGENT = 3
    NORMALNORMAL = 4
    NORMALTANGENT = 5

#* --- TAGS FOR OUTPUT FILE TYPE --- *#
class OutputFileType(Enum):
    TEXT = 1
    CSV = 2
