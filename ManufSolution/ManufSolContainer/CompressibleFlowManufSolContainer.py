from Common.MMSTags import FluidSolutionTags as fst
from Common.MMSTags import CompressibleFlowVarSetTags, VarSetTags, SolutionTags
from ..GeneralManufSol import GeneralManufSol 
from ManufSolution.ManufSolContainer.GeneralManufSolContainer import GeneralManufSolContainerSinglePhase

class CompressibleFlowManufSolContainer(GeneralManufSolContainerSinglePhase):
    """ TODO """
    def __init__(self, manuf_sol_lists: list[list[GeneralManufSol]], dim: int, spatial_sym_var: list, var_set: CompressibleFlowVarSetTags, temporal_sym_var: list = []):
        """ TODO """
        super().__init__(manuf_sol_lists, dim, spatial_sym_var, var_set, temporal_sym_var)
        if not isinstance(var_set, CompressibleFlowVarSetTags):
            raise ValueError("The variable sets choice must be of type 'CompressibleFlowVarSetTags'.")

    def getSolsTagsFromVarSetChoice(self, var_set_choice: VarSetTags) -> list[SolutionTags]:
        if var_set_choice == CompressibleFlowVarSetTags.CONSERVATIVE:
            if self.dim == 2:
                # return ["rho", "rhou", "rhov", "rhoE"]
                return [fst.DENSITY, fst.MOMENTUM_X, fst.MOMENTUM_Y, fst.RHOTOTALENERGY]
            elif self.dim == 3:
                # return ["rho", "rhou", "rhov", "rhow", "rhoE"]
                return [fst.DENSITY, fst.MOMENTUM_X, fst.MOMENTUM_Y, fst.MOMENTUM_Z, fst.RHOTOTALENERGY]
        elif var_set_choice == CompressibleFlowVarSetTags.PRIMITIVE_PVT:
            if self.dim == 2:
                # return ["P", "u", "v", "T"]
                return [fst.PRESSURE, fst.VELOCITY_X, fst.VELOCITY_Y, fst.TEMPERATURE]
            elif self.dim == 3:
                # return ["P", "u", "v", "w", "T"]
                return [fst.PRESSURE, fst.VELOCITY_X, fst.VELOCITY_Y, fst.VELOCITY_Z, fst.TEMPERATURE]
        elif var_set_choice == CompressibleFlowVarSetTags.PRIMITIVE_PVRHO:
            if self.dim == 2:
                # return ["P", "u", "v", "T"]
                return [fst.PRESSURE, fst.VELOCITY_X, fst.VELOCITY_Y, fst.DENSITY]
            elif self.dim == 3:
                # return ["P", "u", "v", "w", "T"]
                return [fst.PRESSURE, fst.VELOCITY_X, fst.VELOCITY_Y, fst.VELOCITY_Z, fst.DENSITY]
        else:
            raise ValueError("This type of compressible flow variable set is not available.")
