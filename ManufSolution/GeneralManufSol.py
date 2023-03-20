from Common.MMSTags import VarSetTags, DefaultVarSetTags, SolutionTags, DefaultSolutionTags
from Common.Utilities import *

class GeneralManufSol:
    """ Abstract class defining the common functions required to assemble a manufactured solution.

    Attributes:
        sym_params: list of symbolic parameters that will be computed based on the conditions imposed on the manuf solution.
        sym_variables: list of symbolic independent variables like spatial coordinates or time.
    """
    def __init__(self, sym_params: list, sym_variables: list):
        """ Inits attributes of GeneralManufSol class. """
        self.sym_params = sym_params
        # self.num_params = num_p
        self.sym_variables = sym_variables
        self.sym_sol = None

    def getSymIndependentVar(self) -> list:
        """ Returns:
            The list of symbolic independent variables associated to the manuf sol. 
        """
        return self.sym_variables

    def getSymUnknownParam(self) -> list:
        """ Returns:
            The list of symbolic unknown parameters associated to the manuf sol. 
        """
        return self.sym_params

    def getSymSol(self):
        return self.sym_sol

    def buildSymSol(self):
        """ Returns:
            The symbolic expression of the specific manuf solution shape selected. 
            Has to be implemented in derived classes. """
        raise NotImplementedError()

    def setSymSol(self, provided_sym_sol):
        self.sym_sol = provided_sym_sol

    # def symbolicGradSol(self):
    #     raise NotImplementedError()

    # def derivSymSol(self, param_to_derive_wrt, deriv_order):
    #     return sp.diff(self.symbolicSol(),param_to_derive_wrt,deriv_order)

def assembleManufSolListIntoTagDic(manuf_sol_list: list, tags: list[SolutionTags]) -> dict:
    """ TODO """
    nb_fields = len(manuf_sol_list)
    if nb_fields != len(tags):
        raise ValueError("")
    manuf_sol_dic_per_tag = dict()
    for tag, isol in zip(tags,range(nb_fields)):
        if tag in manuf_sol_dic_per_tag:
            manuf_sol_dic_per_tag[tag].append(manuf_sol_list[isol])
        else:
            manuf_sol_dic_per_tag[tag] = manuf_sol_list[isol]
    return manuf_sol_dic_per_tag
