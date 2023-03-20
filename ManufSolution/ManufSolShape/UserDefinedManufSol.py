from ..GeneralManufSol import GeneralManufSol

class UserDefinedManufSol(GeneralManufSol):
    """ Class for user defined manuf solution, i.e. no explicit expression is available here, 
        it should be specified by the user outside of the class.

    Attributes:
        user_defined_sol: the symbolic solution constructed outside 
                      of the class by the user and passed as input to the class
    """
    def __init__(self, sym_params, sym_variables, user_defined_sol):
        """ Inits attributes of UserDefinedManufSol class. """
        super().__init__(sym_params, sym_variables)
        self.user_defined_sol = user_defined_sol
        self.buildSymSol()

    def buildSymSol(self):
        """ The symbolic solution is directly the one given in the constructor in this case since
            the user is supposed to build its own manuf sol outside of the class. """
        self.sym_sol = self.user_defined_sol