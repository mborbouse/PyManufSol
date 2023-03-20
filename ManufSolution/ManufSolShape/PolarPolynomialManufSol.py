from ..GeneralManufSol import GeneralManufSol

class PolarPolynomialManufSol(GeneralManufSol):
    """ Class to build polynomial manuf sol in polar coordinates of the form 
            c_0 + c_1 * r + c_2 * r^2 + c_3 * r^3 + ...

    Attributes:
        num_params: numeric values of the known coefficients in the polynomial
        num_params_indices: indices of those coefficients w.r.t. the form given above 
        highest_order: order of the polynomial, i.e. order of the last term here
    """

    def __init__(self, sym_params, sym_variables, num_params, num_params_indices, highest_order):
        """ Inits attributes of PolarPolynomialManufSol class. """
        super().__init__(sym_params, sym_variables)
        self.num_params = num_params
        self.num_params_indices = num_params_indices
        self.highest_order = highest_order
        if len(num_params) != len(num_params_indices):
            raise ValueError("Number of numerical parameters (%d) != number of their indices (%d)"%(len(num_params),len(num_params_indices)))
        
        self.number_term = self.highest_order+1

        sym_params_indices = list(range(0, self.number_term))
        self.sym_params_indices = [i for i in sym_params_indices if i not in self.num_params_indices]
        self.buildSymSol()

    def buildSymSol(self):
        """ Returns the polynomial manuf sol in the coordinate contained in 'sym_variables[0]'. """
        num_part = 0
        j = 0
        for i in self.num_params_indices:
            num_part = num_part + self.num_params[j]*pow(self.sym_variables[0],i)
            j = j+1

        sym_part = 0
        j = 0
        for i in self.sym_params_indices:
            sym_part = sym_part + self.sym_params[j]*pow(self.sym_variables[0],i)
            j = j+1

        self.sym_sol = num_part + sym_part


