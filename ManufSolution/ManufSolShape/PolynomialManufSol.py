from ..GeneralManufSol import GeneralManufSol

class PolynomialManufSol(GeneralManufSol):
    """ Class to build polynomial manuf sol of the form 
            c_0 + 
            c_1  * x^1 + c_2 * y^1 + 
            c_3  * x^2 * y^0 + c_4  * x^1 * y^1 + c_5  * x^0 * y^2 +
            c_6  * x^3 * y^0 + c_7  * x^2 * y^1 + c_8  * x^1 * y^2 + c_9  * x^0 * y^3 +
            c_10 * x^4 * y^0 + c_11 * x^3 * y^1 + c_12 * x^2 * y^2 + c_13 * x^1 * y^3 + c_14 * x^0 * y^4 
            ...

    Attributes:
        num_params: numeric values of the known coefficients in the polynomial
        num_params_indices: indices of those coefficients w.r.t. the form given above 
        highest_order: order of the polynomial, i.e. order of the last term here
    """

    def __init__(self, sym_params, sym_variables, num_params, num_params_indices, highest_order):
        """ Inits attributes of PolynomialManufSol class. """
        super().__init__(sym_params, sym_variables)
        self.num_params = num_params
        self.num_params_indices = num_params_indices
        self.highest_order = highest_order
        if len(num_params) != len(num_params_indices):
            raise ValueError("Number of numerical parameters (%d) != number of their indices (%d)"%(len(num_params),len(num_params_indices)))
        
        #! TODO
        #self.order_p = len(sym_ind)+len(num_ind)
        self.buildSymSol()

    def buildSymSol(self):
        """ !TODO! """
        raise NotImplementedError()

