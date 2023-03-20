from PhysicalModels.GeneralPhysicalModels import GeneralPhysicalModels
from ManufSolution.ManufSolContainer.GeneralManufSolContainer import GeneralManufSolContainerSinglePhase
from Common.MMSTags import CoordinatesSystemType
from Common.Utilities import *
from typing import Optional

class LinearAdvectionDiffusionModels(GeneralPhysicalModels):
    """ 
    Convective fluxes:
    F_x^c = [ [v_11_x v_12_x v_13_x]
              [v_21_x v_22_x v_23_x]    = V_x (nbFiels x nbFields)
              [v_31_x v_32_x v_33_x] ] x 

              [u_1
               u_2
               u_3]

    F_y^c = [ [v_11_y v_12_y v_13_y]
              [v_21_y v_22_y v_23_y]    = V_y (nbFiels x nbFields)
              [v_31_y v_32_y v_33_y] ] x 

              [u_1
               u_2
               u_3]

    ===> convTensor = [V_x, V_y]

    Diffusive fluxes:
    F_x^d = [ [k^x_11_x k^x_12_x k^x_13_x]
              [k^x_21_x k^x_22_x k^x_23_x]   = K^x_x (nbFiels x nbFields)
              [k^x_31_x k^x_32_x k^x_33_x] ] x  
    
              [d(u_1)/dx
               d(u_2)/dx
               d(u_3)/dx]

          + [ [k^y_11_x k^y_12_x k^y_13_x]
              [k^y_21_x k^y_22_x k^y_23_x]   = K^y_x (nbFiels x nbFields)
              [k^y_31_x k^y_32_x k^y_33_x] ] x  
    
              [d(u_1)/dy
               d(u_2)/dy
               d(u_3)/dy]

    F_y^d = [ [k^x_11_y k^x_12_y k^x_13_y]
              [k^x_21_y k^x_22_y k^x_23_y]   = K^x_y (nbFiels x nbFields)
              [k^x_31_y k^x_32_y k^x_33_y] ] x  
    
              [d(u_1)/dx
               d(u_2)/dx
               d(u_3)/dx]

          + [ [k^y_11_y k^y_12_y k^y_13_y]
              [k^y_21_y k^y_22_y k^y_23_y]   = K^y_y (nbFiels x nbFields)
              [k^y_31_y k^y_32_y k^y_33_y] ] x  
    
              [d(u_1)/dy
               d(u_2)/dy
               d(u_3)/dy]

    ===> diffTensor = [K^x_x, K^y_x, K^x_y, K^y_y]
    """
    def __init__(self, tag: str, manuf_sol_container: GeneralManufSolContainerSinglePhase, domain_dim: int, convTensor: list[list[list[sp.Expr]]], diffTensor: list[list[list[sp.Expr]]], coord_system: CoordinatesSystemType = CoordinatesSystemType.CARTESIAN):
        super().__init__(tag, manuf_sol_container, domain_dim, coord_system)
        self.convTensor = convTensor
        self.diffTensor = diffTensor

    def convectiveFlux(self, num_params: Optional[dict] = None):
        F_c = initListOfLists(self.manuf_sol_container.getNumberFields(), self.domain_dim)
        sol = self.solVector(num_params)
        for d in range(self.domain_dim):
            for f in range(self.manuf_sol_container.getNumberFields()):
                val = [a*b for a,b in zip(self.convTensor[d][f],sol)]
                F_c[f][d] = sum(val)
        return F_c

    def diffusiveFlux(self, num_params: Optional[dict] = None, replace_params_before_deriv: bool = False):
        F_d = initListOfLists(self.manuf_sol_container.getNumberFields(), self.domain_dim)
        gradSol = self.gradSolVector(num_params, replace_params_before_deriv)
        nbFields = self.manuf_sol_container.getNumberFields()
        # for d in range(self.domain_dim):
        #     for f1 in range(nbFields):
        #         for f2 in range(nbFields):
        #             val = [a*b for a,b in zip(self.diffTensor[f1+d*nbFields][f2],gradSol[f1])]
        #             F_d[f2][d] = F_d[f2][d] + sum(val)
        for d1 in range(self.domain_dim):
            for d2 in range(self.domain_dim):
                for f1 in range(nbFields):
                    coeff_vec = self.diffTensor[d1+d2][f1]
                    grad_vec = gradSol[f1]
                    mult = 0.0
                    for f2 in range(nbFields):
                        mult = mult + coeff_vec[f2]*grad_vec[f2][d2]
                    F_d[f1][d1] = F_d[f1][d1] + mult
        return F_d

    def sourceTerm(self, num_params: Optional[dict] = None, replace_params_before_deriv: bool = False):
        S = initListOfLists(self.manuf_sol_container.getNumberFields(), 1, 0.0)
        return S


class LinearScalarAdvectionDiffusionModels(LinearAdvectionDiffusionModels):
    """
    Convective fluxes:
    F_x^c = v_x x [u_1
                   u_2
                   u_3]

    F_y^c = v_y x [u_1
                   u_2
                   u_3]

    ===> convVel = [v_x, v_y]

    Diffusive fluxes:
    F_x^d = k_x x [d(u_1)/dx
                   d(u_2)/dx
                   d(u_3)/dx]

    F_y^d = k_y x [d(u_1)/dy
                   d(u_2)/dy
                   d(u_3)/dy]

    ===> diffCoeff = [k_x, k_y]

    """
    def __init__(self, tag: str, manuf_sol_container: GeneralManufSolContainerSinglePhase, domain_dim: int, convVel: list[sp.Expr], diffCoeff: list[sp.Expr], coord_system: CoordinatesSystemType = CoordinatesSystemType.CARTESIAN):
        nb_fields = manuf_sol_container.getNumberFields()
        convTensor = []
        # diffTensor = []
        for d in range(domain_dim):
            tensConv = initListOfLists(nb_fields,nb_fields,0.0)
            # tensDiff = initListOfLists(nb_fields,nb_fields,0.0)
            for m in range(nb_fields):
                for n in range(nb_fields):
                    if m == n: 
                        tensConv[m][n] = convVel[d]
                        # tensDiff[m][n] = diffCoeff[d]
            convTensor.append(tensConv)
            # diffTensor.append(tensDiff)

        diffTensor = []
        for d1 in range(domain_dim):
            for d2 in range(domain_dim):
                tensDiff = initListOfLists(nb_fields,nb_fields,0.0)
                if d1 == d2:
                    for m in range(nb_fields):
                        for n in range(nb_fields):
                            if m == n: 
                                tensDiff[m][n] = diffCoeff[d1]
                diffTensor.append(tensDiff)
        
        super().__init__(tag, manuf_sol_container, domain_dim, convTensor, diffTensor, coord_system)

        # self.advDiffModel = LinearAdvectionDiffusionModels(tag, manuf_sol_container, domain_dim, convTensor, diffTensor, coord_system)

    