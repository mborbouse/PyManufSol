import numpy as np
from scipy.optimize import fsolve
from Common.MMSTags import FluidSolutionTags as fst

class MieGruneisenEOS:
    def __init__(self, gamma: float, Cv: float, dim: int):
        self.gamma = gamma
        self.Cv = Cv
        self.dim = dim
        self.nb_fields = dim+2

    def primitivePVTToConservative(self, U_prim_dic: dict):
        p = U_prim_dic[fst.PRESSURE]
        u = U_prim_dic[fst.VELOCITY_X]
        v = U_prim_dic[fst.VELOCITY_Y]
        w = 0.0
        if self.dim == 3:
            w = U_prim_dic[fst.VELOCITY_Z]
        T = U_prim_dic[fst.TEMPERATURE]
        rho = self.rho_pT(p, T)
        e = self.e_rhoT(rho, T)
        E = e + 0.5 * (u * u + v * v)
        U_cons = [0] * self.nb_fields
        U_cons[0] = rho
        U_cons[1] = rho * u
        U_cons[2] = rho * v
        if self.dim == 3:
            U_cons[3] = rho * w
        U_cons[self.nb_fields-1] = rho * E
        return U_cons

    def primitivePVRHOToConservative(self, U_prim_dic: dict):
        p = U_prim_dic[fst.PRESSURE]
        u = U_prim_dic[fst.VELOCITY_X]
        v = U_prim_dic[fst.VELOCITY_Y]
        w = 0.0
        if self.dim == 3:
            w = U_prim_dic[fst.VELOCITY_Z]
        rho = U_prim_dic[fst.DENSITY]
        T = self.T_rhoP(rho, p)
        e = self.e_rhoT(rho, T)
        E = e + 0.5 * (u * u + v * v)
        U_cons = [0] * self.nb_fields
        U_cons[0] = rho
        U_cons[1] = rho * u
        U_cons[2] = rho * v
        if self.dim == 3:
            U_cons[3] = rho * w
        U_cons[self.nb_fields-1] = rho * E
        return U_cons

    def conservativeToPrimitivePVT(self, U_cons_dic: dict):
        rho = U_cons_dic[fst.DENSITY]
        rhoU = U_cons_dic[fst.MOMENTUM_X]
        rhoV = U_cons_dic[fst.MOMENTUM_Y]
        rhoW = 0.0
        if self.dim == 3:
           rhoW = U_cons_dic[fst.MOMENTUM_Z] 
        rhoE = U_cons_dic[fst.RHOTOTALENERGY]
        u = rhoU / rho
        v = rhoV / rho
        w = 0.0
        if self.dim == 3:
            w = rhoW / rho
        e_kin = 0.5 * (u*u + v*v + w*w)
        e = rhoE/rho - e_kin
        T = self.T_rhoe(rho,e)
        p = self.p_rhoT(rho,T)
        U_prim = [0] * self.nb_fields
        U_prim[0] = p
        U_prim[1] = u
        U_prim[2] = v
        if self.dim == 3:
            U_prim[3] = w
        U_prim[self.nb_fields-1] = T
        return U_prim

    def conservativeToPrimitivePVRHO(self, U_cons_dic: dict):
        rho = U_cons_dic[fst.DENSITY]
        rhoU = U_cons_dic[fst.MOMENTUM_X]
        rhoV = U_cons_dic[fst.MOMENTUM_Y]
        rhoW = 0.0
        if self.dim == 3:
           rhoW = U_cons_dic[fst.MOMENTUM_Z] 
        rhoE = U_cons_dic[fst.RHOTOTALENERGY]
        u = rhoU / rho
        v = rhoV / rho
        w = 0.0
        if self.dim == 3:
            w = rhoW / rho
        e_kin = 0.5 * (u*u + v*v + w*w)
        e = rhoE/rho - e_kin
        T = self.T_rhoe(rho,e)
        p = self.p_rhoT(rho,T)
        rho = self.rho_pT(p, T)
        U_prim = [0] * self.nb_fields
        U_prim[0] = p
        U_prim[1] = u
        U_prim[2] = v
        if self.dim == 3:
            U_prim[3] = w
        U_prim[self.nb_fields-1] = rho
        return U_prim

    def primitivePVRHOToprimitivePVT(self, U_prim_dic: dict):
        p = U_prim_dic[fst.PRESSURE]
        u = U_prim_dic[fst.VELOCITY_X]
        v = U_prim_dic[fst.VELOCITY_Y]
        w = 0.0
        if self.dim == 3:
            w = U_prim_dic[fst.VELOCITY_Z]
        rho = U_prim_dic[fst.DENSITY]
        T = self.T_rhoP(rho, p)
        U_prim = [0] * self.nb_fields
        U_prim[0] = p
        U_prim[1] = u
        U_prim[2] = v
        if self.dim == 3:
            U_prim[3] = w
        U_prim[self.nb_fields-1] = T
        return U_prim

    def primitivePVTToprimitivePVRHO(self, U_prim_dic: dict):
        p = U_prim_dic[fst.PRESSURE]
        u = U_prim_dic[fst.VELOCITY_X]
        v = U_prim_dic[fst.VELOCITY_Y]
        w = 0.0
        if self.dim == 3:
            w = U_prim_dic[fst.VELOCITY_Z]
        T = U_prim_dic[fst.TEMPERATURE]
        rho = self.rho_pT(p, T)
        U_prim = [0] * self.nb_fields
        U_prim[0] = p
        U_prim[1] = u
        U_prim[2] = v
        if self.dim == 3:
            U_prim[3] = w
        U_prim[self.nb_fields-1] = rho
        return U_prim

    def getCv(self,rho):
        return self.Cv

    def getEOSName(self):
        return "MieGruneisen EOS"

    def get_e_ref(self,rho):
        raise NotImplementedError()

    def get_p_ref(self,rho):
        raise NotImplementedError()

    def get_e_th(self,rho,T):
        raise NotImplementedError()

    def p_rhoT(self,rho,T):
        return self.get_p_ref(rho) + rho*(self.gamma-1)*self.get_e_th(rho,T)

    def e_rhoT(self,rho,T):
        return self.get_e_ref(rho) + self.get_e_th(rho,T)

    def nonLinearEqTemperature(self,T,rho):
        fun = self.e_rhoT(rho,T) - self.get_e_ref(rho) - self.get_e_th(rho,T)
        return fun

    def nonLinearEqDensity(self,rho,p,T):
        fun = p - self.p_rhoT(rho,T)
        return fun

    def nonLinearEqTemperatureEnergy(self,T,rho, e):
        fun = e - self.e_rhoT(rho,T)
        return fun

    def T_rhoP(self,rho,P):
        T_guess = np.array([1.0]) # TO BE IMPROVED?
        T = fsolve(self.nonLinearEqTemperature, T_guess, args = (rho))
        return T

    def rho_pT(self,p,T):
        rho_guess = np.array([1.0]) # TO BE IMPROVED?
        rho = fsolve(self.nonLinearEqDensity, rho_guess, args = (p,T))
        return rho

    def T_rhoe(self, rho, e):
        T_guess = np.array([1.0]) # TO BE IMPROVED?
        T = fsolve(self.nonLinearEqTemperatureEnergy, T_guess, args = (rho,e))
        return T

    def nonLinearEqTemperatureSol(self,T, U_sol):
        rho = U_sol[0]
        Ekin = (1/rho)*0.5*(U_sol[1]*U_sol[1]+U_sol[2]*U_sol[2])
        rhoe = U_sol[3]-Ekin
        fun = rho*self.e_rhoT(rho,T) - rhoe
        return fun

    def T_sol(self, U):
        T_guess = np.array([1.0]) # TO BE IMPROVED?
        T = fsolve(self.nonLinearEqTemperatureSol, T_guess, args = (U))
        return T

    def p_sol(self, U):
        rho = U[0]
        T = self.T_sol(U)
        return self.p_rhoT(rho, T)

    def delta_e(self, rho, p, U):
        rho_ref = U[0]
        p_ref = U[2]
        T_ref = self.T_rhoP(rho_ref,p_ref)
        T = self.T_rhoP(rho,p)
        return self.e_rhoT(rho,T) - self.e_rhoT(rho_ref,T_ref)

    def computeDpDrho_e(self, U):
        raise NotImplementedError()  # TO BE GENERALIZED

    def computeDpDe_rho(self, U):
        raise NotImplementedError()  # TO BE GENERALIZED

    def computeDTDrho_e(self, U):
        raise NotImplementedError()  # TO BE GENERALIZED

    def computeDTDe_rho(self, U):
        raise NotImplementedError()  # TO BE GENERALIZED

    def computeSoundSpeed(self, U):
        DpDrho = self.computeDpDrho_e(U)
        DpDe = self.computeDpDe_rho(U)
        p = self.p_sol(U)
        rho = U[0]
        c_2 = DpDrho + (p / (rho * rho)) * DpDe
        return np.sqrt(c_2)

    def computeMachNumber(self, U, n):
        c = self.computeSoundSpeed(U)
        rho = U[0]
        u = U[1] / rho
        v = U[2] / rho
        u_n = u*n[0] + v*n[1]
        return u_n/c

    def computeDpDU(self, U):
        rho = U[0]
        u = U[1]/rho
        v = U[2]/rho
        Ekin = 0.5*(u*u + v*v)
        e = U[3]/rho - Ekin

        dpdrho = self.computeDpDrho_e(U)
        dpde = self.computeDpDe_rho(U)

        dedrho = (Ekin - e) / rho
        dedu = -u / rho
        dedv = -v / rho
        dedrhoE = 1 / rho

        dpdQ = np.zeros(np.size(U))
        dpdQ[0] = dpdrho + dpde * dedrho
        dpdQ[1] = dpde * dedu
        dpdQ[2] = dpde * dedv
        dpdQ[3] = dpde * dedrhoE

        return dpdQ

    def computeDTDU(self, U):
        rho = U[0]
        u = U[1]/rho
        v = U[2]/rho
        Ekin = 0.5*(u*u + v*v)
        e = U[3]/rho - Ekin

        dTdrho = self.computeDTDrho_e(U)
        dTde = self.computeDTDe_rho(U)

        dedrho = (Ekin - e) / rho
        dedu = -u / rho
        dedv = -v / rho
        dedrhoE = 1 / rho

        dTdQ = np.zeros(np.size(U))
        dTdQ[0] = dTdrho + dTde * dedrho
        dTdQ[1] = dTde * dedu
        dTdQ[2] = dTde * dedv
        dTdQ[3] = dTde * dedrhoE

        return dTdQ

    def computeDTDrho(self,rho, T):
        raise NotImplementedError() # TO BE GENERALIZED

    def computeDTDP(self, rho):
        raise NotImplementedError() # TO BE GENERALIZED


# class TaitEOS(MieGruneisenEOS):
#     def __init__(self, rho, T, gamma, Cv, p_0, k_0):
#         MieGruneisenEOS.__init__(self,rho,T,gamma,Cv)
#         self.p_0 = p_0
#         self.k_0 = k_0
#
#     def get_e_ref(self):
#         return (self.k_0-self.p_0)/self.rho
#
#     def get_p_ref(self):
#         return self.p_0-self.k_0
#
#     def get_e_th(self):
#         return self.Cv*self.T

class StiffenedGasEOS(MieGruneisenEOS):
    def __init__(self, gamma, Cv, p_inf, dim):
        MieGruneisenEOS.__init__(self, gamma, Cv, dim)
        self.p_inf = p_inf

    def getEOSName(self):
        return "Stiffened Gas EOS"

    def get_e_ref(self,rho):
        return (self.p_inf) / rho

    def get_p_ref(self,rho):
        return  -self.p_inf

    def get_e_th(self,rho,T):
        return self.Cv * T

    def computeDTDrho(self, rho, T):
        R = (self.gamma-1)*self.Cv
        p = self.p_rhoT(rho,T)
        dTdRho = -(p+self.p_inf)/(rho*rho*R)
        return dTdRho

    def computeDTDP(self, rho):
        R = (self.gamma - 1) * self.Cv
        dTdP = 1/(rho*R)
        return dTdP

    def T_rhoP(self, rho, P):
        R = (self.gamma - 1) * self.Cv
        return (P + self.p_inf)/ (rho * R)

    def rho_pT(self, P, T):
        R = (self.gamma-1)*self.Cv
        return (P + self.p_inf)/ (R * T)

    def T_rhoe(self, rho, e):
        T = (1/self.Cv) * (e - (self.p_inf/rho))
        return T

    def computeDpDrho_e(self, U):
        rho = U[0]
        T = self.T_sol(U)
        return (self.gamma-1)*self.e_rhoT(rho, T)

    def computeDpDe_rho(self, U):
        rho = U[0]
        return rho*(self.gamma - 1)

    def computeDTDrho_e(self, U):
        rho = U[0]
        p_inf = self.p_inf
        Cv = self.Cv
        return (p_inf/Cv)*(1/(rho*rho))

    def computeDTDe_rho(self, U):
        Cv = self.Cv
        return (1.0/Cv)

    def getGammaFromPTRho(self, p, T, rho):
        gamma = ((p + self.p_inf) / (rho*self.Cv*T)) + 1
        return gamma


class PerfectGasEOS(MieGruneisenEOS):
    def __init__(self, gamma, Cv, dim):
        MieGruneisenEOS.__init__(self, gamma, Cv, dim)

    def getEOSName(self):
        return "Perfect Gas EOS"

    def get_e_ref(self,rho):
        return 0.0

    def get_p_ref(self,rho):
        return  0.0

    def get_e_th(self,rho,T):
        return self.Cv * T

    def computeDTDrho(self, rho, T):
        R = (self.gamma-1)*self.Cv
        p = self.p_rhoT(rho,T)
        dTdRho = -(p)/(rho*rho*R)
        return dTdRho

    def computeDTDP(self, rho):
        R = (self.gamma - 1) * self.Cv
        dTdP = 1/(rho*R)
        return dTdP

    # def delta_e(self,rho,p,U):
    #     rho_ref = U[0]
    #     p_ref = U[2]
    #     T_ref = self.T_rhoP(rho_ref,p_ref)
    #     T = self.T_rhoP(rho,p)
    #     return self.Cv * (T - T_ref)

    def T_rhoP(self,rho,P):
        R = (self.gamma-1)*self.Cv
        return P/(rho*R)

    def rho_pT(self, P, T):
        R = (self.gamma-1)*self.Cv
        return P/(R*T)

    def T_rhoe(self, rho, e):
        T = (1/self.Cv) * e
        return T