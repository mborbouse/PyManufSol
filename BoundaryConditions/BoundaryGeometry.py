
from Common.Utilities import *
import sympy as sp
from math import sqrt, fabs
from Common.MMSTags import ProjectionType
from sympy.core.sympify import kernS

class GeneralBoundaryGeometry:

    def __init__(self, sym_variables: list, tag: str):
        self.sym_variables = sym_variables
        self.domain_dim = len(sym_variables)
        self.tag = tag
        self.boundary_eq = None

    def getSymVariables(self):
        return self.sym_variables
    
    def setSymVariables(self,new_sym_var):
        self.sym_variables = new_sym_var

    def getBoundaryTag(self):
        return self.tag

    def getBoundaryEquation(self):
        return self.boundary_eq

    def setBoundaryEquation(self, bnd_eq):
        self.boundary_eq = bnd_eq

    def setSymVar(self, sym_var):
        self.sym_variables = sym_var

    def buildBoundaryEquation(self):
        raise NotImplementedError()

    def getDomainDim(self):
        return self.domain_dim

    def evalBoundaryAt(self, num_coord: dict):
        return subsNumParams(self.getBoundaryEquation(), num_coord)

    def getNormalToBoundary(self):
        normal = [0.] * self.domain_dim
        for i in range(self.domain_dim):
            normal[i] = sp.diff(self.getBoundaryEquation(), self.sym_variables[i])
        unit_normal = normalizeVec(normal)
        return unit_normal

    def getTangentsToBoundary(self):
        normal = self.getNormalToBoundary()
        nb_row = self.domain_dim-1
        n_col = self.domain_dim
        tangents = initListOfLists(nb_row, n_col)
        if self.domain_dim == 2:
            tangents[0][0] = -normal[1]
            tangents[0][1] = normal[0]
        elif self.domain_dim == 3:
            for i in range(3):
                if (fabs(normal[i]) > 1.0/sqrt(3.0)):
                    ip1 = (i+1)%3
                    ip2 = (i+2)%3
                    tangents[0][i] = normal[ip1]
                    tangents[0][ip1] = -normal[i]
                    tangents[0][ip2] = 0.
                    break
            tangents[0] = normalizeVec(tangents[0])

            tangents[1][1] = -(tangents[0][1]*normal[2] - tangents[0][2]*normal[1])
            tangents[1][2] = -(tangents[0][2]*normal[0] - tangents[0][0]*normal[2])
            tangents[1][3] = -(tangents[0][0]*normal[1] - tangents[0][1]*normal[0])
            tangents[1] = normalizeVec(tangents[1])
        return tangents

    def evalNormalAt(self, num_coord: dict):
        return subsNumParams(self.getNormalToBoundary(), num_coord)

class BoundaryGeometryFromEquation(GeneralBoundaryGeometry):

    def __init__(self, sym_variables: list, boundary_eq_coeffs: list, tag: str):
        super().__init__(sym_variables, tag)
        self.boundary_eq_coeffs = boundary_eq_coeffs
        self.buildBoundaryEquation()

    def buildBoundaryEquation(self):
        """ a*x + b*y + c*z = d """
        boundary_equa = sp.Mul(self.boundary_eq_coeffs[0], self.sym_variables[0], evaluate= False)
        for i in range(1,self.domain_dim):
            boundary_equa =  sp.Add(boundary_equa, sp.Mul(self.boundary_eq_coeffs[i], self.sym_variables[i], evaluate=False), evaluate=False)
        self.boundary_eq = sp.Add(boundary_equa, -1.0*self.boundary_eq_coeffs[-1], evaluate= False)

class TimeDependentBoundaryGeometryFromEquation(BoundaryGeometryFromEquation):

    def __init__(self, sym_variables: list, boundary_eq_coeffs: list, tag: str, temporal_sym_var: list, temporal_boundary_coeffs: list):
        super().__init__(sym_variables, boundary_eq_coeffs, tag)
        self.boundary_eq = sp.Add(self.boundary_eq, sp.Mul(temporal_boundary_coeffs[0], temporal_sym_var[0], evaluate=False), evaluate= False)

class BoundaryGeometryFromPoints(GeneralBoundaryGeometry):

    def __init__(self, sym_variables: list, passage_pts: list[list[float]], tag: str):
        super().__init__(sym_variables, tag)
        self.passage_pts = passage_pts
        nb_pts = len(passage_pts)
        self.order = nb_pts - self.domain_dim + 1
        if self.order > 2 and self.domain_dim == 2:
            raise ValueError(" One-D surface based on coordinates of passage points of order > 2 not implemented yet.")
        if self.order > 1 and self.domain_dim == 3:
            raise ValueError(" Two-D surface based on coordinates of passage points of order > 1 not implemented yet.")
        self.buildBoundaryEquation()

    def buildBoundaryEquation(self):
        """ n_x*(x-x_0) + n_y*(y-y_0) + n_z*(z-z_0) = 0 """
        n = self.computeNormalVectorFromPts()
        eq = n[0] * (self.sym_variables[0]-self.passage_pts[0][0])
        for i in range(1,self.domain_dim):
            eq = eq + n[i] * (self.sym_variables[i]-self.passage_pts[0][i])
        self.boundary_eq = eq

    def computeNormalVectorFromPts(self):
        n = [0.0] * self.domain_dim
        if self.domain_dim == 2:
            # n_x = -(y_1 - y_0), n_y = (x_1 - x_0)
            n[0] = -(self.passage_pts[1][1]-self.passage_pts[0][1])
            n[1] = (self.passage_pts[1][0]-self.passage_pts[0][0])
        elif self.domain_dim == 3:
            # t_x_1 = (x_1 - x_0), t_y_1 = (y_1 - y_0)
            # t_x_2 = (x_2 - x_0), t_y_2 = (y_2 - y_0)
            t1 = [0.0] * self.domain_dim
            t2 = [0.0] * self.domain_dim
            for i in range(self.domain_dim):
                t1[i] = self.passage_pts[1][i]-self.passage_pts[0][i]
                t2[i] = self.passage_pts[2][i]-self.passage_pts[0][i]
            # n = t1 x t2 (cross-product)
            n[0] = t1[1]*t2[2] - t1[2]*t2[1]
            n[1] = -(t1[0]*t2[2] - t1[2]*t2[0])
            n[2] = t1[0]*t2[1] - t1[1]*t2[0]
        return n
    

def projectBoundaryVarAlongDir(var, boundary_geo: GeneralBoundaryGeometry, projectToDir: ProjectionType):
    if not isinstance(var, list) or projectToDir == ProjectionType.NOPROJECTION:
        return var
    else:
        n_boundary = boundary_geo.getNormalToBoundary()
        t_boundary = boundary_geo.getTangentsToBoundary()
        size_var = len(var)
        if projectToDir == ProjectionType.NORMAL:
            # if size_var != len(n_boundary):
            #     raise ValueError("Size of vector to be projected different from size of normal vector.")
            return projectThisVecOrTensorAlongThisDirection(var, n_boundary)
        elif projectToDir == ProjectionType.TANGENT:
            if size_var != len(t_boundary[0]):
                raise ValueError("Size of vector to be projected different from size of tangent vector.")
            if boundary_geo.getDomainDim() == 2:
                return projectThisVecOrTensorAlongThisDirection(var, t_boundary[0])
            elif boundary_geo.getDomainDim() == 3:
                return [projectThisVecOrTensorAlongThisDirection(var, t_boundary[0]), projectThisVecOrTensorAlongThisDirection(var, t_boundary[1])]
        elif projectToDir == ProjectionType.NORMALNORMAL:
            return projectThisTensorAlongThoseDirections(var, n_boundary, n_boundary)
        elif projectToDir == ProjectionType.NORMALTANGENT: 
            if boundary_geo.getDomainDim()== 2:
                return projectThisTensorAlongThoseDirections(var, n_boundary, t_boundary[0])
            elif boundary_geo.getDomainDim() == 3:
                return [projectThisTensorAlongThoseDirections(var, n_boundary, t_boundary[0]), projectThisTensorAlongThoseDirections(var, n_boundary, t_boundary[1])]
        else:
            raise ValueError("Unknown projection type on boundaries.")
