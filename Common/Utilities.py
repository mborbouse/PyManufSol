from typing import Optional
import sympy as sp
from Common.MMSTags import CoordinatesSystemType, OutputFileType
import sympy as sp
import math
import pathlib
import csv
import os
import numpy as np
from matplotlib import pyplot as plt

#* --- OPERATIONS ON VECTORS ------------------------------------------------ *#

def initListOfLists(nb_row, nb_col, init_val: float = 0.0):
    return [[init_val] * nb_col for i in range(nb_row)]

def copyListOfLists(input_list, output_list):
    nb_row = len(input_list)
    if nb_row != len(output_list):
        raise ValueError("Nb of rows of both lists of lists not compatible for copy.")
    for i in range(nb_row):
        nb_col = len(input_list[i])
        if nb_col != len(output_list[i]):
            raise ValueError("Nb of cols of both lists of lists not compatible for copy.")
        for j in range(nb_col):
            output_list[i][j] = input_list[i][j]

def flattenNestedList(input_list: list[list]) -> list:
    return [item for sublist in input_list for item in sublist]

def multiplyNestedListsWithList(tens: list[list], vec: list) -> list:
    nb_row = len(tens)
    nb_col = len(tens[0])
    if not all(len(sub) == nb_col for sub in tens):
        raise ValueError("Lengths of lists composing the tensor must be the same.")
    if nb_col != len(vec):
        raise ValueError("Number of columns of tensor not compatible with size of vector for matrix vector multiplication.")
    result = [0.0] * nb_row
    for i in range(nb_row):
        result[i] = sum([a*b for a,b in zip(tens[i],vec)])
    return result

def computeNorm(vec: list):
    norm = 0.0
    for i in vec:
        norm = norm + sp.Pow(i,2)
    return sp.Pow(norm,0.5)

def normalizeVec(vec: list):
    norm = computeNorm(vec)
    return [i/norm for i in vec]

def subsNumParams(var, num_params: Optional[dict] = None):
    sub_params = isinstance(num_params,dict)
    if sub_params and num_params:
        if isinstance(var, list):
            # return [i.subs(num_params) for i in var]
            return [subsNumParams(i,num_params) for i in var]
        else:
            if isinstance(var, sp.Expr): 
                return var.subs(num_params)
            else:
                return var
    else:
        return var

def subsNumParamsBis(var: sp.Expr, sym_var: list[sp.Expr], num_var: list):
    if len(sym_var) != len(num_var):
        raise ValueError("Sizes of sym variables and numeric values must be the same.")
    dic_val = dict()
    for i,j in zip(sym_var,num_var):
        dic_val[i] = j
    return subsNumParams(var, dic_val)

def is_linear(expr: sp.Expr, vars: list[sp.Expr]):
    for x in vars:
        for y in vars:
            try: 
                if not sp.Eq(sp.diff(expr, x, y), 0):
                    return False
            except TypeError:
                return False
    return True

#* --- PROJECTION OF VECTORS AND TENSORS ALONG CHOSEN DIRECTIONS ------------ *#

def projectThisVecAlongThisDirection(vec: list, dir: list):
    nb_var = len(vec)
    if nb_var != len(dir):
        raise ValueError("Size of vector to be projected different from direction of projection.")
    return sum([vec[i]*dir[i] for i in range(nb_var)])

def projectThisTensorAlongThisDirection(tens: list[list], dir: list) -> list:
    nb_lists = len(tens)
    return [projectThisVecAlongThisDirection(tens[i], dir) for i in range(nb_lists)]

def projectThisVecOrTensorAlongThisDirection(var, dir: list):
    if all(isinstance(el, list) for el in var):
        return projectThisTensorAlongThisDirection(var, dir)
    elif isinstance(var, list):
        return projectThisVecAlongThisDirection(var, dir)
    else:
        raise ValueError("The quantity to be projected must be a list or a list of lists.")

def projectThisTensorAlongThoseDirections(tens: list[list], dir1: list, dir2: list):
    if not all(isinstance(el, list) for el in tens):
        raise ValueError("The variable to be projected along two directions must be a list of lists.")
    vec_proj = projectThisTensorAlongThisDirection(tens, dir1)
    return projectThisVecAlongThisDirection(vec_proj, dir2)


#* --- DERIVATIVES COMPUTATION IN VARIOUS REFERENCE FRAMES ------------------ *#
# From "Viscous Fluid Flow" by Papanastasiou, Georgiou, Alexandrou, Tables 1.4, 1.6

# scal => grad(scal) = [d(scal)/dx d(scal)/dy d(scal)/dz] (in Cartesian coordinates)
def computeGradOfScalar(var: sp.Expr, indep_vars: list, coord: CoordinatesSystemType) -> list[sp.Expr]:
    if not isinstance(var, sp.Expr):
        raise ValueError("The scalar on which the gradient wants to be computed must be a symbolic expression scalar type.")
    dim = len(indep_vars)
    grad = [0.0] * dim
    for i in range(dim):
            grad[i] = sp.diff(var, indep_vars[i])
    if coord == CoordinatesSystemType.CYLINDRICAL:
        r = indep_vars[0]
        grad[1] = (1./r) * grad[1]
    return grad

# vec = [v1, v2, v3] 
# => grad(vec) = [ [dv1/dx dv1/dy dv1/dz],
#                  [dv2/dx dv2/dy dv2/dz],
#                  [dv3/dx dv3/dy dv3/dz] ] (in Cartesian coordinates)
def computeGradOfVector(var: list[sp.Expr], indep_vars: list, coord: CoordinatesSystemType) -> list[list[sp.Expr]]:
    if not all(isinstance(el, sp.Expr) for el in var):
        raise ValueError("The vector on which the gradient wants to be computed must be a list of symbolic expressions type.")
    nb_row = len(var)
    nb_col = len(indep_vars)
    grad = initListOfLists(nb_row, nb_col)
    for i in range(nb_row):
            for j in range(nb_col):
                grad[i][j] = sp.diff(var[i],indep_vars[j])
    if coord == CoordinatesSystemType.CYLINDRICAL:
        r = indep_vars[0]
        grad[0][1] = (1./r) * (grad[0][1] - var[1])
        grad[1][1] = (1./r) * (grad[1][1] + var[0])
        if nb_row == 3:
            grad[2][1] = (1./r) * grad[2][1]
    return grad

# vec = [v1, v2, v3]
# => div(vec) = dv1/dx + dv2/dy + dv3/dz (in Cartesian coordinates)
def computeDivergenceOfVector(var: list[sp.Expr], indep_vars: list, coord: CoordinatesSystemType) -> sp.Expr:
    if not all(isinstance(el, sp.Expr) for el in var):
        raise ValueError("The vector on which the divergence wants to be computed must be a list of symbolic expressions type.")
    nb_row = len(var)
    if nb_row != len(indep_vars):
        raise ValueError("Size of vector must be equal to number of indep vars for divergence computation.")
    div = 0.0
    for i in range(nb_row):
        if coord == CoordinatesSystemType.CYLINDRICAL:
            r = indep_vars[0]
            if i == 0:
                div = div + (1./r) * sp.diff(r * var[i],indep_vars[i])
            elif i == 1:
                div = div + (1./r) * sp.diff(var[i],indep_vars[i])
        elif coord == CoordinatesSystemType.CARTESIAN:
            div = div + sp.diff(var[i],indep_vars[i])
    return div

# tens = = [ [t_xx t_xy t_xz],
#            [t_yx t_yy t_yz],
#            [t_zx t_zy t_zz] ]
# div(tens) = [ (d(t_xx)/dx + d(t_yx)/dy + d(t_zx)/dz) 
#               (d(t_xy)/dx + d(t_yy)/dy + d(t_zy)/dz)
#               (d(t_xz)/dx + d(t_yz)/dy + d(t_zy)/dz) ] 
# (in Cartesian coordinates)
def computeDivergenceOfTensor(var: list[list[sp.Expr]], indep_vars: list, coord: CoordinatesSystemType) -> list[sp.Expr]:
    for l in var:
        if not all(isinstance(el, sp.Expr) for el in l):
            raise ValueError("The tensor on which the divergence wants to be computed must be a nested list of symbolic expressions type.")
    nb_row = len(var)
    if nb_row != len(indep_vars):
        raise ValueError("Number of rows of tensor must be equal to number of indep vars for divergence computation.")
    nb_col = len(var[0])
    if not all(len(sub) == nb_col for sub in var):
        raise ValueError("Lengths of lists composing the tensor must be the same.")
    div = [0.0] * nb_col
    if coord == CoordinatesSystemType.CYLINDRICAL:
        r = indep_vars[0]
        div[0] = 1./r * (sp.diff(r*var[0][0],indep_vars[0]) + sp.diff(var[1][0],indep_vars[1]) - var[1][1])
        div[1] = 1./(r*r) * (sp.diff(r*r*var[0][1],indep_vars[0])) + 1./r * (sp.diff(var[1][1],indep_vars[1]) - (var[1][0]-var[0][1]))
        if len(indep_vars) == 3:
            div[2] =  1./r * (sp.diff(r*var[0][2],indep_vars[0]) + sp.diff(var[1][2],indep_vars[1]))
            for j in range(nb_col):
                div[j] = div[j] + sp.diff(var[2][j],indep_vars[2])
    elif coord == CoordinatesSystemType.CARTESIAN:
        for j in range(nb_col):
            for i in range(nb_row):
                div[j] = div[j] + sp.diff(var[i][j],indep_vars[i])

    return div

# vec = [v1, v2, v3]
# curl(vec) = [dv3/dy-dv2/dz dv1/dz-dv3/dx dv2/dx-dv1/dy] 
# (in Cartesian coordinates)
def computeCurlOfVector(var: list[sp.Expr], indep_vars: list, coord: CoordinatesSystemType) -> list[sp.Expr]:
    if not all(isinstance(el, sp.Expr) for el in var):
        raise ValueError("The vector on which the curl wants to be computed must be a list of symbolic expressions type.")
    size_vec = len(var)
    if size_vec != len(indep_vars):
        raise ValueError("Size of vector must be equal to size of indep variables for curl computation.")
    size_curl = 3 if size_vec == 3 else 1
    curl = [0.0] * size_vec
    if size_curl == 1:
        if coord == CoordinatesSystemType.CYLINDRICAL:
            r = indep_vars[0]
            curl[0] = 1./r * (sp.diff(r * var[1],indep_vars[0]) - sp.diff(var[0],indep_vars[1]))
        elif coord == CoordinatesSystemType.CARTESIAN:
            curl[0] = sp.diff(var[1],indep_vars[0]) - sp.diff(var[0],indep_vars[1])
    else:
        if coord == CoordinatesSystemType.CYLINDRICAL:
            curl[0] = 1./r * (sp.diff(var[2],indep_vars[1])) - sp.diff(var[1],indep_vars[2])
            curl[1] = sp.diff(var[0],indep_vars[2]) - sp.diff(var[2],indep_vars[0])
            curl[2] = 1./r * (sp.diff(r * var[1],indep_vars[0]) - sp.diff(var[0],indep_vars[1]))
        elif coord == CoordinatesSystemType.CARTESIAN:
            list_ind = [i for i in range(size_curl)]
            for i in range(size_curl): 
                ind = [ele for ele in list_ind if ele not in [i]]
                curl[i] = ((-1)**i) * (sp.diff(var[ind[-1]],indep_vars[ind[0]]) - sp.diff(var[ind[0]],indep_vars[ind[-1]]))


#* --- CHANGE OF COORDINATES ------------------------------------------------ *#

def currentSpatialVarToNewCoordSystem(currentCoordSystem: CoordinatesSystemType, newSpatialVar: list[sp.Expr], newCoordSystem: CoordinatesSystemType) -> list[sp.Expr]:
    nb_dim = len(newSpatialVar)
    # if nb_dim != len(newSpatialVar):
    #     raise ValueError("Current and new coordinates system must have the same number of spatial coordinates.")
    new_coord = [0.0] * nb_dim
    if currentCoordSystem == newCoordSystem:
        new_coord = newSpatialVar
    elif currentCoordSystem == CoordinatesSystemType.CARTESIAN:
        if newCoordSystem == CoordinatesSystemType.CYLINDRICAL:
            r = newSpatialVar[0]
            theta = newSpatialVar[1]
            x = r * sp.cos(theta)
            y = r * sp.sin(theta)
            if nb_dim == 2:
                new_coord = [x, y]
            else:
                z = newSpatialVar[3]
                new_coord = [x, y, z]
        else:
            raise ValueError("Change of spatial variables from Cartesian to this other reference system not available.")
    elif currentCoordSystem == CoordinatesSystemType.CYLINDRICAL:
        if newCoordSystem == CoordinatesSystemType.CARTESIAN:
            x = newSpatialVar[0]
            y =  newSpatialVar[1]
            r = sp.sqrt(x*x + y*y)
            theta = sp.atan2(y, x)
            if nb_dim == 2:
                new_coord = [r, theta]
            else:
                z = newSpatialVar[3]
                new_coord = [r, theta, z]
        else:
            raise ValueError("Change of spatial variables from Cylindrical to this other reference system not available.")
    else:
        raise ValueError("Unknown reference coordinates system type.")     
    return new_coord

def translateCurrentSpatialVar(newSpatialVar: list[sp.Expr], translation: list[float]) -> list[sp.Expr]:
    dim = len(newSpatialVar)
    # new_coord = newSpatialVar
    new_coord = [0] * dim
    if isinstance(translation, list):
        if len(translation) != dim:
            raise ValueError("Size of translation vector must correspond to spatial dimension.")
        for i in range(dim):
            new_coord[i] = newSpatialVar[i] - translation[i]
    return new_coord 

def transfoMatrixRotationAroundOneAxis(dim: int, rotation: float, ind_axis_rotation: int = 2) -> list[sp.Expr]:
    mat = initListOfLists(dim, dim, 0.0)
    theta = math.radians(rotation)
    ind_rot = ind_axis_rotation if dim == 3 else 2
    ind_coord = [ele for ele in [0, 1, 2] if ele not in [ind_rot]]
    sign_sin = (-1)**ind_rot
    mat[ind_coord[0]][ind_coord[0]] = mat[ind_coord[1]][ind_coord[1]] = sp.cos(theta)
    if dim == 3:
        mat[ind_rot][ind_rot] = 1.0
    mat[ind_coord[0]][ind_coord[1]] = -sign_sin * sp.sin(theta)
    mat[ind_coord[1]][ind_coord[0]] = sign_sin * sp.sin(theta)
    return mat

def rotateCurrentSpatialVar(newSpatialVar: list[sp.Expr], rotation: list[list[int,float]]) -> list[sp.Expr]:
    dim = len(newSpatialVar)
    if dim == 2 and len(rotation) != 1:
        raise ValueError("In 2-D, only a single rotation axis must be provided.")
    # new_coord = newSpatialVar
    # theta = math.radians(rotation)
    # ind_rot = ind_axis_rotation if dim == 3 else 2
    # ind_coord = [ele for ele in [0, 1, 2] if ele not in [ind_rot]]
    # sign_sin = (-1)**ind_coord
    # new_coord[ind_coord[0]] = new_coord[ind_coord[0]]*sp.cos(theta) - sign_sin*new_coord[ind_coord[1]]*sp.sin(theta)
    # new_coord[ind_coord[1]] = sign_sin*new_coord[ind_coord[0]]*sp.sin(theta) + new_coord[ind_coord[1]]*sp.cos(theta)
    new_coord = newSpatialVar
    for i in range(len(rotation)):
        transfo_mat = transfoMatrixRotationAroundOneAxis(dim, rotation[i][1], rotation[i][0])
        new_coord = multiplyNestedListsWithList(transfo_mat, new_coord)
    return new_coord

def transfoMatrixFromCurrentToNewCoordSystem(currentCoordSystem: CoordinatesSystemType, newSpatialVar: list[sp.Expr], newCoordSystem: CoordinatesSystemType) -> list[list[sp.Expr]]:
    nb_dim = len(newSpatialVar)
    mat = initListOfLists(nb_dim, nb_dim, 0.0)
    if currentCoordSystem == newCoordSystem:
        for i in range(nb_dim):
            mat[i][i] = 1.0
    elif currentCoordSystem == CoordinatesSystemType.CARTESIAN:
        if newCoordSystem == CoordinatesSystemType.CYLINDRICAL:
            theta = newSpatialVar[1]
            mat[0][0] = sp.cos(theta)
            mat[0][1] = sp.sin(theta)
            mat[1][0] = -sp.sin(theta)
            mat[1][1] = sp.cos(theta)
            if nb_dim == 3:
                mat[2][2] = 1.0
        else:
            raise ValueError("Transfo matrix from Cartesian to this other reference system not available.")
    elif currentCoordSystem == CoordinatesSystemType.CYLINDRICAL:
        if newCoordSystem == CoordinatesSystemType.CARTESIAN:
            theta = sp.atan2(newSpatialVar[1], newSpatialVar[0])
            mat[0][0] = sp.cos(theta)
            mat[0][1] = -sp.sin(theta)
            mat[1][0] = sp.sin(theta)
            mat[1][1] = sp.cos(theta)
            if nb_dim == 3:
                mat[2][2] = 1.0
        else:
            raise ValueError("Transfo matrix from Cylindrical to this other reference system not available.")
    else:
        raise ValueError("Unknown reference coordinates system type for transformation matrix computation.")
    return mat

def expressScalarInNewCoordSystem(currentScalar: sp.Expr, currentSpatialVar: list[sp.Expr], currentCoordSystem: CoordinatesSystemType, newSpatialVar: list[sp.Expr], newCoordSystem: CoordinatesSystemType) -> sp.Expr:
    nb_dim = len(currentSpatialVar)
    newCoordSystemVar = currentSpatialVarToNewCoordSystem(currentCoordSystem, newSpatialVar, newCoordSystem)
    newScalar = currentScalar
    for i in range(nb_dim):
        newScalar = newScalar.subs({currentSpatialVar[i]: newCoordSystemVar[i]})
    # if nb_dim == 2:
    #     newScalar = currentScalar.subs({currentSpatialVar[0]: newCoordSystemVar[0], currentSpatialVar[1]: newCoordSystemVar[1]})
    # else:
    #     newScalar = currentScalar.subs({currentSpatialVar[0]: newCoordSystemVar[0], currentSpatialVar[1]: newCoordSystemVar[1], currentSpatialVar[2]: newCoordSystemVar[2]})
    return newScalar

def translateScalarInNewCoordSystem(currentScalar: sp.Expr, currentSpatialVar: list[sp.Expr], newSpatialVar: list[sp.Expr], translation: list[float]) -> sp.Expr:
    nb_dim = len(currentSpatialVar)
    newCoordSystemVar = translateCurrentSpatialVar(newSpatialVar, translation)
    newScalar = currentScalar
    for i in range(nb_dim):
        newScalar = newScalar.subs({currentSpatialVar[i]: newCoordSystemVar[i]})
    return newScalar

def rotateScalarInNewCoordSystem(currentScalar: sp.Expr, currentSpatialVar: list[sp.Expr], newSpatialVar: list[sp.Expr], rotation: list[list[int,float]]) -> sp.Expr:
    nb_dim = len(currentSpatialVar)
    newCoordSystemVar = rotateCurrentSpatialVar(newSpatialVar, rotation)
    newScalar = currentScalar
    for i in range(nb_dim):
        newScalar = newScalar.subs({currentSpatialVar[i]: newCoordSystemVar[i]})
    return newScalar

def expressVectorInNewCoordSystem(currentVector: list[sp.Expr], currentSpatialVar: list[sp.Expr], currentCoordSystem: CoordinatesSystemType, newSpatialVar: list[sp.Expr], newCoordSystem: CoordinatesSystemType) -> list[sp.Expr]:
    transfo_mat = transfoMatrixFromCurrentToNewCoordSystem(currentCoordSystem, newSpatialVar, newCoordSystem)
    currentVecInNewCoord = [expressScalarInNewCoordSystem(v, currentSpatialVar, currentCoordSystem, newSpatialVar, newCoordSystem) for v in currentVector]
    newVectorInNewCoord = multiplyNestedListsWithList(transfo_mat, currentVecInNewCoord)
    return newVectorInNewCoord

def translateVectorInNewCoordSystem(currentVector: list[sp.Expr], currentSpatialVar: list[sp.Expr], newSpatialVar: list[sp.Expr], translation: list[float]) -> list[sp.Expr]:
    return [translateScalarInNewCoordSystem(v, currentSpatialVar, newSpatialVar, translation) for v in currentVector]

def rotateVectorInNewCoordSystem(currentVector: list[sp.Expr], currentSpatialVar: list[sp.Expr], newSpatialVar: list[sp.Expr], rotation: list[list[int,float]]) -> list[sp.Expr]:
    currentVecInNewCoord = [rotateScalarInNewCoordSystem(v, currentSpatialVar, newSpatialVar, rotation) for v in currentVector]
    new_vec = currentVecInNewCoord
    for i in range(len(rotation)):
        transfo_mat = transfoMatrixRotationAroundOneAxis(len(newSpatialVar), -rotation[i][1], rotation[i][0])
        new_vec = multiplyNestedListsWithList(transfo_mat, new_vec)
    return new_vec


def shiftSolution(list_sol, symb_var, shift_values):
    if len(shift_values) == 0:
        return list_sol
    list_sol_shifted = []
    for i in list_sol:
        sol_shifted = i
        for j in range(len(symb_var)):
            sol_shifted = sol_shifted.subs({symb_var[j]: symb_var[j]-shift_values[j]})
        list_sol_shifted.append(sol_shifted)
    return list_sol_shifted


#* --- PRINTING IN FILE ----------------------------------------------------- *#
def defaultFilePath() -> str:
    return str(pathlib.Path().resolve())+"/"
    # return str(pathlib.Path(__file__).parent.resolve())+"/"

def openFile(fileName, fileType: OutputFileType = OutputFileType.TEXT) -> list:
    os.makedirs(os.path.dirname(fileName), exist_ok=True)
    f = open(fileName,'w')
    if fileType == OutputFileType.TEXT:
        return [f]
    elif fileType == OutputFileType.CSV:
        return [f, csv.writer(f, delimiter=',')]

def closeFile(file: list) -> None:
    file[0].close()

def writeOneLineInFile(file: list, val_str: str, fileType: OutputFileType = OutputFileType.TEXT) -> None:
    if fileType == OutputFileType.TEXT:
        file[0].write(val_str)
        file[0].write('\n')
    elif fileType == OutputFileType.CSV:
        file[1].writerow([val_str])

def extensionName(fileType: OutputFileType = OutputFileType.TEXT):
    if fileType == OutputFileType.TEXT:
        return ".txt"
    elif fileType == OutputFileType.CSV:
        return ".csv"


#* --- COMPUTATIONS ON GEOMETRICAL ENTITIES --------------------------------- *#

def findIntersectionBetweenTwo2DLines(curve_1: sp.Expr, curve_2: sp.Expr, sym_variables: list[sp.Expr]) -> list[sp.Expr]:
    sol = sp.solve([curve_1, curve_2], tuple(sym_variables))
    result = []
    for s in sym_variables:
        if s in sol:
            result.append(sol[s])
        # else:
        #     result.append(sp.Mul(0.0, s))
    return result
    # ind_solve = 1
    # y_1 = sp.solve(curve_1, sym_variables[ind_solve])
    # y_2 = sp.solve(curve_2, sym_variables[ind_solve])
    # if len(y_1) > 0 and all(var not in y_1[0].free_symbols for var in sym_variables) and len(y_2) == 0:
    #     y_2 = sp.solve(curve_2, sym_variables[1-ind_solve])
    #     if len(y_2) > 0 and all(var not in y_2[0].free_symbols for var in sym_variables):
    #         return [y_2[0], y_1[0]]
    # elif len(y_2) > 0 and all(var not in y_2[0].free_symbols for var in sym_variables) and len(y_1) == 0:
    #     y_1 = sp.solve(curve_1, sym_variables[1-ind_solve])
    #     if len(y_1) > 0 and all(var not in y_1[0].free_symbols for var in sym_variables):
    #         return [y_1[0], y_2[0]]
    # y_1 = [0.0] if len(y_1) == 0 else y_1
    # y_2 = [0.0] if len(y_2) == 0 else y_2
    # for i in sym_variables:
    #     y_1[0] = sp.Add(y_1[0], sp.Mul(0.0,i,evaluate=False), evaluate=False)
    #     y_2[0] = sp.Add(y_2[0], sp.Mul(0.0,i,evaluate=False), evaluate=False)
    # x = sp.solve(sp.Add(y_1[0],-1.0*y_2[0],evaluate=False), sym_variables[1-ind_solve])
    # if len(x) > 0:
    #     y = y_1[0].subs({sym_variables[ind_solve]: x})
    #     return [x[0], y]
    # else:
    #     return []

def sampleLineBtwTwoPtsOnCurve(curve: sp.Expr, pt_1: list[float], pt_2: list[float], nb_pts: int, sym_variables: list[sp.Expr], sampleDir: int) -> list[list[float]]:
    ind = sampleDir
    sample_pts = []
    pt_1_save = pt_1
    if pt_1[ind] > pt_2[ind]:
        pt_1 = pt_2
        pt_2 = pt_1_save
    x_sample = np.linspace(float(pt_1[ind]), float(pt_2[ind]), nb_pts)
    y_sol = sp.solve(curve, sym_variables[1-ind])
    for i in range(nb_pts):
        if ind == 0:
            sample_pts.append([x_sample[i], y_sol[0].subs({sym_variables[ind]: x_sample[i]})])
        elif ind == 1:
            sample_pts.append([y_sol[0].subs({sym_variables[ind]: x_sample[i]}), x_sample[i]])
    return sample_pts, x_sample

def plotVarBtwTwoPts(var: sp.Expr, curve: sp.Expr, pt_1: list[float], pt_2: list[float], nb_pts: int, sym_variables: list[sp.Expr], sampleDir: int, name_var: str, line_color: Optional[str] = None, line_type: Optional[str] = None):
    sample_pts, x_sample = sampleLineBtwTwoPtsOnCurve(curve, pt_1, pt_2, nb_pts, sym_variables, sampleDir)
    var_eval = np.zeros(nb_pts)
    for i in range(nb_pts):
        val = subsNumParamsBis(var, [sym_variables[0], sym_variables[1]], [sample_pts[i][0], sample_pts[i][1]])
        if isinstance(val, list):
            val = val[0]
        var_eval[i] = val
    plt.plot(x_sample, var_eval, line_type, linewidth=3, color=line_color, label="%s"%(name_var))

# -------
# def cartesianToCylindricalScalar(symb_scalar: sp.Expr, symb_var, cart_symb_var):
#     x = symb_var[0]*sp.cos(symb_var[1])
#     y = symb_var[0]*sp.sin(symb_var[1])
#     var_ra = symb_scalar.subs({cart_symb_var[0]: x, cart_symb_var[1]: y})
#     return var_ra

# def cylindricalToCartesianScalar(symb_scalar, symb_var, cylind_symb_var):
#     theta = sp.atan2(symb_var[1], symb_var[0])
#     r = sp.sqrt(symb_var[0]*symb_var[0]+symb_var[1]*symb_var[1])
#     var_xy = symb_scalar.subs({cylind_symb_var[0]: r, cylind_symb_var[1]: theta})
#     return var_xy

# def rotateAndTranslateCartesianScalar(symb_scalar, symb_var, new_cartesian_symb_var, operations):
#     theta = sp.atan2(symb_var[1], symb_var[0])
#     r = sp.sqrt(symb_var[0]*symb_var[0]+symb_var[1]*symb_var[1])
#     theta = operations[0]
#     shift_x = operations[1]
#     shift_y = operations[2]
#     x_new = (symb_var[0]-shift_x)*sp.cos(theta) - (symb_var[1]-shift_y)*sp.sin(theta)
#     y_new = (symb_var[0]-shift_x)*sp.sin(theta) + (symb_var[1]-shift_y)*sp.cos(theta)
#     var_xy = symb_scalar.subs({new_cartesian_symb_var[0]: x_new, new_cartesian_symb_var[1]: y_new})
#     return var_xy

# def cartesianToCylindricalVector(symb_vec, symb_var, cart_symb_var):
#     v_x_ra = cartesianToCylindricalScalar(symb_vec[0], symb_var, cart_symb_var)
#     v_y_ra = cartesianToCylindricalScalar(symb_vec[1], symb_var, cart_symb_var)
#     v_r = v_x_ra*sp.cos(symb_var[1]) + v_y_ra*sp.sin(symb_var[1])
#     v_theta = -v_x_ra*sp.sin(symb_var[1]) + v_y_ra*sp.cos(symb_var[1])
#     return [v_r, v_theta]

# def cylindricalToCartesianVector(symb_vec, symb_var, cylind_symb_var):
#     theta = sp.atan2(symb_var[1], symb_var[0])
#     v_r_xy = cylindricalToCartesianScalar(symb_vec[0], symb_var, cylind_symb_var)
#     v_theta_xy = cylindricalToCartesianScalar(symb_vec[1], symb_var, cylind_symb_var)
#     v_x = v_r_xy*sp.cos(theta) - v_theta_xy*sp.sin(theta)
#     v_y = v_r_xy*sp.sin(theta) + v_theta_xy*sp.cos(theta)
#     return [v_x, v_y]

# def rotateAndTranslateCartesianVector(symb_vec, symb_var, new_cartesian_symb_var, operations):
#     v_x_new = rotateAndTranslateCartesianScalar(symb_vec[0], symb_var, new_cartesian_symb_var, operations)
#     v_y_new = rotateAndTranslateCartesianScalar(symb_vec[1], symb_var, new_cartesian_symb_var, operations)
#     return [v_x_new, v_y_new]