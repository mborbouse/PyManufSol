import os
from PhysicalModels.ThermoPhysicalModels.EOS_factory import StiffenedGasEOS
from PhysicalModels.ThermoPhysicalModels.TransportProperties_factory import ConstantFluidTransportProperty
from PhysicalModels.LinearAdvectionDiffusionModels import LinearScalarAdvectionDiffusionModels
from ManufSolution.ManufSolShape.PolarPolynomialManufSol import PolarPolynomialManufSol
from ManufSolution.ManufSolShape.UserDefinedManufSol import UserDefinedManufSol
from ManufSolution.ManufSolContainer.GeneralManufSolContainer import GeneralManufSolContainerSinglePhase
from ManufSolution.ManufSolContainer.GeneralManufSolContainer import GeneralManufSolContainerMultiplePhases
from Common.MMSTags import DefaultSolutionTags, ProjectionType
from BoundaryConditions.BoundaryGeometry import TimeDependentBoundaryGeometryFromEquation
from BoundaryConditions.GeneralBoundaryConditions import GeneralBoundaryConditions
from Common.MMSTags import CoordinatesSystemType, OutputFileType, DefaultVarSetTags
from ProblemSolving.GeneralProblemHandler import GeneralProblemHandler, QuantityInfoForPlot
from Outputs.PlotOverArea import *
import sympy as sp

current_directory = os.getcwd()
output_folder_name = "UnsteadyScalarAdvectionDiffusionOutput/"
output_absolute_path = os.path.join(current_directory,output_folder_name)

#* -------------------------------------------------------------------------- *#
#* --- 0/ DEFINITION OF PROBLEM PARAMETERS ---
#* -------------------------------------------------------------------------- *#

domain_dim = 2

u_gamma_0 = 1.0 # velocity of the interface
x_init_gamma = 0.0 # initial position of the interface

jump_sol = -1.0
jump_grad = -1.0

left_bc_sol_val = 0.0
right_bc_grad_val = 1.0

left_bc_pos = -1.0 
right_bc_pos = 1.0

freq = sp.pi/4.0

#* -------------------------------------------------------------------------- *#
#* --- 1/ MANUFACTURED SOLUTION ASSEMBLY ---
#* -------------------------------------------------------------------------- *#

sym_variables = sp.symbols('x, y, t')

u_0_coeff = list(sp.symbols('a_0, b_0'))
u_0_user = u_0_coeff[0] + u_0_coeff[1] * sp.sin(sym_variables[0]*freq)
u_0 = UserDefinedManufSol(u_0_coeff, sym_variables, u_0_user)

u_1_coeff = list(sp.symbols('a_1, b_1'))
u_1_user = u_1_coeff[0] + u_1_coeff[1] * sp.sin(sym_variables[0]*freq)
u_1 = UserDefinedManufSol(u_1_coeff, sym_variables, u_1_user)

# Assembly into dict
manuf_sol_list_0 = [[u_0]]
manuf_sol_list_1 = [[u_1]]
manuf_sol_cont_0 = GeneralManufSolContainerSinglePhase(manuf_sol_list_0,domain_dim, sym_variables[0:2], DefaultVarSetTags.MYVARSETTAG, [sym_variables[-1]])
manuf_sol_cont_1 = GeneralManufSolContainerSinglePhase(manuf_sol_list_1,domain_dim, sym_variables[0:2], DefaultVarSetTags.MYVARSETTAG, [sym_variables[-1]])
dict_symb_sols = GeneralManufSolContainerMultiplePhases([manuf_sol_cont_0, manuf_sol_cont_1], [0, 1])

print(dict_symb_sols.getDicManufSolPerPhaseAndTag())

#* -------------------------------------------------------------------------- *#
#* --- 2/ DEFINITION OF THE GOVERNING EQUATIONS ---
#* -------------------------------------------------------------------------- *#
conv_coeffs_0 = [1.0, 1.0]
diff_coeffs_0 = [1.0, 1.0]

conv_coeffs_1 = [1.0, 1.0]
diff_coeffs_1 = [1.0, 1.0]

name_model_0 = "phase0"
eqs_0 = LinearScalarAdvectionDiffusionModels(name_model_0, manuf_sol_cont_0, domain_dim, conv_coeffs_0, diff_coeffs_0, CoordinatesSystemType.CARTESIAN)

name_model_1 = "phase1"
eqs_1 = LinearScalarAdvectionDiffusionModels(name_model_1, manuf_sol_cont_1, domain_dim, conv_coeffs_1, diff_coeffs_1, CoordinatesSystemType.CARTESIAN)

#* -------------------------------------------------------------------------- *#
#* --- 3/ DEFINITION OF THE BOUNDARY AND JUMPS CONDITIONS (OPTIONAL!) ---
#* -------------------------------------------------------------------------- *#

# Definition of boundaries geometry
vel_interface = 1.0
iso0_geo = TimeDependentBoundaryGeometryFromEquation(list(sym_variables[0:2]), [1.0, 0.0, x_init_gamma], "iso0", [sym_variables[-1]], [-vel_interface])

left_bnd_geo = BoundaryGeometryFromEquation(list(sym_variables[0:2]), [-1.0, 0.0, -left_bc_pos], "left_bnd")

right_bnd_geo = BoundaryGeometryFromEquation(list(sym_variables[0:2]), [1.0, 0.0, right_bc_pos], "right_bnd")


# === Jumps on interface ===
iso0_jump = GeneralBoundaryConditions([eqs_0, eqs_1], iso0_geo)
coord_to_subsituted = sym_variables[0]

# -> Dirichlet jumps on solutions
iso0_jump.addDirichletBCAndSubsituteBoundaryEq(DefaultSolutionTags.VAR_1, [jump_sol], coord_to_subsituted)

# -> Neumann jumps on fluxes
iso0_jump.addNeumannBCAndSubsituteBoundaryEq(DefaultSolutionTags.VAR_1, [jump_grad], coord_to_subsituted, ProjectionType.NORMAL)

# === Boundary conditions at the left boundary ===
left_bc = GeneralBoundaryConditions([eqs_0], left_bnd_geo)
left_bc.addDirichletBCAndSubsituteBoundaryEq(DefaultSolutionTags.VAR_1, [left_bc_sol_val], coord_to_subsituted)

# === Boundary conditions at the left boundary ===
right_bc = GeneralBoundaryConditions([eqs_1], right_bnd_geo)
right_bc.addNeumannBCAndSubsituteBoundaryEq(DefaultSolutionTags.VAR_1, [right_bc_grad_val], coord_to_subsituted, ProjectionType.NORMAL)

print(" ")
print("======== Printing list of conditions ========")
jumps_print = iso0_jump.getImposedConditions()
for i in jumps_print:
    print(sp.simplify(i))
jumps_print = left_bc.getImposedConditions()
for i in jumps_print:
    print(sp.simplify(i))
jumps_print = right_bc.getImposedConditions()
for i in jumps_print:
    print(sp.simplify(i))

#* -------------------------------------------------------------------------- *#
#* --- 4/ PROBLEM HANDLING ---
#* -------------------------------------------------------------------------- *#
my_problem = GeneralProblemHandler([eqs_0, eqs_1], [iso0_jump, left_bc, right_bc])
params_sol = my_problem.getParametricSolCoeff()

print(" ")
print("======== Printing values of unknown parameters ========")
for key, value in params_sol.items() :
    print ("%s =  %s"%(str(key), str(value)))
print("========================================================")

print(" ")
print("======== Printing solutions vectors of both phases ========")
print("*** Phase 0 ***")
for i in my_problem.getThisPhysicalModel(name_model_0).getSolVectorFromTags([],[], params_sol):
    print(i)
print("*** Phase 0 ***")
for i in my_problem.getThisPhysicalModel(name_model_1).getSolVectorFromTags([],[], params_sol):
    print(i)

#* -------------------------------------------------------------------------- *#
#* --- 5/ OUTPUTS PRINTING AND PLOTS ---
#* -------------------------------------------------------------------------- *#

# Print in files
output_path = output_absolute_path
my_problem.printMMSSourceTermInFile(output_path, OutputFileType.TEXT, False, True)
my_problem.printSolutionVectorInFile([], DefaultVarSetTags.MYVARSETTAG, output_path)

# Area over which to plot
y_dim = 0.1
time_plot = 0.7
pt_1 = [left_bc_pos, -y_dim]
pt_2 = [left_bc_pos, y_dim]
pt_3 = [right_bc_pos, y_dim]
pt_4 = [right_bc_pos, -y_dim]
iso0_plot_geo = my_problem.getThisBoundary("iso0").getBoundaryGeometry()
iso0_plot_geo.setBoundaryEquation(iso0_plot_geo.getBoundaryEquation().subs({sym_variables[-1]: time_plot}))
plot_area = PlotOver2DAreaWithStraightBoundaries(sym_variables[0:2], [iso0_plot_geo], [pt_1, pt_2, pt_3, pt_4])
line_plot = [BoundaryGeometryFromEquation(sym_variables[0:2], [0.0, 1.0, 0.0], "plot_line")]

u_plot = QuantityInfoForPlot(DefaultSolutionTags.VAR_1, False, ProjectionType.NOPROJECTION, time_plot)
grad_u_plot = QuantityInfoForPlot(DefaultSolutionTags.VAR_1, True, ProjectionType.NORMAL, time_plot)
quantities_plot = [u_plot, grad_u_plot]

# Sides of interface
side_interface_plot = dict()
side_interface_plot[name_model_0] = -1
side_interface_plot[name_model_1] = 1

my_problem.plotQuantitiesAlongOneDimLinesOverThisArea(quantities_plot, line_plot, plot_area, side_interface_plot, 100, 0)

my_problem.plotQuantitiesOverThisTwoDimArea(quantities_plot, plot_area, side_interface_plot, 30)
