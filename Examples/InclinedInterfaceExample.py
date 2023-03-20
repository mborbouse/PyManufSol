import sympy as sp
from PhysicalModels.ThermoPhysicalModels.EOS_factory import StiffenedGasEOS
from PhysicalModels.ThermoPhysicalModels.TransportProperties_factory import ConstantFluidTransportProperty
from PhysicalModels.CompressibleNavierStokesModels import CompressibleNavierStokesModels
from ManufSolution.ManufSolShape.PolarPolynomialManufSol import PolarPolynomialManufSol
from ManufSolution.ManufSolShape.UserDefinedManufSol import UserDefinedManufSol
from ManufSolution.ManufSolContainer.CompressibleFlowManufSolContainer import CompressibleFlowManufSolContainer
from ManufSolution.ManufSolContainer.GeneralManufSolContainer import GeneralManufSolContainerMultiplePhases
from Common.MMSTags import CompressibleFlowVarSetTags, FluidSolutionTags, FluidSolutionGradientTags, ProjectionType
from BoundaryConditions.BoundaryGeometry import BoundaryGeometryFromEquation, BoundaryGeometryFromPoints
from BoundaryConditions.GeneralBoundaryConditions import GeneralBoundaryConditions
from Common.MMSTags import CoordinatesSystemType, OutputFileType
from ProblemSolving.GeneralProblemHandler import GeneralProblemHandler, QuantityInfoForPlot
from Outputs.PlotOverArea import *

#* -------------------------------------------------------------------------- *#
#* --- 0/ DEFINITION OF PROBLEM PARAMETERS ---
#* -------------------------------------------------------------------------- *#
pressureNormalStressDecoupled = True
pressureHeatFluxDecoupled = False

# Domain geometry
L = [0.0,1.0]
H = [0.0,1.0]
x_dim = L[1]-L[0]
y_dim = H[1]-H[0]
y_interface_pos = H[0] + y_dim/3
domain_dim = 2

# Jumps conditions parameters
m_dot = 0.0 # using m_dot != 0 leads to a non-solvable system...
deltaTangentStress = 0.2
deltaTangentVel = -0.5
delta_p = -0.2
deltaTemp = 0.5
heatSource = 1.0

#* -------------------------------------------------------------------------- *#
#* --- 1/ MANUFACTURED SOLUTION ASSEMBLY ---
#* -------------------------------------------------------------------------- *#
sym_variables = sp.symbols('X, Y')

# Temperatures
T_0_user = 1.0*sp.cos(0.75*sp.pi*sym_variables[1]) + 2.45
T_0 = UserDefinedManufSol([], sym_variables, T_0_user)

T_coeff = list(sp.symbols('bT_1, bT_2'))
T_1_user = T_coeff[0]*sp.cos(5*sp.pi*sym_variables[1]/(4*y_dim)) + T_coeff[1] + 2.25
T_1 = UserDefinedManufSol(T_coeff, sym_variables, T_1_user)

# Horizontal velocities
u_0_user = 1.0*sp.cos(0.75*sp.pi*sym_variables[1]/y_dim) + 2.45
u_0 = UserDefinedManufSol([], sym_variables, u_0_user)

u_coeff = list(sp.symbols('bU_1, bU_2'))
u_1_user = u_coeff[0]*sp.cos(5*sp.pi*sym_variables[1]/(4*y_dim)) + u_coeff[1] + 2.25
u_1 = UserDefinedManufSol(u_coeff, sym_variables, u_1_user)

# Vertical velocities
v_0_user = 1.0*sp.cos(0.75*sp.pi*sym_variables[1]/y_dim) + 2.45
v_0 = UserDefinedManufSol([], sym_variables, v_0_user)

v_coeff = list(sp.symbols('bV_1, bV_2'))
v_1_user = v_coeff[0]*sp.cos(5*sp.pi*sym_variables[1]/(4*y_dim)) + v_coeff[1] + 2.25
v_1 = UserDefinedManufSol(v_coeff, sym_variables, v_1_user)

# Pressures
p_in = 1.1
y_top = H[1]
y_int = y_interface_pos
p_0_in = p_in
p_0_bnd = 1.1*p_0_in#1.2*p_0_in
p_1_in = p_0_in - delta_p
p_1_bnd = 0.9*p_1_in#0.8*p_1_in
a_0 = p_0_bnd
b_0 = 0.5*sp.pi*(1/y_top)
c_0 = (p_0_in - p_0_bnd) / (sp.cos(b_0*y_int))
b_1 = 1.0
c_1 = (p_1_in-p_1_bnd) / (sp.cos(b_1*y_int)-1.0)
a_1 = p_1_bnd-c_1
p_0_user = c_0*sp.cos(b_0*sym_variables[1])+a_0
p_0 = UserDefinedManufSol([], sym_variables, p_0_user)

p_coeff = []
if pressureNormalStressDecoupled:
    p_coeff = [sp.symbols('bP_1')]
    p_1_user = c_1*sp.cos(b_1*sym_variables[1])+a_1+p_coeff[0]
else:
    p_1_user = c_1*sp.cos(b_1*sym_variables[1])+a_1
p_1 = UserDefinedManufSol(p_coeff, sym_variables, p_1_user)

# Assembly into dict
manuf_sol_list_0 = [[p_0], [u_0, v_0], [T_0]]
manuf_sol_list_1 = [[p_1], [u_1, v_1], [T_1]]
manuf_sol_cont_0 = CompressibleFlowManufSolContainer(manuf_sol_list_0, domain_dim, sym_variables, CompressibleFlowVarSetTags.PRIMITIVE_PVT, [])
manuf_sol_cont_1 = CompressibleFlowManufSolContainer(manuf_sol_list_1, domain_dim, sym_variables, CompressibleFlowVarSetTags.PRIMITIVE_PVT, [])
dict_symb_sols = GeneralManufSolContainerMultiplePhases([manuf_sol_cont_0, manuf_sol_cont_1], [0, 1])

print(dict_symb_sols.getDicManufSolPerPhaseAndTag())

#* -------------------------------------------------------------------------- *#
#* --- 2/ DEFINITION OF THE GOVERNING EQUATIONS ---
#* -------------------------------------------------------------------------- *#

# Equation of states
gamma_0 = 1.4
Cv_0 = 2.5
p_inf_0 = 0.0
EOS_phase_0 = StiffenedGasEOS(gamma_0, Cv_0, p_inf_0, domain_dim)

gamma_1 = 1.6
Cv_1 = 2.0
p_inf_1 = 0.0
EOS_phase_1 = StiffenedGasEOS(gamma_1, Cv_1, p_inf_1, domain_dim) 

# Transport properties
mu_0 = 1.0
lambda_0 = 1.0
transp_prop_0 = ConstantFluidTransportProperty(mu_0, lambda_0)

mu_1 = 0.25
lambda_1 = 0.25
transp_prop_1 = ConstantFluidTransportProperty(mu_1, lambda_1)

# Equations models
name_model_0 = "phase0"
eqs_0 = CompressibleNavierStokesModels(name_model_0, manuf_sol_cont_0, domain_dim, EOS_phase_0, transp_prop_0, CoordinatesSystemType.CARTESIAN)

name_model_1 = "phase1"
eqs_1 = CompressibleNavierStokesModels(name_model_1, manuf_sol_cont_1, domain_dim, EOS_phase_1, transp_prop_1, CoordinatesSystemType.CARTESIAN)

#* -------------------------------------------------------------------------- *#
#* --- 3/ DEFINITION OF THE BOUNDARY AND JUMPS CONDITIONS (OPTIONAL!) ---
#* -------------------------------------------------------------------------- *#

# Definition of boundaries geometry
iso0_geo = BoundaryGeometryFromEquation(sym_variables, [0.0, -1.0, -y_interface_pos], "iso0")
coord_to_subsituted = sym_variables[1]

# === Jumps on interface ===
iso0_jump = GeneralBoundaryConditions([eqs_0, eqs_1], iso0_geo)

# -> Dirichlet jumps on solutions

# - Normal velocity
jump_v_n = 0.0
rho_0 = eqs_0.getSolTag(FluidSolutionTags.DENSITY)
rho_1 = eqs_1.getSolTag(FluidSolutionTags.DENSITY)
if m_dot != 0.0:
    jump_inv_rho = ((1.0/rho_0) - (1.0/rho_1))
    jump_v_n = m_dot * jump_inv_rho
iso0_jump.addDirichletBCAndSubsituteBoundaryEq(FluidSolutionTags.VELOCITY_VEC, [jump_v_n], coord_to_subsituted, ProjectionType.NORMAL)

# - Tangential velocity
iso0_jump.addDirichletBCAndSubsituteBoundaryEq(FluidSolutionTags.VELOCITY_VEC, [deltaTangentVel], coord_to_subsituted, ProjectionType.TANGENT)

# - Temperature
iso0_jump.addDirichletBCAndSubsituteBoundaryEq(FluidSolutionTags.TEMPERATURE, [deltaTemp], coord_to_subsituted)

# - Pressure
if pressureNormalStressDecoupled:
    jump_pressure = delta_p - m_dot * jump_v_n
    iso0_jump.addDirichletBCAndSubsituteBoundaryEq(FluidSolutionTags.PRESSURE, [jump_pressure], coord_to_subsituted)

# -> Neumann jumps on fluxes

# - Normal shear stress
p_0 = eqs_0.getSolTag(FluidSolutionTags.PRESSURE)
p_1 = eqs_1.getSolTag(FluidSolutionTags.PRESSURE) 
if pressureNormalStressDecoupled: # -[tau_nn] = 0.0 -> [p] = sigma*kappa - m_dot * [u_n]
    jump_tau_nn = 0.0
else: # -[tau_nn] = -m_dot*[u_n] - [p] + sigmaKappa
    jump_tau_nn = -m_dot * jump_v_n - (p_0-p_1) + delta_p
iso0_jump.addNeumannBCAndSubsituteBoundaryEq(FluidSolutionGradientTags.SHEARSTRESS, [-jump_tau_nn], coord_to_subsituted, ProjectionType.NORMALNORMAL)

# - Tangential shear stress
jump_tau_nt = -m_dot * deltaTangentVel + deltaTangentStress
iso0_jump.addNeumannBCAndSubsituteBoundaryEq(FluidSolutionGradientTags.SHEARSTRESS, [-jump_tau_nt], coord_to_subsituted, ProjectionType.NORMALTANGENT)

# - Heat flux minus viscous dissipation
E_0 = eqs_0.getSolTag(FluidSolutionTags.TOTALENERGY)
E_1 = eqs_1.getSolTag(FluidSolutionTags.TOTALENERGY)
vel_0 = eqs_0.getSolTag(FluidSolutionTags.VELOCITY_VEC)
vel_1 = eqs_1.getSolTag(FluidSolutionTags.VELOCITY_VEC)
n_gamma = iso0_jump.getBoundaryGeometry().getNormalToBoundary()
t_gamma = iso0_jump.getBoundaryGeometry().getTangentsToBoundary()
vel_0_n_gamma = sum([a*b for a,b in zip(vel_0,n_gamma)])
vel_1_n_gamma = sum([a*b for a,b in zip(vel_1,n_gamma)])
vel_0_t_gamma = sum([a*b for a,b in zip(vel_0,t_gamma[0])])
vel_1_t_gamma = sum([a*b for a,b in zip(vel_1,t_gamma[0])])
u_n_ave = 0.5 * (vel_0_n_gamma+vel_1_n_gamma)
overRho_ave = 0.5 * (1.0/rho_0 + 1.0/rho_1)
u_gamma_n = u_n_ave - m_dot * overRho_ave 
u_gamma_t = 0.5 * (vel_0_t_gamma+vel_1_t_gamma)   
if pressureHeatFluxDecoupled:
    jump_q_minus_tau = u_gamma_t * deltaTangentStress + heatSource
else:
    jump_q_minus_tau = -m_dot * (E_0-E_1) - (p_0*vel_0_n_gamma-p_1*vel_1_n_gamma) + delta_p * u_gamma_n + u_gamma_t * deltaTangentStress + heatSource
iso0_jump.addNeumannBCAndSubsituteBoundaryEq(FluidSolutionGradientTags.HEATFLUX_MINUS_VISCOUSDISSIPATION, [jump_q_minus_tau], coord_to_subsituted, ProjectionType.NORMAL)

print(" ")
print("======== Printing list of conditions ========")
jumps_print = iso0_jump.getImposedConditions()
for i in jumps_print:
    print(sp.simplify(i))

#* -------------------------------------------------------------------------- *#
#* --- 4/ PROBLEM HANDLING ---
#* -------------------------------------------------------------------------- *#

my_problem = GeneralProblemHandler([eqs_0, eqs_1], [iso0_jump])
params_sol = my_problem.getParametricSolCoeff()

tag_pvT = [FluidSolutionTags.PRESSURE, FluidSolutionTags.VELOCITY_X, FluidSolutionTags.VELOCITY_Y, FluidSolutionTags.TEMPERATURE]

print(" ")
print("======== Printing values of unknown parameters ========")
for key, value in params_sol.items() :
    print ("%s =  %s"%(str(key), str(value)))
print("========================================================")

print(" ")
print("======== Printing (p, v, T) solutions vectors of both phases ========")
print("*** Phase 0 ***")
for i in tag_pvT:
    print(my_problem.getThisPhysicalModel(name_model_0).getSolTag(i, params_sol))
print("*** Phase 1 ***")
for i in tag_pvT:
    print(my_problem.getThisPhysicalModel(name_model_1).getSolTag(i, params_sol))
print(" ")
print("======== Printing normal to boundaries ========")
print(my_problem.getThisBoundary("iso0").getBoundaryGeometry().getNormalToBoundary())

#* -------------------------------------------------------------------------- *#
#* --- 5/ CHANGE OF COORDINATES SYSTEM (OPTIONAL!) ---
#* -------------------------------------------------------------------------- *#

new_sym_variables = sp.symbols('x, y')
angle_rotation = -15
my_new_problem = my_problem.copyAndChangeCoordinatesSystem(list(new_sym_variables), CoordinatesSystemType.CARTESIAN, None, [[2,angle_rotation]])
params_sol_new = my_new_problem.getParametricSolCoeff()

print(" ")
print("======== Printing (p, v, T) solutions vectors in new rotated reference frame ========")
print("*** Phase 0 ***")
for i in tag_pvT:
    print(my_new_problem.getThisPhysicalModel(name_model_0).getSolTag(i, params_sol_new))
print("*** Phase 1 ***")
for i in tag_pvT:
    print(my_new_problem.getThisPhysicalModel(name_model_1).getSolTag(i, params_sol_new))

print(" ")
print("======== Printing normal to boundaries ========")
print(my_new_problem.getThisBoundary("iso0").getBoundaryGeometry().getBoundaryEquation())
new_normal = my_new_problem.getThisBoundary("iso0").getBoundaryGeometry().getNormalToBoundary()
print([sp.simplify(i) for i in new_normal])

#* -------------------------------------------------------------------------- *#
#* --- 6/ OUTPUTS PRINTING AND PLOTS ---
#* -------------------------------------------------------------------------- *#

# Print in files
output_path = "/Users/henneauxd/Softwares/myMMS_Solver/Examples/InclinedInterfaceOutput/"
my_new_problem.printMMSSourceTermInFile(output_path, OutputFileType.TEXT)
my_new_problem.printSolutionVectorInFile([], CompressibleFlowVarSetTags.PRIMITIVE_PVT, output_path)

# Area over which to plot
pt_1 = [L[0], H[0]]
pt_2 = [L[0], H[1]]
pt_3 = [L[1], H[1]]
pt_4 = [L[1], H[0]]
plot_area_cart = PlotOver2DAreaWithStraightBoundaries(new_sym_variables, [my_new_problem.getThisBoundary("iso0").getBoundaryGeometry()], [pt_1, pt_2, pt_3, pt_4])

plot_area_cart_original = PlotOver2DAreaWithStraightBoundaries(sym_variables, [my_problem.getThisBoundary("iso0").getBoundaryGeometry()], [pt_1, pt_2, pt_3, pt_4])

# 1-D lines along which to plot
nb_lines_plot = 4
coords_plot = np.linspace(L[0]+0.01, L[1]-0.01, nb_lines_plot)
lines_plot = []
for i in range(nb_lines_plot):
    lines_plot.append(BoundaryGeometryFromEquation(new_sym_variables, [1.0, 0.0, coords_plot[i]], "plot_line_"+str(i)))
lines_plot_original = []
for i in range(nb_lines_plot):
    lines_plot_original.append(BoundaryGeometryFromEquation(sym_variables, [1.0, 0.0, coords_plot[i]], "plot_line_"+str(i)))

# Quantities to plot
u_plot = QuantityInfoForPlot(FluidSolutionTags.VELOCITY_X)
v_plot = QuantityInfoForPlot(FluidSolutionTags.VELOCITY_Y)
vel_n_plot = QuantityInfoForPlot(FluidSolutionTags.VELOCITY_VEC, False, ProjectionType.NORMAL)
vel_t_plot = QuantityInfoForPlot(FluidSolutionTags.VELOCITY_VEC, False, ProjectionType.TANGENT)
T_plot = QuantityInfoForPlot(FluidSolutionTags.TEMPERATURE)
p_plot = QuantityInfoForPlot(FluidSolutionTags.PRESSURE)
rho_plot = QuantityInfoForPlot(FluidSolutionTags.DENSITY)
tau_nn_plot = QuantityInfoForPlot(FluidSolutionGradientTags.SHEARSTRESS, True, ProjectionType.NORMALNORMAL)
tau_nt_plot = QuantityInfoForPlot(FluidSolutionGradientTags.SHEARSTRESS, True, ProjectionType.NORMALTANGENT)
q_n_plot = QuantityInfoForPlot(FluidSolutionGradientTags.HEATFLUX, True, ProjectionType.NORMAL)
q_minus_tauu_n_plot = QuantityInfoForPlot(FluidSolutionGradientTags.HEATFLUX_MINUS_VISCOUSDISSIPATION, True, ProjectionType.NORMAL)
quantities_plot = [u_plot, v_plot, vel_n_plot, vel_t_plot, T_plot, p_plot, rho_plot, tau_nn_plot, tau_nt_plot, q_n_plot, q_minus_tauu_n_plot]

# Sides of interface
side_interface_plot = dict()
side_interface_plot[name_model_0] = -1
side_interface_plot[name_model_1] = 1

# Plot everything
# my_problem.plotQuantitiesAlongOneDimLinesOverThisArea(quantities_plot, lines_plot_original, plot_area_cart_original, side_interface_plot, 20, 1)

my_new_problem.plotQuantitiesAlongOneDimLinesOverThisArea(quantities_plot, lines_plot, plot_area_cart, side_interface_plot, 20, 1)

my_new_problem.plotQuantitiesOverThisTwoDimArea(quantities_plot, plot_area_cart, side_interface_plot, 30)
