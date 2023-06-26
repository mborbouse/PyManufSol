import os
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

current_directory = os.getcwd()
output_folder_name = "Examples/StaticBubbleOutput/"
output_absolute_path = os.path.join(current_directory,output_folder_name)

#* -------------------------------------------------------------------------- *#
#* --- 0/ DEFINITION OF PROBLEM PARAMETERS ---
#* -------------------------------------------------------------------------- *#
Inviscid = False
CstPressure = True
pressureNormalStressDecoupled = True
pressureHeatFluxDecoupled = True

varySigma = False
# Geometry definition
domain_dim = 2
L = [0.0,1.0]
H = [0.0,1.0]
domain_length = L[1]-L[0]
center = [0.0, 0.0]
radius = 0.2#0.385
diameter = 2.0*radius
mesh = 16
delta_x = domain_length/mesh
p_order = 3

# Properties
if Inviscid:
    kappa = 0.0
else:
    kappa = 1.0
Cv = 2.5
p_inf = 0.0

sigma = 1.0*radius
mu = 1.0

# Refs values
vel_ref = 0.0
p_ref = 1.0
rho_ref = 1.0
T_ref = 1.0
EOS_interm = StiffenedGasEOS(0.0, Cv, p_inf, 2)
gamma = EOS_interm.getGammaFromPTRho(p_ref, T_ref, rho_ref)
print("==>> gamma = %2.16f"%gamma)

# Non-dim numbers
if Inviscid:
    La = math.inf 
    T_nu = math.inf
    mu = 0.0
else:
    La = 12000 #* <<<<==== TO BE SET
    if varySigma:
        sigma = (La * mu * mu) / (rho_ref*diameter)
    else:
        mu = math.sqrt((rho_ref*diameter*sigma)/La)
    nu = mu/rho_ref
    T_nu = (diameter**2) / nu

U_sigma = math.sqrt((sigma)/(rho_ref*diameter))
T_sigma = math.sqrt((rho_ref*diameter**3)/sigma)
delta_t = math.sqrt((rho_ref*(delta_x/(p_order+1))**3)/(math.pi*sigma))

# Jumps
m_dot = 0.0 # using m_dot != 0 leads to a non-solvable system...
deltaTangentStress = 0.0
deltaTangentVel = 0.0
deltaTemp = 0.0
heatSource = 0.0
jump_p = sigma / radius

print("========= INPUTS INFOS =========")
print("-> Viscosity = %2.16f"%mu)
print("-> Sigma = %2.16f"%sigma)
print("-> T_nu = %2.16f"%T_nu)
print("-> T_sigma = %2.16f"%T_sigma)
print("-> U_sigma = %2.16f"%U_sigma)
print("-> delta_t = %2.16f"%delta_t)
print("-> Delta p = %2.16f"%jump_p)
print("-> La = %2.3f"%(sigma*rho_ref*diameter/(mu*mu)))
print("-> Oh = %2.3f"%(1/math.sqrt(La)))
print("================================")
print("")

#* -------------------------------------------------------------------------- *#
#* --- 1/ MANUFACTURED SOLUTION ASSEMBLY ---
#* -------------------------------------------------------------------------- *#

sym_variables = sp.symbols('r, theta')
known_ind = [3, 4]

# # Density
# sym_param_rho_0 = []#[sp.symbols('aR')]
# rho_0_user = rho_ref #+ sym_param_rho_0[0]
# rho_0_user = sp.Add(rho_0_user, sp.Mul(0.0,sym_variables[0],evaluate=False), evaluate=False)
# rho_0 = UserDefinedManufSol(sym_param_rho_0, sym_variables, rho_0_user)

# # sym_param_rho_1 = [sp.symbols('bR')]
# # rho_1_user = rho_ref + sym_param_rho_1[0]
# # rho_1_user = sp.Add(rho_1_user, sp.Mul(0.0,sym_variables[0],evaluate=False), evaluate=False)
# # rho_1 = UserDefinedManufSol(sym_param_rho_1, sym_variables, rho_1_user)

# sym_param_T1 = sp.symbols('bT_2, bT_3, bT_4')
# rho_1 = PolarPolynomialManufSol(sym_param_T1, sym_variables, [0, 0], known_ind, 4)

# # Radial velocities
# sym_param_u_0 = []#[sp.symbols('aU')]
# u_0_user = vel_ref# + sym_param_u_0[0]
# u_0_user = sp.Add(u_0_user, sp.Mul(0.0,sym_variables[0],evaluate=False), evaluate=False)
# u_0 = UserDefinedManufSol(sym_param_u_0, sym_variables, u_0_user)

# # sym_param_u_1 = [sp.symbols('bU')]
# # u_1_user = vel_ref + sym_param_u_1[0]
# # u_1_user = sp.Add(u_1_user, sp.Mul(0.0,sym_variables[0],evaluate=False), evaluate=False)
# # u_1 = UserDefinedManufSol(sym_param_u_1, sym_variables, u_1_user)

# sym_param_u_r_1 = sp.symbols('bU_2, bU_3, bU_4')
# u_1 = PolarPolynomialManufSol(sym_param_u_r_1, sym_variables, [0, 0],known_ind, 4)

# # Tangential velocities
# sym_param_v_0 = []#[sp.symbols('aV')]
# v_0_user = vel_ref# + sym_param_v_0[0]
# v_0_user = sp.Add(v_0_user, sp.Mul(0.0,sym_variables[0],evaluate=False), evaluate=False)
# v_0 = UserDefinedManufSol(sym_param_v_0, sym_variables, v_0_user)

# # sym_param_v_1 = [sp.symbols('bV')]
# # v_1_user = vel_ref + sym_param_v_1[0]
# # v_1_user = sp.Add(v_1_user, sp.Mul(0.0,sym_variables[0],evaluate=False), evaluate=False)
# # v_1 = UserDefinedManufSol(sym_param_v_1, sym_variables, v_1_user)

# sym_param_u_theta_1 = sp.symbols('bV_2, bV_3, bV_4')
# v_1 = PolarPolynomialManufSol(sym_param_u_theta_1, sym_variables, [0, 0], known_ind, 4)

# # Pressures
# sym_param_p_0 = []#sp.symbols('aP_1, aP_2')
# # p_0_user = p_ref + sym_param_p_0[0] + sym_param_p_0[1]*sp.Pow(sym_variables[0]/radius,2)
# p_0_user = p_ref + sp.Pow((sym_variables[0])/radius,3)
# p_0 = UserDefinedManufSol(sym_param_p_0, sym_variables, p_0_user)

# sym_param_p_1 = [sp.symbols('bP_1')]
# # p_1_user = p_ref + sym_param_p_1[0] + sym_param_p_1[1]*sp.Pow((sym_variables[0]-radius)/radius,2)
# p_1_user = p_ref + sym_param_p_1[0]*sp.Pow((sym_variables[0]-0.5*radius)/radius,3)
# p_1 = UserDefinedManufSol(sym_param_p_1, sym_variables, p_1_user)



# Density
sym_param_rho_0 = []#[sp.symbols('aR')]
rho_0_user = rho_ref #+ sym_param_rho_0[0]
rho_0_user = sp.Add(rho_0_user, sp.Mul(0.0,sym_variables[0],evaluate=False), evaluate=False)
rho_0 = UserDefinedManufSol(sym_param_rho_0, sym_variables, rho_0_user)

sym_param_rho_1 = []#[sp.symbols('aR')]
rho_1_user = rho_ref #+ sym_param_rho_0[0]
rho_1_user = sp.Add(rho_1_user, sp.Mul(0.0,sym_variables[0],evaluate=False), evaluate=False)
rho_1 = UserDefinedManufSol(sym_param_rho_1, sym_variables, rho_1_user)

# Radial velocities
sym_param_u_0 = []#[sp.symbols('aU')]
u_0_user = vel_ref# + sym_param_u_0[0]
u_0_user = sp.Add(u_0_user, sp.Mul(0.0,sym_variables[0],evaluate=False), evaluate=False)
u_0 = UserDefinedManufSol(sym_param_u_0, sym_variables, u_0_user)

sym_param_u_1 = []#[sp.symbols('aU')]
u_1_user = vel_ref# + sym_param_u_0[0]
u_1_user = sp.Add(u_1_user, sp.Mul(0.0,sym_variables[0],evaluate=False), evaluate=False)
u_1 = UserDefinedManufSol(sym_param_u_1, sym_variables, u_1_user)

# Tangential velocities
sym_param_v_0 = []#[sp.symbols('aV')]
v_0_user = vel_ref# + sym_param_v_0[0]
v_0_user = sp.Add(v_0_user, sp.Mul(0.0,sym_variables[0],evaluate=False), evaluate=False)
v_0 = UserDefinedManufSol(sym_param_v_0, sym_variables, v_0_user)

sym_param_v_1 = []#[sp.symbols('aV')]
v_1_user = vel_ref# + sym_param_v_0[0]
v_1_user = sp.Add(v_1_user, sp.Mul(0.0,sym_variables[0],evaluate=False), evaluate=False)
v_1 = UserDefinedManufSol(sym_param_v_1, sym_variables, v_1_user)

# Pressures
sym_param_p_0 = []#sp.symbols('aP_1, aP_2')
if CstPressure:
    p_0_user = sp.Add(p_ref + jump_p, sp.Mul(0.0,sym_variables[0],evaluate=False), evaluate=False)
else:
    p_0_user = p_ref + jump_p + sp.Pow((sym_variables[0]-radius)/radius,3) 
p_0 = UserDefinedManufSol(sym_param_p_0, sym_variables, p_0_user)

sym_param_p_1 = []
if CstPressure:
    p_1_user = sp.Add(p_ref, sp.Mul(0.0,sym_variables[0],evaluate=False), evaluate=False)
else:
    p_1_user = p_ref + sp.Pow((sym_variables[0]-radius)/radius,3)
p_1 = UserDefinedManufSol(sym_param_p_1, sym_variables, p_1_user)

# Assembly into dict
manuf_sol_list_0 = [[p_0], [u_0, v_0], [rho_0]]
manuf_sol_list_1 = [[p_1], [u_1, v_1], [rho_1]]
manuf_sol_cont_0 = CompressibleFlowManufSolContainer(manuf_sol_list_0, domain_dim, sym_variables, CompressibleFlowVarSetTags.PRIMITIVE_PVRHO, [])
manuf_sol_cont_1 = CompressibleFlowManufSolContainer(manuf_sol_list_1, domain_dim, sym_variables, CompressibleFlowVarSetTags.PRIMITIVE_PVRHO, [])
dict_symb_sols = GeneralManufSolContainerMultiplePhases([manuf_sol_cont_0, manuf_sol_cont_1], [0, 1])

print(dict_symb_sols.getDicManufSolPerPhaseAndTag())

#* -------------------------------------------------------------------------- *#
#* --- 2/ DEFINITION OF THE GOVERNING EQUATIONS ---
#* -------------------------------------------------------------------------- *#

# Equation of states
EOS_phase_0 = StiffenedGasEOS(gamma, Cv, p_inf, domain_dim)
EOS_phase_1 = StiffenedGasEOS(gamma, Cv, p_inf, domain_dim) 

# Transport properties
transp_prop_0 = ConstantFluidTransportProperty(mu, kappa)
transp_prop_1 = ConstantFluidTransportProperty(mu, kappa)

# Equations models
name_model_0 = "phase0"
eqs_0 = CompressibleNavierStokesModels(name_model_0, manuf_sol_cont_0, domain_dim, EOS_phase_0, transp_prop_0, CoordinatesSystemType.CYLINDRICAL)

name_model_1 = "phase1"
eqs_1 = CompressibleNavierStokesModels(name_model_1, manuf_sol_cont_1, domain_dim, EOS_phase_1, transp_prop_1, CoordinatesSystemType.CYLINDRICAL)

#* -------------------------------------------------------------------------- *#
#* --- 3/ DEFINITION OF THE BOUNDARY AND JUMPS CONDITIONS (OPTIONAL!) ---
#* -------------------------------------------------------------------------- *#

# Definition of boundaries geometry
iso0_geo = BoundaryGeometryFromEquation(sym_variables, [1.0, 0.0, radius], "iso0")
far_field_geo = BoundaryGeometryFromEquation(sym_variables, [1.0, 0.0, domain_length], "farField")
coord_to_subsituted = sym_variables[0]

# # === Jumps on interface ===
iso0_jump = GeneralBoundaryConditions([eqs_0, eqs_1], iso0_geo)

# # -> Dirichlet jumps on solutions

# # - Normal velocity
# jump_v_n = 0.0
# rho_0 = eqs_0.getSolTag(FluidSolutionTags.DENSITY)
# rho_1 = eqs_1.getSolTag(FluidSolutionTags.DENSITY)
# if m_dot != 0.0:
#     jump_inv_rho = ((1.0/rho_0) - (1.0/rho_1))
#     jump_v_n = m_dot * jump_inv_rho
# iso0_jump.addDirichletBCAndSubsituteBoundaryEq(FluidSolutionTags.VELOCITY_VEC, [jump_v_n], coord_to_subsituted, ProjectionType.NORMAL)

# # - Tangential velocity
# iso0_jump.addDirichletBCAndSubsituteBoundaryEq(FluidSolutionTags.VELOCITY_VEC, [deltaTangentVel], coord_to_subsituted, ProjectionType.TANGENT)

# # - Temperature
# iso0_jump.addDirichletBCAndSubsituteBoundaryEq(FluidSolutionTags.DENSITY, [deltaTemp], coord_to_subsituted)

# # - Pressure
# if pressureNormalStressDecoupled:
#     jump_pressure = jump_p - m_dot * jump_v_n
#     iso0_jump.addDirichletBCAndSubsituteBoundaryEq(FluidSolutionTags.PRESSURE, [jump_pressure], coord_to_subsituted)

# # -> Neumann jumps on fluxes

# # - Normal shear stress
# p_0 = eqs_0.getSolTag(FluidSolutionTags.PRESSURE)
# p_1 = eqs_1.getSolTag(FluidSolutionTags.PRESSURE) 
# if pressureNormalStressDecoupled: # -[tau_nn] = 0.0 -> [p] = sigma*kappa - m_dot * [u_n]
#     jump_tau_nn = 0.0
# else: # -[tau_nn] = -m_dot*[u_n] - [p] + sigmaKappa
#     jump_tau_nn = -m_dot * jump_v_n - (p_0-p_1) + jump_p
# iso0_jump.addNeumannBCAndSubsituteBoundaryEq(FluidSolutionGradientTags.SHEARSTRESS, [-jump_tau_nn], coord_to_subsituted, ProjectionType.NORMALNORMAL)

# # - Tangential shear stress
# jump_tau_nt = -m_dot * deltaTangentVel + deltaTangentStress
# iso0_jump.addNeumannBCAndSubsituteBoundaryEq(FluidSolutionGradientTags.SHEARSTRESS, [-jump_tau_nt], coord_to_subsituted, ProjectionType.NORMALTANGENT)

# # - Heat flux minus viscous dissipation
# E_0 = eqs_0.getSolTag(FluidSolutionTags.TOTALENERGY)
# E_1 = eqs_1.getSolTag(FluidSolutionTags.TOTALENERGY)
# vel_0 = eqs_0.getSolTag(FluidSolutionTags.VELOCITY_VEC)
# vel_1 = eqs_1.getSolTag(FluidSolutionTags.VELOCITY_VEC)
# n_gamma = iso0_jump.getBoundaryGeometry().getNormalToBoundary()
# t_gamma = iso0_jump.getBoundaryGeometry().getTangentsToBoundary()
# vel_0_n_gamma = sum([a*b for a,b in zip(vel_0,n_gamma)])
# vel_1_n_gamma = sum([a*b for a,b in zip(vel_1,n_gamma)])
# vel_0_t_gamma = sum([a*b for a,b in zip(vel_0,t_gamma[0])])
# vel_1_t_gamma = sum([a*b for a,b in zip(vel_1,t_gamma[0])])
# u_n_ave = 0.5 * (vel_0_n_gamma+vel_1_n_gamma)
# overRho_ave = 0.5 * (1.0/rho_0 + 1.0/rho_1)
# u_gamma_n = u_n_ave - m_dot * overRho_ave 
# u_gamma_t = 0.5 * (vel_0_t_gamma+vel_1_t_gamma)   
# if pressureHeatFluxDecoupled:
#     jump_q_minus_tau = u_gamma_t * deltaTangentStress + heatSource
# else:
#     jump_q_minus_tau = -m_dot * (E_0-E_1) - (p_0*vel_0_n_gamma-p_1*vel_1_n_gamma) + jump_p * u_gamma_n + u_gamma_t * deltaTangentStress + heatSource
# # iso0_jump.addNeumannBCAndSubsituteBoundaryEq(FluidSolutionGradientTags.HEATFLUX_MINUS_VISCOUSDISSIPATION, [jump_q_minus_tau], coord_to_subsituted, ProjectionType.NORMAL)
# iso0_jump.addNeumannBCAndSubsituteBoundaryEq(FluidSolutionGradientTags.HEATFLUX, [jump_q_minus_tau], coord_to_subsituted, ProjectionType.NORMAL)

# # === Boundary conditions in the far field ===
far_field_bc = GeneralBoundaryConditions([eqs_1], far_field_geo)

# far_field_bc.addDirichletBCAndSubsituteBoundaryEq(FluidSolutionTags.DENSITY, [rho_ref], coord_to_subsituted)
# far_field_bc.addDirichletBCAndSubsituteBoundaryEq(FluidSolutionTags.VELOCITY_X, [vel_ref], coord_to_subsituted)
# far_field_bc.addDirichletBCAndSubsituteBoundaryEq(FluidSolutionTags.VELOCITY_Y, [vel_ref], coord_to_subsituted)

# print(" ")
# print("======== Printing list of conditions ========")
# jumps_print = iso0_jump.getImposedConditions()
# for i in jumps_print:
#     print(sp.simplify(i))
# jumps_print = far_field_bc.getImposedConditions()
# for i in jumps_print:
#     print(sp.simplify(i))
# print("===============================================")

#* -------------------------------------------------------------------------- *#
#* --- 4/ PROBLEM HANDLING ---
#* -------------------------------------------------------------------------- *#

my_problem = GeneralProblemHandler([eqs_0, eqs_1], [iso0_jump, far_field_bc])
# my_problem = GeneralProblemHandler([eqs_0, eqs_1], [iso0_geo, far_field_geo])
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
center_new_domain = [0.0, 0.0]
my_new_problem = my_problem.copyAndChangeCoordinatesSystem(list(new_sym_variables), CoordinatesSystemType.CARTESIAN, center_new_domain, [[2,0.0]])
params_sol_new = my_new_problem.getParametricSolCoeff()

print(" ")
print("======== Printing (p, v, T) solutions vectors in new Cartesian reference frame ========")
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
output_path = output_absolute_path
my_new_problem.printMMSSourceTermInFile(output_path, OutputFileType.TEXT)
my_new_problem.printSolutionVectorInFile([], CompressibleFlowVarSetTags.PRIMITIVE_PVT, output_path)

# Area over which to plot
R_far = domain_length
plot_area_cylind = PlotOver2DCircularAreaCenteredAtOrigin(sym_variables, [my_problem.getThisBoundary("iso0").getBoundaryGeometry()], R_far)

pt_1 = [center_new_domain[0]-R_far, center_new_domain[1]-R_far]
pt_2 = [center_new_domain[0]-R_far, center_new_domain[1]+R_far]
pt_3 = [center_new_domain[0]+R_far, center_new_domain[1]+R_far]
pt_4 = [center_new_domain[0]+R_far, center_new_domain[1]-R_far]
plot_area_cart = PlotOver2DAreaWithStraightBoundaries(new_sym_variables, [my_new_problem.getThisBoundary("iso0").getBoundaryGeometry()], [pt_1, pt_2, pt_3, pt_4])

# 1-D lines along which to plot
nb_lines_plot = 4
coords_plot = np.linspace(0.0, 360.0, nb_lines_plot)
lines_plot = []
for i in range(nb_lines_plot):
    lines_plot.append(BoundaryGeometryFromEquation(sym_variables, [0.0, 1.0, coords_plot[i]], "plot_line_"+str(i)))


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
# quantities_plot = [p_plot, T_plot, rho_plot]

# Sides of interface
side_interface_plot = dict()
side_interface_plot[name_model_0] = -1
side_interface_plot[name_model_1] = 1

# Plot everything
my_problem.plotQuantitiesAlongOneDimLinesOverThisArea(quantities_plot, lines_plot, plot_area_cylind, side_interface_plot, 50, 0)

# my_problem.plotQuantitiesOverThisTwoDimArea(quantities_plot, plot_area_cylind, side_interface_plot, 10)

# my_new_problem.plotQuantitiesOverThisTwoDimArea(quantities_plot, plot_area_cart, side_interface_plot, 20)

# my_new_problem.plotQuantitiesAlongOneDimLinesOverThisArea(quantities_plot, lines_plot, plot_area_cart, side_interface_plot, 20)