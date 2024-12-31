SUBPROGRAM One (T_PV : T_PTC)
$COMMON G_GHI, A_PTC, alpha_PTC, sigma, epsilon_PV, A_PV, h_PTC, epsilon_PTC, T_sky, T_a

" ---- ---- ---- ---- Calculation of heat fluxes ---- ---- ---- ---- "

q_dot_sol_PTC = G_GHI * A_PTC * alpha_PTC

q_dot_rad_PV = sigma * epsilon_PV * A_PV * (T_PV^4 - T_PTC^4)
q_dot_rad_PTC = sigma * epsilon_PTC * A_PTC * (T_PTC^4 - T_sky^4)

q_dot_conv_front_PTC = A_PTC * h_PTC * (T_PTC - T_a)
q_dot_conv_back_PTC = A_PTC * h_PTC * (T_PTC - T_a)

" ---- ---- ---- ---- Estimation of T_PTC "

q_dot_sol_PTC + q_dot_rad_PV = q_dot_conv_front_PTC + q_dot_conv_back_PTC + q_dot_rad_PTC

END

" ****************************************************************************************************************** "

SUBPROGRAM Two (T_a, T_sky, T_in_HTF, T_PV, T_PTC : T_abs, T_sub, T_out_HTF, P_PV, q_dot_HTF)
$COMMON G_GHI, G_DNI, A_PV, A_abs, A_PTC, alpha_PV, alpha_abs, alpha_PTC, eta_opt, eta_PV, IAM_th, IAM_elec, CR_PTC, sigma, epsilon_PV, epsilon_abs, epsilon_PTC, m_dot_HTF, C_p_HTF, E, R_conv_PV, R_conv_abs, R_cond_abs, R_cond_sub

" ---- ---- ---- ---- Calculation of heat fluxes ---- ---- ---- ---- "

q_dot_sol_PV = G_DNI * A_PV * CR_PTC * eta_opt  * IAM_th * alpha_PV
q_dot_sol_abs = G_GHI * A_abs * alpha_abs

q_dot_rad_PV = sigma * epsilon_PV * A_PV * (T_PV^4 - T_PTC^4)
q_dot_rad_abs = sigma * epsilon_abs * A_abs * (T_abs^4 - T_sky^4)

q_dot_conv_PV = (T_PV - T_a) * R_conv_PV
q_dot_conv_abs = (T_abs - T_a) * R_conv_abs

q_dot_cond_abs_x_sub = (T_abs - T_sub) / (R_cond_abs + R_cond_sub)

q_dot_HTF = m_dot_HTF  * C_p_HTF * (T_out_HTF - T_in_HTF)
P_PV = G_DNI * A_PV * CR_PTC * eta_opt * IAM_elec * eta_PV

" ---- ---- ---- ---- System of three equations for: T_out_HTF, T_abs and T_sub "

q_dot_sol_PV - P_PV - q_dot_conv_PV - q_dot_rad_PV + q_dot_sol_abs - q_dot_conv_abs - q_dot_rad_abs = q_dot_HTF
q_dot_sol_abs - q_dot_conv_abs - q_dot_rad_abs = q_dot_cond_abs_x_sub

m_dot_HTF * C_p_HTF * (T_out_HTF - T_in_HTF) = E * m_dot_HTF  * C_p_HTF * (T_sub - T_in_HTF)

END

" ****************************************************************************************************************** "

SUBPROGRAM Three (T_abs, T_sub, T_in_HTF, T_out_HTF : T_estimated_PV)
$COMMON m_dot_HTF, C_p_HTF, R_cond_PV, R_cond_abs, R_cond_sub

" ---- ---- ---- ---- Calculation of heat fluxes ---- ---- ---- ---- "

q_dot_cond_PV_x_sub = (T_estimated_PV - T_sub) / (R_cond_PV + R_cond_sub)
q_dot_cond_abs_x_sub = (T_abs - T_sub) / (R_cond_abs + R_cond_sub)

q_dot_HTF = m_dot_HTF  * C_p_HTF * (T_out_HTF - T_in_HTF)

" ---- ---- ---- ---- Estimation of T_estimated_PV "

q_dot_cond_PV_x_sub + q_dot_cond_abs_x_sub = q_dot_HTF

END

" ****************************************************************************************************************** "

" ---- ---- ---- ---- Boundary conditions ---- ---- ---- ---- "

G_GHI =1000                               "Global horizontal irradiance"
G_DNI = 800                               "Direct normal irradiance"
u_air = 5                                 "Air velocity"

T_a = 25 + 273.15                         "Ambient temperature"
T_in_HTF = 70 + 273.15                    "HTF inlet temperature"
T_sky = 25 + 273.15                       "Sky temperature"
T_PV = 355.45                             "PV initial temperature"

" ---- ---- ---- ---- Geometry of the PTC and SRC-PVT ---- ---- ---- ---- "

D = 0.03                                  "Pipe diameter"
W_PV = 0.12                               "Width of the 2 faces of the SRC-PVT, in cross-section"
W_abs = 0.06                              "Width of the absorber, in cross-section"
W_PTC = 1.2                               "Width of PTC, in cross-section"

L_PTC = 10                                "Concentrator length"
L_SRC_PVT = L_PTC                         "SRC-PVT length"

A_PV = W_PV * L_SRC_PVT                   "PV area"
A_abs = W_abs * L_SRC_PVT                 "Absorber area"
A_ap = 1.2 * L_PTC                        "PTC aperture area"
A_PTC = 3 * L_PTC                         "PTC area"
CR_PTC = A_ap / A_PV                      "PTC concentration ratio"

" ---- ---- ---- ---- Design parameters ---- ---- ---- ---- "

alpha_PV = 0.97                           "PV admittance"
alpha_abs = 0.9                           "Absorber admittance"
alpha_PTC =0.03                           "PTC admittance"
eta_opt = 0.83                            "Optical efficiency"
IAM_th = 0.86                             "Angle of thermal incidence"
IAM_elec = 0.72                           "Angle of electric incidence"
epsilon_PV = 0.2                          "PV emissivity"
epsilon_abs = 0.2                         "Absorber emissivity"
epsilon_PTC = 0.3                         "PTC emissivity"
sigma = 5.67*10^(-8)                      "Stefan-Boltzmann constant"
P_air = 1.01325                           "Air pressure in bar"
P_HTF = 0.3119                            "Water pressure in bar"
m_dot_HTF = 0.15                          "Water mass flow"

th_PV = 0.003                             "PV thickness"
th_abs = 0.003                            "Absorber thickness"
k_PV = 50                                 "Coefficient of thermal conductivity of PV "
k_abs = 205                               "Coefficient of thermal conductivity of absorber"
k_sub = 250                               "Coefficient of thermal conductivity of substrate"

" ****************************************************************************************************************** "

" ---- ---- ---- ---- Calculation of thermal conductive resistances ---- ---- ---- ---- "

R_cond_PV = th_PV / (k_PV * A_PV)
R_cond_abs = th_abs / (k_abs * A_abs)
R_cond_sub = th_sub / (k_sub * A_PV)

" ---- ---- ---- ---- Calculation of convective heat resistances ---- ---- ---- ---- "

R_conv_PV = 1 / (h_PV * A_PV)
R_conv_abs = 1 / (h_abs * A_abs)

" ---- ---- Calculation of convection coefficients ---- ---- "

" ---- ---- Internal Convection "

" -- Thermal properties of water at 70 °C "
rho_water = density(Water,T=T_in_HTF,P=P_HTF)       "Density"
mu_water = viscosity(Water,T=T_in_HTF,P=P_HTF)      "Dynamic viscosity"
Pr_water = prandtl(Water,T=T_in_HTF,P=P_HTF)        "Prandtl number"
k_water = conductivity(Water,T=T_in_HTF,P=P_HTF)    "Thermal conductivity"

" -- Convection coefficient of the HTF "

Re_HTF=(rho_water*u_water*D)/mu_water               "Calculation of the Reynolds number"
N_HTF=0.023*Re_HTF^(4/5)*Pr_water^(2/5)             "Calculation of the Nusselt number"
h_HTF=(N_HTF*k_water)/D                             "Calculation of the convection coefficient"

" ---- ---- External Convection "

" -- Thermal properties of air at 25 °C "
rho_air = density(Air_ha,T=T_a,P=P_air)
mu_air = viscosity(Air_ha,T=T_a,P=P_air)
Pr_air = prandtl(Air_ha,T=T_a,P=P_air)
k_air = conductivity(Air_ha,T=T_a,P=P_air)

" -- Convection coefficient of the PV "
Re_PV=(rho_air*u_air*W_PV)/mu_air
N_PV=0.664*Pr_air^(1/3)*Re_PV^(1/2)
h_PV=(N_PV*k_air)/W_PV

" -- Convection coefficient of the Absorber "
Re_abs=(rho_air*u_air*W_abs)/mu_air
N_abs=0.664*Pr_air^(1/3)*Re_abs^(1/2)
h_abs=(N_abs*k_air)/W_abs

" -- Convection coefficient of the PTC "
Re_PTC=(rho_air*u_air*W_PTC)/mu_air
N_PTC=0.664*Pr_air^(1/3)*(Re_PTC)^(1/2)
h_PTC=(N_PTC*k_air)/W_PTC

" ****************************************************************************************************************** "

" ---- ---- ---- ---- Water (HTF) velocity clearing ---- ---- ---- ---- "

A_cs_pipe = (3.1416 * (D/2)^2)                       "Pipe cross-sectional area"
u_water = ( m_dot_HTF  ) / ( rho_water * A_cs_pipe ) "Water velocity"

" ****************************************************************************************************************** "

" ---- ---- ---- ---- Calculation of the NUT number ---- ---- ---- ---- "

C_p_HTF = cp(Water,T=T_in_HTF,P=P_HTF)               "Specific heat of water at 70 °C"
A_hx = 3.1416 * D * 10                               "Heat exchanger area"
NUT = ( ( 1 / ( ( 1 / h_HTF ) + R_cond_sub ) ) * A_hx ) / ( m_dot_HTF * ( C_p_HTF ) )  "NUT number"
E = 1-exp(-NUT)                                      "Effectiveness of heat transfer"

" ****************************************************************************************************************** "

" ---- ---- ---- ---- Calculation of PV efficiency as a function of PTC geometry ---- ---- ---- ---- "

eta_PV = 0.298 + 0.0142 * ( ln(CR_PTC) ) + ( - 0.000715 +0.0000697 * ( ln(CR_PTC) ) ) * (T_estimated_PV - 298)

" ****************************************************************************************************************** "

" ---- ---- ---- ---- Calculation of the thickness of the substrate using a circle inscribed in a triangle ---- ---- ---- ---- "

r_pipe = D / 2
h_triangle = ( sqrt(3) / 2 ) * W_abs
h_inscribed = h_triangle / 3

" ---- ---- ---- ---- Minimum and maximum thickness "

th_min = h_inscribed - r_pipe
th_max = h_triangle - r_pipe

" ---- ---- ---- ---- Average thickness for one face of SRC-PVT "

th_sub = (th_min + 2 * th_max) / 3                   "Substrate thickness"

" ****************************************************************************************************************** "

" ---- ---- ---- ---- Calculation of electrical and thermal efficiency ---- ---- ---- ---- "

eta_elec = ( P_PV ) / ( G_DNI * A_ap )
eta_th = ( q_dot_HTF ) / ( G_DNI * A_ap )

" ****************************************************************************************************************** "

" ---- ---- ---- ---- Obtaining the model temperatures ---- ---- ---- ---- "

CALL One (T_PV : T_PTC)
CALL Two (T_a, T_sky, T_in_HTF, T_PV, T_PTC : T_abs, T_sub, T_out_HTF, P_PV, q_dot_HTF)
CALL Three (T_abs, T_sub, T_in_HTF, T_out_HTF : T_estimated_PV) 

Tc_pv = T_estimated_PV - 273.15
Tc_abs = T_abs - 273.15
Tc_sub = T_sub - 273.15
Tc_out_HTF = T_out_HTF - 273.15