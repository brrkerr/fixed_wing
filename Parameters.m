% Creator: Burak ER 
% Aerosonde UAV Parameters from 
% Small Unmanned Aircraft Theory and Practice, Randal W. Beard, Timothy
% W. McLain

% Aerodinamic parameter tables are created from these derivatives 
% for different alpha and beta


clear all;
close all;
clc;

model_sampling_time = 1/1000;

mass = 13.5;
Jx = 0.8244;
Jy = 1.135;
Jz = 1.759;
Jxz = 0.1204;
Jxy = 0;
Jyz = 0;
Inertia_Tensor = [Jx -Jxy -Jxz;
                -Jxy Jy -Jyz;
                -Jxz -Jyz Jz];

wing_span = 0.55;
chord = 0.1894;
S = 2.8956;
rho_msl = 1.225;

v_body_initial = [30; 0; 0];
pqr_initial = [0; 0; 0];
euler_initial = [0; 0; 0];
quat_initial = eul2quat(euler_initial');

mat_fun_initial = zeros(9,1);
mat_fun_initial(1:3) = v_body_initial;
mat_fun_initial(4:6) = pqr_initial;
mat_fun_initial(7:9) = euler_initial;


da_lim = deg2rad(45);
de_lim = deg2rad(45);
dr_lim = deg2rad(45);
dth_lim = 100;


h_ref = 0;
lat_ref_deg = 40;
long_ref_deg = 20;

%% 
Motor.Sprop = 0.2027;
Motor.k_motor = 80;
Motor.k_Tp = 0;
Motor.k_omega = 0;
Motor.Cprop = 1;

%% AERODYNAMIC COEFFICIENTS
CL_0 = 0.28;
CL_alpha = 3.45;
CL_p = 0;
CL_q = 0;
CL_r = 0;
CL_da = 0;
CL_de = -0.36;
CL_dr = 0;
CL_alpha_dot = 0;
CL_u = 0;

%
CY_0 = 0;
CY_beta = -0.98;
CY_p = 0;
CY_q = 0;
CY_r = 0;
CY_da = 0;
CY_de = 0;
CY_dr = -0.17;
CY_alpha_dot = 0;
CY_u = 0;

%
CD_0 = 0.03;
CD_alpha = 0.30;
CD_p = 0.0437;
CD_q = 0;
CD_r = 0;
CD_da = 0;
CD_de = 0;
CD_dr = 0;
CD_alpha_dot = 0;
CD_u = 0;

%
Cl_0 = 0;
Cl_beta = -0.12;
Cl_p = -0.26;
Cl_q = 0;
Cl_r = 0.14;
Cl_da = 0.08;
Cl_de = 0;
Cl_dr = 0.105;
Cl_alpha_dot = 0;
Cl_u = 0;

%
Cm_0 = -0.02338;
Cm_alpha = -0.38;
Cm_p = 0;
Cm_q = -3.6;
Cm_r = 0;
Cm_da = 0;
Cm_de = -0.5;
Cm_dr = 0;
Cm_alpha_dot = 0;
Cm_u = 0;

Cn_0 = 0;
Cn_beta = 0.25;
Cn_p = 0.022;
Cn_q = 0;
Cn_r = -0.35;
Cn_da = 0.06;
Cn_de = 0;
Cn_dr = -0.032;
Cn_alpha_dot = 0;
Cn_u = 0;

run("create_table.m");
close all
%% AERODYNAMIC COEFFICIENTS

Aerodynamics.Alpha = Alpha;
Aerodynamics.Beta = Beta;
Aerodynamics.da = da;
Aerodynamics.de = de;
Aerodynamics.dr = dr;

Aerodynamics.CL_alpha_beta = CL_alpha_beta;
Aerodynamics.CLp = CLp;
Aerodynamics.CLq = CLq;
Aerodynamics.CLr = CLr;
Aerodynamics.CLda = CLda;
Aerodynamics.CLde = CLde;
Aerodynamics.CLdr = CLdr;


%
Aerodynamics.CY_alpha_beta = CY_alpha_beta;
Aerodynamics.CYp = CYp;
Aerodynamics.CYq = CYq;
Aerodynamics.CYr = CYr;
Aerodynamics.CYda = CYda;
Aerodynamics.CYde = CYde;
Aerodynamics.CYdr = CYdr;

%
Aerodynamics.CD_alpha_beta = CD_alpha_beta;
Aerodynamics.CDp = CDp;
Aerodynamics.CDq = CDq;
Aerodynamics.CDr = CDr;
Aerodynamics.CDda = CDda;
Aerodynamics.CDde = CDde;
Aerodynamics.CDdr = CDdr;

%
Aerodynamics.Cl_alpha_beta = Cl_alpha_beta;
Aerodynamics.Clp = Clp;
Aerodynamics.Clq = Clq;
Aerodynamics.Clr = Clr;
Aerodynamics.Clda = Clda;
Aerodynamics.Clde = Clde;
Aerodynamics.Cldr = Cldr;

%
Aerodynamics.Cm_alpha_beta = Cm_alpha_beta;
Aerodynamics.Cmp = Cmp;
Aerodynamics.Cmq = Cmq;
Aerodynamics.Cmr = Cmr;
Aerodynamics.Cmda = Cmda;
Aerodynamics.Cmde = Cmde;
Aerodynamics.Cmdr = Cmdr;

%
Aerodynamics.Cn_alpha_beta = Cn_alpha_beta;
Aerodynamics.Cnp = Cnp;
Aerodynamics.Cnq = Cnq;
Aerodynamics.Cnr = Cnr;
Aerodynamics.Cnda = Cnda;
Aerodynamics.Cnde = Cnde;
Aerodynamics.Cndr = Cndr;

save aircraft_params
