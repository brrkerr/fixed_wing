function [X_dot] = aircraft_model(X,u, params)

X_dot = zeros(9,1);

v_body = X(1:3);
w_body = X(4:6);
euler = wrapToPi(X(7:9));

da_ = u(1);
de_ = u(2);
dr_ = u(3);
dth_ = u(4);

u = v_body(1);
v = v_body(2);
w = v_body(3);

phi = euler(1);
theta = euler(2);
psi = euler(3);

Alpha_ = atan2(w, u);
Va = sqrt(u^2 + v^2 + w^2);
Beta_ = asin(v / Va);


%%

DCM_e2b = [cos(theta)*cos(psi), cos(theta)*sin(psi), -sin(theta);
    sin(phi)*sin(theta)*cos(psi)-cos(phi)*sin(psi), sin(phi)*sin(theta)*sin(psi)+cos(phi)*cos(psi), sin(phi)*cos(theta);
    cos(phi)*sin(theta)*cos(psi)+sin(phi)*sin(psi), cos(phi)*sin(theta)*sin(psi)-sin(phi)*cos(psi), cos(phi)*cos(theta)];

DCM_b2s = [cos(Alpha_) , 0 , sin(Alpha_);
            0 , 1 , 0;
            -sin(Alpha_) , 0 , cos(Alpha_)];

DCM_s2w = [cos(Beta_) , sin(Beta_) , 0;
            -sin(Beta_) , cos(Beta_) , 0;
            0 , 0 , 1];

DCM_b2w = DCM_s2w*DCM_b2s;

Rot_pqr2euler = [1 , sin(phi)*tan(theta) , cos(phi)*tan(theta);
                0 , cos(phi) , sin(phi);
                0 , sin(phi)/cos(theta) , cos(phi)/cos(theta)];

%% Parameters

mass = params.mass;
Jx = params.Jx;
Jy = params.Jy;
Jz = params.Jz;
Jxz = params.Jxz;
Jxy = params.Jxy;
Jyz = params.Jyz;
Inertia_Tensor = [Jx -Jxy -Jxz;
                -Jxy Jy -Jyz;
                -Jxz -Jyz Jz];

wing_span = params.wing_span;
chord = params.chord;
S = params.S;
rho_msl = params.rho_msl;

da_lim = params.da_lim;
de_lim = params.de_lim;
dr_lim = params.dr_lim;
dth_lim = params.dth_lim;

if da_ >= da_lim
    da_ = da_lim;
elseif da_ <= -da_lim
    da_ = da_lim;
else
    da_ = da_;
end

if de_ >= de_lim
    de_ = de_lim;
elseif de_ <= -de_lim
    de_ = de_lim;
else
    de_ = de_;
end

if dr_ >= dr_lim
    dr_ = dr_lim;
elseif dr_ <= -dr_lim
    dr_ = dr_lim;
else
    dr_ = dr_;
end

if dth_ >= dth_lim
    dth_ = dth_lim;
elseif dth_ <= -dth_lim
    dth_ = dth_lim;
else
    dth_ = dth_;
end

Sprop = params.Motor.Sprop;
k_motor = params.Motor.k_motor;
k_Tp = params.Motor.k_Tp;
k_omega = params.Motor.k_omega;
Cprop = params.Motor.Cprop;

Alpha = params.Aerodynamics.Alpha;
Beta = params.Aerodynamics.Beta;
da = params.Aerodynamics.da;
de = params.Aerodynamics.de;
dr = params.Aerodynamics.dr;
%% AERODYNAMIC COEFFICIENTS
CL_ab = interp2(Beta, Alpha, params.Aerodynamics.CL_alpha_beta, Beta_, Alpha_);
CL_p = interp1(Alpha, params.Aerodynamics.CLp, Alpha_);
CL_q = interp1(Alpha, params.Aerodynamics.CLq, Alpha_);
CL_r = interp1(Alpha, params.Aerodynamics.CLr, Alpha_);
CL_da = interp2(Alpha, da, params.Aerodynamics.CLda, Alpha_, da_);
CL_de = interp2(Alpha, de, params.Aerodynamics.CLde, Alpha_, de_);
CL_dr = interp2(Beta, dr, params.Aerodynamics.CLdr, Alpha_, dr_);


%
CY_ab = interp2(Beta, Alpha, params.Aerodynamics.CY_alpha_beta, Beta_, Alpha_);
CY_p = interp1(Alpha, params.Aerodynamics.CYp, Alpha_);
CY_q = interp1(Alpha, params.Aerodynamics.CYq, Alpha_);
CY_r = interp1(Alpha, params.Aerodynamics.CYr, Alpha_);
CY_da = interp2(Alpha, da, params.Aerodynamics.CYda, Alpha_, da_);
CY_de = interp2(Alpha, de, params.Aerodynamics.CYde, Alpha_, de_);
CY_dr = interp2(Beta, dr, params.Aerodynamics.CYdr, Alpha_, dr_);

%
CD_ab = interp2(Beta, Alpha, params.Aerodynamics.CD_alpha_beta, Beta_, Alpha_);
CD_p = interp1(Alpha, params.Aerodynamics.CDp, Alpha_);
CD_q = interp1(Alpha, params.Aerodynamics.CDq, Alpha_);
CD_r = interp1(Alpha, params.Aerodynamics.CDr, Alpha_);
CD_da = interp2(Alpha, da, params.Aerodynamics.CDda, Alpha_, da_);
CD_de = interp2(Alpha, de, params.Aerodynamics.CDde, Alpha_, de_);
CD_dr = interp2(Beta, dr, params.Aerodynamics.CDdr, Alpha_, dr_);

%
Cl_ab = interp2(Beta, Alpha, params.Aerodynamics.Cl_alpha_beta, Beta_, Alpha_);
Cl_p = interp1(Alpha, params.Aerodynamics.Clp, Alpha_);
Cl_q = interp1(Alpha, params.Aerodynamics.Clq, Alpha_);
Cl_r = interp1(Alpha, params.Aerodynamics.Clr, Alpha_);
Cl_da = interp2(Alpha, da, params.Aerodynamics.Clda, Alpha_, da_);
Cl_de = interp2(Alpha, de, params.Aerodynamics.Clde, Alpha_, de_);
Cl_dr = interp2(Beta, dr, params.Aerodynamics.Cldr, Alpha_, dr_);

%
Cm_ab = interp2(Beta, Alpha, params.Aerodynamics.Cm_alpha_beta, Beta_, Alpha_);
Cm_p = interp1(Alpha, params.Aerodynamics.Cmp, Alpha_);
Cm_q = interp1(Alpha, params.Aerodynamics.Cmq, Alpha_);
Cm_r = interp1(Alpha, params.Aerodynamics.Cmr, Alpha_);
Cm_da = interp2(Alpha, da, params.Aerodynamics.Cmda, Alpha_, da_);
Cm_de = interp2(Alpha, de, params.Aerodynamics.Cmde, Alpha_, de_);
Cm_dr = interp2(Beta, dr, params.Aerodynamics.Cmdr, Alpha_, dr_);

%
Cn_ab = interp2(Beta, Alpha, params.Aerodynamics.Cn_alpha_beta, Beta_, Alpha_);
Cn_p = interp1(Alpha, params.Aerodynamics.Cnp, Alpha_);
Cn_q = interp1(Alpha, params.Aerodynamics.Cnq, Alpha_);
Cn_r = interp1(Alpha, params.Aerodynamics.Cnr, Alpha_);
Cn_da = interp2(Alpha, da, params.Aerodynamics.Cnda, Alpha_, da_);
Cn_de = interp2(Alpha, de, params.Aerodynamics.Cnde, Alpha_, de_);
Cn_dr = interp2(Beta, dr, params.Aerodynamics.Cndr, Alpha_, dr_);

%% GRAVITY
Gravity_Forces_body = DCM_e2b*[0; 0; mass*9.81];

%% PROPULSION

Propulsion_Forces_body = [1/2*Sprop*Cprop*rho_msl*((k_motor*dth_)^2 - Va^2); 0; 0];
Propulsion_Moments_body = [(k_omega*dth_)^2*(-k_Tp); 0; 0];


%% AERODYNAMIC FORCES

q_dyn = 1/2*rho_msl*Va^2;

Lift_Force = (CL_ab + CL_p*wing_span/(2*Va)*w_body(1) + ...
    CL_q*chord/(2*Va)*w_body(2) + CL_r*wing_span/(2*Va)*w_body(3) + ...
    CL_da + CL_de + CL_dr)*q_dyn*S;

Drag_Force = (CD_ab + CD_p*wing_span/(2*Va)*w_body(1) + ...
    CD_q*chord/(2*Va)*w_body(2) + CD_r*wing_span/(2*Va)*w_body(3) + ...
    CD_da + CD_de + CD_dr)*q_dyn*S;

Side_Force = (CY_ab + CY_p*wing_span/(2*Va)*w_body(1) + ...
    CY_q*chord/(2*Va)*w_body(2) + CY_r*wing_span/(2*Va)*w_body(3) + ...
    CY_da + CY_de + CY_dr)*q_dyn*S;


Aero_Roll_Moment = (Cl_ab + Cl_p*wing_span/(2*Va)*w_body(1) + ...
    Cl_q*chord/(2*Va)*w_body(2) + Cl_r*wing_span/(2*Va)*w_body(3) + ...
    Cl_da + Cl_de + Cl_dr)*wing_span*q_dyn*S;

Aero_Pitch_Moment = (Cm_ab + Cm_p*wing_span/(2*Va)*w_body(1) + ...
    Cm_q*chord/(2*Va)*w_body(2) + Cm_r*wing_span/(2*Va)*w_body(3) + ...
    Cm_da + Cm_de + Cm_dr)*chord*q_dyn*S;

Aero_Yaw_Moment = (Cn_ab + Cn_p*wing_span/(2*Va)*w_body(1) + ...
    Cn_q*chord/(2*Va)*w_body(2) + Cn_r*wing_span/(2*Va)*w_body(3) + ...
    Cn_da + Cn_de + Cn_dr)*wing_span*q_dyn*S;


Aero_Forces_body = DCM_b2w'*[-Drag_Force; 0; -Lift_Force] + [0; Side_Force; 0];
Aero_Moments_body = [Aero_Roll_Moment; Aero_Pitch_Moment; Aero_Yaw_Moment];

%% TOTAL FORCES AND MOMENTS

Total_Moments_body = Aero_Moments_body + Propulsion_Moments_body;
Total_Forces_body = Aero_Forces_body + Propulsion_Forces_body + Gravity_Forces_body;


%%

w_body_dot = pinv(Inertia_Tensor)*(Total_Moments_body - cross(w_body, Inertia_Tensor*w_body));
v_body_dot = 1/mass*(Total_Forces_body - cross(w_body, mass*v_body));
euler_rate = Rot_pqr2euler*w_body;

%%
X_dot(1:3) = v_body_dot;
X_dot(4:6) = w_body_dot;
X_dot(7:9) = euler_rate;

end

