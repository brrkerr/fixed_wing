clear all;
close all;
clc;

params = load('aircraft_params.mat');
load trim_straight_level.mat

Xdoto = aircraft_model(XStar, UStar, params);
Xo = XStar;
Uo = UStar;

dxdot_matrix = 10e-12*ones(9,9);
dx_matrix = 10e-12*ones(9,9);
du_matrix = 10e-12*ones(9,5);

[E, Ap, Bp] = FindLinearizedModel(Xdoto, Xo, Uo, ...
    dxdot_matrix, dx_matrix, du_matrix, params);

A = -pinv(E)*Ap;
B = -pinv(E)*Bp;

%% [Long; Lat; Psi]
Tinv = zeros(size(A));
Tinv(1,1) = 1;
Tinv(2,3) = 1;
Tinv(3,5) = 1;
Tinv(4,8) = 1;
Tinv(5,2) = 1;
Tinv(6,4) = 1;
Tinv(7,6) = 1;
Tinv(8,7) = 1;
Tinv(9,9) = 1;

A_long_lat_psi = Tinv*A*pinv(Tinv);
B_long_lat_psi = Tinv*B;

%% Longitudinal (u,w,q,theta) -> (x1,x3,x5,x8)

A_long = A_long_lat_psi(1:4,1:4);
B_long = B_long_lat_psi(1:4,:);

C = eye(size(A_long));
D = zeros(size(B_long));

sys_long = ss(A_long, B_long, C, D);

[num_long, den_long] = ss2tf(A_long, B_long, C, D, 2);

tf_q_de = tf(num_long(3,:), den_long);

%% Lateral (v,p,r,phi) -> (x2,x4,x6,x7)

A_lat = A_long_lat_psi(5:8,5:8);
B_lat = B_long_lat_psi(5:8,:);

sys_lat = ss(A_lat, B_lat, C, D);

[num_lat, den_lat] = ss2tf(A_lat, B_lat, C, D, 1);

tf_p_da = tf(num_lat(2,:), den_lat);

%% TF_servo
tf_servo = tf([10], [1 10]);


save linear_models ...
    A B C D ...
    A_long_lat_psi B_long_lat_psi ...
    A_long B_long A_lat B_lat ...
    tf_q_de tf_p_da tf_servo
