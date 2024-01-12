function [c,ceq] = nonlinear_const(Z)
X = Z(1:9);
U = Z(10:13);

xdot = aircraft_model(X, U);
theta = X(8);
Va = sqrt(X(1)^2 + X(2)^2 + X(3)^2);
Alpha = atan2(X(3), X(1));
Gamma = theta - Alpha;

ceq = [xdot;Gamma; Va-30];
c = [];

end

