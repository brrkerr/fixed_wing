function [F0] = cost_straight_level(Z, params)

X = Z(1:9);
U = Z(10:13);

xdot = aircraft_model(X, U, params);
theta = X(8);
Va = sqrt(X(1)^2 + X(2)^2 + X(3)^2);
Alpha = atan2(X(3), X(1));
Gamma = theta - Alpha;

Q = [xdot;
    Va - 30;
    Gamma;
    X(2);
    X(7)];

H = diag(ones(1,13));

F0 = Q'*H*Q;


end

