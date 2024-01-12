function [E, A_p, B_p] = FindLinearizedModel(XDOTo, Xo, Uo, DXDOT, DX, DU, params)

n = length(XDOTo);
m = length(Uo);


%%
E = zeros(n,n);

for i = 1:n
    for j = 1:n

        dxdot = DXDOT(i,j);

        xdot_plus = XDOTo;
        xdot_minus = XDOTo;

        xdot_plus(j) = xdot_plus(j) + dxdot;
        xdot_minus(j) = xdot_minus(j) - dxdot;

        F = aircraft_model(Xo, Uo, params) - xdot_plus;
        F_plus = F(i);

        F = aircraft_model(Xo, Uo, params) - xdot_minus;
        F_minus = F(i);

        E(i,j) = (F_plus - F_minus)/(2*dxdot);
    end
end

%%
A_p = zeros(n,n);

for i = 1:n
    for j = 1:n

        dx = DX(i,j);

        x_plus = Xo;
        x_minus = Xo;

        x_plus(j) = x_plus(j) + dx;
        x_minus(j) = x_minus(j) - dx;

        F = aircraft_model(x_plus, Uo, params) - XDOTo;
        F_plus = F(i);

        F = aircraft_model(x_minus, Uo, params) - XDOTo;
        F_minus = F(i);

        A_p(i,j) = (F_plus - F_minus)/(2*dx);
    end
end

%%
B_p = zeros(n,m);

for i = 1:n
    for j = 1:m

        du = DU(i,j);

        u_plus = Uo;
        u_minus = Uo;

        u_plus(j) = u_plus(j) + du;
        u_minus(j) = u_minus(j) - du;

        F = aircraft_model(Xo, u_plus, params) - XDOTo;
        F_plus = F(i);

        F = aircraft_model(Xo, u_minus, params) - XDOTo;
        F_minus = F(i);

        B_p(i,j) = (F_plus - F_minus)/(2*du);
    end
end


end

