clear all
clc
close all

initialization = 1;

if (initialization == 0)
    Z_guess = zeros(13,1);
    Z_guess(1) = 30;
else
    load trim_straight_level
    Z_guess = [XStar;UStar];
end

params = load('aircraft_params.mat');

[ZStar, f0] = fminsearch(@(Z)cost_straight_level(Z, params), Z_guess, ...
    optimset('TolX', 1e-10, 'MaxFunEvals', 10000, 'MaxIter', 10000));


XStar = ZStar(1:9)
UStar = ZStar(10:13)
f0

XdotStar = aircraft_model(XStar, UStar, params)
VaStar = sqrt(XStar(1)^2 + XStar(2)^2 + XStar(3)^2)
AlphaStar = rad2deg(atan2(XStar(3), XStar(1)))
GammaStar = rad2deg(XStar(8)) - AlphaStar
vStar = XStar(2)
phiStar = rad2deg(XStar(7))
psiStar = rad2deg(XStar(9))

%%
save trim_straight_level ...
    XStar UStar
