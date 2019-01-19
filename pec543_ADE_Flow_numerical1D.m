% pec543_Flow_numerical1D.m
% Author: Peter Carlson
% Date: 4/7/2015
% Solves advection of dissolved chemical constituents
% using the theta method with theta = 1, theta = 0, theta = 0.5
% Calculates numerical diffusion for theta = 1.



close all
clear all
clc

g = 981; %cm/s^2
rho = 1; %g/cm^3
mu = 0.01; %g/cm-s
Length = 45; %Length cm
r = 2.5; %inner radius of column cm
d = 0.023; %diameter of beads cm
phi = 0.4; %porosity
Q = 10/60; %cm^3/sec


hr = 0; %atmospheric pressure

K = (rho*g/mu)*(phi^3/(1-phi)^2)*(d^2/180)


Grid.xmin = 0; Grid.xmax = Length; Grid.Nx = 200;
Grid = build_grid(Grid);

qL_BC = Q/(pi*r^2)

h_analytic = @(x) -qL_BC * (x-Length) / K;
q_analytic = ones(Grid.Nfx,1)*qL_BC;
Q_analytic = ones(Grid.Nfx,1)*Q;

hR_BC = h_analytic(Grid.xc(Grid.Nx));

[D,G,I] = build_ops(Grid);
L = -D*(K*G);
f = zeros(Grid.Nx,1);
f(1) = f(1) + qL_BC/Grid.dx;
B = I([Grid.Nx],:);
g = [hR_BC]';
N = spnull(B);
h = solve_LBVP(L,f,B,g,N);

q = comp_flux(D,K,G,h,f,Grid);
q(1) = qL_BC;

figure

subplot(2,1,1)
plot(Grid.xc, h_analytic(Grid.xc),...
    Grid.xc, h)
title('Head')
xlabel('cm')
ylabel('cm')
legend('analytic','numeric')

subplot(2,1,2)
plot(Grid.xf, q_analytic,...
    Grid.xf, q)
title('q')
xlabel('cm')
ylabel('cm/s')
legend('analytic','numeric')