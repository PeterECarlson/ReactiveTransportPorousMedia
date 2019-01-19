% pec543_NeumannBC.m
% author: Peter Carlson
% date: 16 Feb 2015
% Description
% Solves Laboratory Core Flood through Column Problem

%Problem 2:3

close all
clear all
clc

Length = 45; % cm
r = 2.5; % cm
Q = 10; % mL/ min
phi = 0.4; % []
d = 0.023; % cm
ro = 1; % g/cm^2
g = -981; % cm/s^2 
mu = 1; % g/(cm*s)
K = ro*g/mu*(phi)^3/(1-phi)^2*d^2/180; 

Grid.xmin = 0; Grid.xmax = Length; Grid.Nx = 10;
Grid = build_grid(Grid);

qL_BC = Q/(pi*r^2);

h_analytic = @(x) -qL_BC * x / K;
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
Q = q*pi*r^2;

figure

subplot(3,1,1)
plot(Grid.xc, h_analytic(Grid.xc), Grid.xc, h)
title('Head')
xlabel('cm')
ylabel('cm')
legend('analytic','numeric')

subplot(3,1,2)
plot(Grid.xf, q_analytic, Grid.xf, q)
title('q')
xlabel('cm/s')
ylabel('cm')
legend('analytic','numeric')

subplot(3,1,3)
plot(Grid.xf, Q_analytic, Grid.xf, Q)
title('Q')
xlabel('mL/s')
ylabel('cm')
legend('analytic','numeric')