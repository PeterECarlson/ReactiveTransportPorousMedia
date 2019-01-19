% pec543_OperatorsPolarSpherical.m
% author: Peter Carlson
% date: 24 Feb 2015
% Description
% Solves LBVP in cartesian, polar, and spherical coordinates.
% Problem 2:5

close all
clear all
clc


K = 10E-5; %m/s
Qo = 10; %m^3/s
A = 1; %m^2
H = 1; %m
xo = 100; %m
xw = 0.1; %m
ho = 0; %set this so you can calculate dh

%%analytical

h_analytic_cartesian = @(x) -Qo./(K*A).*(x-xo) + ho;

h_analytic_polar = @(x) -Qo./(2*pi*H*K).*log(x/xo) + ho;

h_analytic_spherical = @(x) Qo./(4*pi*K).*(1./x-1/(xo)) + ho;



Grid.xmin = xw;
Grid.xmax = xo;
Grid.Nx = 50;



% Cartesian


Grid.geom = 1;
Grid = build_grid(Grid);
[D,G,I] = build_ops(Grid);
L = -D*(K*G);

hL_BC = h_analytic_cartesian(Grid.xc(1));
hR_BC = h_analytic_cartesian(Grid.xc(Grid.Nx));

f = zeros(Grid.Nx,1);
B = I([1, Grid.Nx],:);
g = [hL_BC, hR_BC]';
N = spnull(B);
h_cartesian = solve_LBVP(L,f,B,g,N);

%Polar
Grid.geom = 2;
Grid = build_grid(Grid);
[D,G,I] = build_ops(Grid);
L = -D*(K*G);

hL_BC = h_analytic_polar(Grid.xc(1));
hR_BC = h_analytic_polar(Grid.xc(Grid.Nx));

f = zeros(Grid.Nx,1);
B = I([1, Grid.Nx],:);
g = [hL_BC, hR_BC]';
N = spnull(B);
h_polar = solve_LBVP(L,f,B,g,N);

%Spherical

Grid.geom = 3;
Grid = build_grid(Grid);
[D,G,I] = build_ops(Grid);
L = -D*(K*G);

hL_BC = h_analytic_spherical(Grid.xc(1));
hR_BC = h_analytic_spherical(Grid.xc(Grid.Nx));

f = zeros(Grid.Nx,1);
B = I([1, Grid.Nx],:);
g = [hL_BC, hR_BC]';
N = spnull(B);
h_spherical = solve_LBVP(L,f,B,g,N);



figure

subplot(3,1,1)
plot(Grid.xc, h_analytic_cartesian(Grid.xc),...
    Grid.xc, h_cartesian);
title('Cartesian')
xlabel('m')
ylabel('m')
legend('Analytic', 'Numeric')

subplot(3,1,2)
plot(Grid.xc,h_analytic_polar(Grid.xc),...
Grid.xc,h_polar);
title('Polar')
xlabel('m')
ylabel('m')

subplot(3,1,3)
plot(Grid.xc, h_analytic_spherical(Grid.xc),...
    Grid.xc, h_spherical);
title('Spherical')
xlabel('m')
ylabel('m')