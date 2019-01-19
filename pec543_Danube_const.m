% pec543_Danube_const.m
% author: Peter Carlson
% date: 10 Feb 2015
% Description
% Tests function solve_LBVP.m
% Tests against analytical solution of Danube River problem.

%Problem 1:3

close all
clear all
clc

Danube_properties




Grid.xmin = 0; Grid.xmax = Length; Grid.Nx = 50;
Grid = build_grid(Grid);

h_analytic = @(x) -(qp/(2*K*b))*x.^2 + ...
    ((hR-hL)/Length + Length*qp/(2*K*b))*x + hL;
hL_BC = h_analytic(Grid.xc(1));
hR_BC = h_analytic(Grid.xc(Grid.Nx));
[D,G,I] = build_ops(Grid);
L = -D*(K*G);
f = ones(Grid.Nx,1)*qp/b;
B = I([1,Grid.Nx],:);
g = [hL_BC,hR_BC]';
N = spnull(B);
h = solve_LBVP(L,f,B,g,N);

figure
title('Head')
plot(Grid.xf, h_analytic(Grid.xf), Grid.xc, h,...
    dist_gw, gw, dist_topo, topo)
xlabel('meters')
ylabel('meters')
legend('analytic','numeric','actual','topography')