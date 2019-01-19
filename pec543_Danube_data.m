% pec543_Danube_data.m
% author: Peter Carlson
% date: 11 Feb 2015
% Description
% Tests function solve_LBVP.m
% Tests against analytical solution of Danube River problem.

%problem 1:

close all
clear all
clc

Danube_properties


Grid.xmin = 0; Grid.xmax = Length; Grid.Nx = 50;
Grid = build_grid(Grid);


topo_ic = interp1(dist_topo, topo, Grid.xc);
total_rain = qp/b*(Length-Grid.dx);
scaling = total_rain/trapz(Grid.xc,topo_ic);
rain_topo = scaling*topo_ic;

h_analytic = @(x) -(qp/(2*K*b))*x.^2 + ((hR-hL)/Length + ...
    Length*qp/(2*K*b))*x + hL;
hL_BC = h_analytic(Grid.xc(1));
hR_BC = h_analytic(Grid.xc(Grid.Nx));
[D,G,I] = build_ops(Grid);
L = -D*(K*G);
f_const = ones(Grid.Nx,1)*qp/b;
f_data = rain_topo';
B = I([1,Grid.Nx],:);
g = [hL_BC,hR_BC]';
N = spnull(B);
h_const = solve_LBVP(L,f_const,B,g,N);

h_data = solve_LBVP(L,f_data,B,g,N);

%The distributions of the rainfall do not strongly affect the final
%solutions.


figure
subplot(2,1,1)
plot(Grid.xc, f_const, Grid.xc, f_data)
legend('uniform rainfall','elevation-based rainfall')
xlabel('meters')
ylabel('meters')
subplot(2,1,2)
plot(Grid.xf, h_analytic(Grid.xf), Grid.xc, h_data, dist_gw, gw)
legend('analytic','numeric','actual')
xlabel('meters')
ylabel('meters')