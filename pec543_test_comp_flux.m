% pec543_test_comp_flux.m
% author: Peter Carlson
% date: 2/11/2015
% Description:
% Tests comp_flux.m with the Danube/Tisza river systems
% Problem 2:2

close all
clear all
clc

Danube_properties


Grid.xmin = 0; Grid.xmax = Length; Grid.Nx = 100;
Grid = build_grid(Grid);

h_analytic = @(x) -(qp/(2*K*b))*x.^2 + ((hR-hL)/Length + ...
    Length*qp/(2*K*b))*x + hL;
q_analytic = @(x) x*qp/b - K*(hR-hL)/Length - Length*qp/(2*b);
q_analytic(Grid.xf);

hL_BC = h_analytic(Grid.xc(1));
hR_BC = h_analytic(Grid.xc(Grid.Nx));
[D,G,I] = build_ops(Grid);
L = -D*(K*G);
f = ones(Grid.Nx,1)*qp/b;
B = I([1,Grid.Nx],:);
g = [hL_BC,hR_BC]';
N = spnull(B);
h = solve_LBVP(L,f,B,g,N);

q = comp_flux(D,K,G,h,f,Grid);

figure

subplot(2,1,1)
title('Groundwater heads')
plot(Grid.xf, h_analytic(Grid.xf), Grid.xc, h, dist_gw, gw)
legend('analytic','numeric','actual')
xlabel('meters')
ylabel('meters')

subplot(2,1,2)
title('Fluxes')
plot(Grid.xf, q_analytic(Grid.xf),Grid.xf,q)
legend('analytic', 'numeric')
xlabel('meters')
ylabel('meters/second')