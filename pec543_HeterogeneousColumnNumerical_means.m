% pec543_HeterogeneousColumnNumerical_means.m
% author: Peter Carlson
% date: 16 Feb 2015
% Description
% Solves flow through heterogenous aquifer using arithmetic and harmonic
% averaging methods
% Problem 3:3


close all
clear all
clc

Length = 100; % m
b1 = 20; % m
b2 = 40; % m

K1 = 5E-3; % m/s
K2 = 5E-5; % m/s

hL = 120; % m
hR = 100; % m

%%Analytical

C1 = (hR-hL)/((b1-b2)*(1-K1/K2)+Length);

h1_analytic = @(x) C1*x + hL;

h2_analytic = @(x) C1*(K1/K2*x + (1-K1/K2)*b1) + hL;

h3_analytic = @(x) C1*(x - Length) + hR;

q_analytic = @(x) -C1*K1;


Grid.xmin = 0; Grid.xmax = Length; Grid.Nx = 5;
Grid = build_grid(Grid);

hL_BC = h1_analytic(Grid.xc(1));
hR_BC = h3_analytic(Grid.xc(Grid.Nx));

% harmonic
K_harm = comp_mean([K1,K2,K1,K1,K1]',-1,Grid); 

[D,G,I] = build_ops(Grid);
L_harm = -D*(K_harm*G);
f = zeros(Grid.Nx,1);
B = I([1, Grid.Nx],:);
g = [hL_BC, hR_BC]';
N = spnull(B);
h_harm = solve_LBVP(L_harm,f,B,g,N);
q_harm = comp_flux(D,K_harm,G,h_harm,f,Grid);

%arithmetic

K_arith = comp_mean([K1,K2,K1,K1,K1]',1,Grid); 

L_arith = -D*(K_arith*G);
h_arith = solve_LBVP(L_arith,f,B,g,N);
q_arith = comp_flux(D,K_arith,G,h_arith,f,Grid);


x1 = 0:b1;
x2 = b1:b2;
x3 = b2:Length;


figure

subplot(2,1,1)
plot([x1 x2 x3], [h1_analytic(x1) h2_analytic(x2) h3_analytic(x3)],...
    Grid.xc,h_harm,...
    Grid.xc, h_arith)
title('Head')
xlabel('m')
ylabel('m')
legend('Analytic','Harmonic', 'Arithmetic')

subplot(2,1,2)
plot([x1 x2 x3], q_analytic([x1 x2 x3]),...
    Grid.xf, q_harm,...
    Grid.xf, q_arith);
title('Flow')
xlabel('m')
ylabel('m/s')
legend('Analytic','Harmonic', 'Arithmetic')