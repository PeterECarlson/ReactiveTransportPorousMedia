% pec543_HeterogeneousColumnNumerical_convergence.m
% author: Peter Carlson
% date: 16 Feb 2015
% Description
% Solves flow through heterogenous aquifer using arithmetic mean methods
% Problem 3:4


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

Grid.xmin = 0; Grid.xmax = Length;
Nxs = [5 10 15 20];
h = zeros (Nxs(4),4);
q = zeros (Nxs(4)+1,4);
Gridxclist = zeros (20,4);
Gridxflist = zeros (21,4);

for i = 1:length(Nxs)
    Grid.Nx = Nxs(i);
    Grid = build_grid(Grid);
    hL_BC = h1_analytic(Grid.xc(1));
    hR_BC = h3_analytic(Grid.xc(Grid.Nx));
    [D,G,I] = build_ops(Grid);
    f = zeros(Grid.Nx,1);
    B = I([1, Grid.Nx],:);
    g = [hL_BC, hR_BC]';
    N = spnull(B);
    Kc = zeros(Grid.Nx,1);
    for j = 1:Grid.Nx
       if Grid.xc(j) < b1 
           Kc(j) = K1;
       
       elseif Grid.xc(j) <b2
           Kc(j) = K2;
       else
           Kc(j) = K1;
       end
    end    
    K = comp_mean(Kc,1,Grid); 
    L = -D*(K*G);
    h(1:Grid.Nx,i) = solve_LBVP(L,f,B,g,N);
    q(1:Grid.Nfx,i) = comp_flux(D,K,G,h(1:Grid.Nx,i),f,Grid);
    Gridxclist(1:Grid.Nx,i) = Grid.xc;
    Gridxflist(1:Grid.Nfx,i) = Grid.xf;
end    



x1 = 0:b1;
x2 = b1:b2;
x3 = b2:Length;


figure

subplot(2,1,1)
plot([x1 x2 x3], [h1_analytic(x1) h2_analytic(x2) h3_analytic(x3)],...
 Gridxclist(1:Nxs(1),1), h(1:Nxs(1),1),...
 Gridxclist(1:Nxs(2),2), h(1:Nxs(2),2),...
 Gridxclist(1:Nxs(3),3), h(1:Nxs(3),3),...
 Gridxclist(1:Nxs(4),4), h(1:Nxs(4),4))
title('Head')
xlabel('m')
ylabel('m')
legend('Analytic','Arithmetic 5','Arithmetic 10','Arithmetic 15',...
    'Arithmetic 20')

subplot(2,1,2)
plot([x1 x2 x3], q_analytic([x1 x2 x3]),...
    Gridxflist(1:Nxs(1)+1,1), q(1:Nxs(1)+1,1),...
    Gridxflist(1:Nxs(2)+1,2), q(1:Nxs(2)+1,2),...
    Gridxflist(1:Nxs(3)+1,3), q(1:Nxs(3)+1,3),...
    Gridxflist(1:Nxs(4)+1,4), q(1:Nxs(4)+1,4))
title('Flow')
xlabel('m')
ylabel('m/s')
legend('Analytic','Arithmetic 5','Arithmetic 10','Arithmetic 15',...
    'Arithmetic 20')