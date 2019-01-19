% pec543_ADEnumerical1D.m
% Author: Peter Carlson
% Date: 4/2/2015
% Solves advection of dissolved chemical constituents
% using the theta method with theta = 1, theta = 0, theta = 0.5
% Calculates numerical diffusion for theta = 1.



close all
clear all
clc

g = 981; %cm/s^2
rho = 1; %g/cm^3
mu = 1; %g/cm/s
Length = 45; %Length cm
r = 2.5; %inner radius of column cm
d = 0.023; %diameter of beads cm
phi = 0.4; %porosity
Q = 10/60; %cm^3/sec
hr = 0; %atmospheric pressure
Dm = 2e-4; %cm^2/sec
alphaL = 0.1; %cm
tau = sqrt(2);
CFL = 0.9;


K = (rho*g/mu)*(phi^3/(1-phi)^2)*(d^2/180)


Grid.xmin = 0; Grid.xmax = Length; Grid.Nx = 100;
Grid = build_grid(Grid);

qL_BC = Q/(pi*r^2);

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


%%Theta Method

cL_BC = 1;
cR_BC = 0;
Dh = Dm/tau + abs(alphaL*q(1)/phi);;
T = 15*60;
dt = abs(0.5*phi*Grid.dx./q(1));
[D,G,I] = build_ops(Grid);
qn = min(q(1:Grid.Nx),0);
qp = max(q(2:Grid.Nfx),0);
A = full(spdiags([qp,qn],[1,0],Grid.Nx,Grid.Nfx))';
Lc = D*(A-Dh*speye(Grid.Nfx)*G);
Bc = I([1,Grid.Nx],:);
gc = [cL_BC, cR_BC]';
Nc = spnull(Bc);
M = phi*I;


theta = 0;
Im = M + dt*(1-theta)*Lc;
Ex = M - dt*theta*Lc;
c = zeros(Grid.Nx,1);

for i = 1:1:ceil(T/dt)-1
    c = solve_LBVP(Im,Ex*c,Bc,gc,Nc);
end
c_reverse = c;


theta = 1;
Im = M + dt*(1-theta)*Lc;
Ex = M - dt*theta*Lc;
c = zeros(Grid.Nx,1);

for i = 1:1:ceil(T/dt)-1
    c = solve_LBVP(Im,Ex*c,Bc,gc,Nc);
end
c_forward = c;


theta = 0.5;
Im = M + dt*(1-theta)*Lc;
Ex = M - dt*theta*Lc;
c = zeros(Grid.Nx,1);

for i = 1:1:ceil(T/dt)-1
    c = solve_LBVP(Im,Ex*c,Bc,gc,Nc);
end
c_cranknicholson = c;

%%Analytic
x = linspace(0,1);
t = (15*60)/(phi*Length/(q(1)));
D_num = q(1)*Grid.dx/(2*phi)*(1-q(1)/phi*dt/Grid.dx)
Pe = Length*q(1)/phi/D_num
c_analytic = ADEanalytic(x, t, Pe);




figure
plot(Grid.xc, c_reverse,...
    Grid.xc, c_forward,...
    Grid.xc, c_cranknicholson,...
    x*Length, c_analytic)
xlabel('x')
ylabel('c')
legend('Reverse',...
    'Forward',...
    'Crank-Nicholson',...
    'Analytic Numerical Diffusion')