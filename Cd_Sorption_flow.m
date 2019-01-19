% Cd_Sorption_flow.m
% Author: Peter Carlson
% Date: 4/28/2015
% Solves advection of dissolved chemical constituents
% using the theta method with theta = 1, theta = 0, theta = 0.5
% Calculates numerical diffusion for theta = 1.



K = abs((rhob*g/mu)*(phi^3/(1-phi)^2)*(d^2/180))



Grid.xmin = 0; Grid.xmax = Length; Grid.Nx = Nx;
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