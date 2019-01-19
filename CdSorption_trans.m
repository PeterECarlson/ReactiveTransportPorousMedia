% CdSorption_trans.m
% Author: Peter Carlson
% Date: 4/28/2015
% Driver script
% Question 2.3


tol = 1E-9;
imax = 10;
alphaL = 0.1; %cm
CFL = 100;
T = 15*3600;
theta = 0;

Dh = spdiags(Dm/tau + alphaL*q/phi, 0, Grid.Nfx, Grid.Nfx);



%Determine the time step and total time run
dt = abs(CFL*phi * Grid.dx/q(1));
Nt =ceil(T/dt);

%Calculate the explicit and implicit linear matricies
M = phi*I;
[Lim, Lex] = build_ade_ops(Grid,I,D,G,q,M,Dh,dt);

%Build the LVBP conditions
dof_dir = 1;
B = I(dof_dir,:);
g = 0;
N = I;
N(:,dof_dir) = [];

%Set up initial conditions
c = ones(Grid.Nx, 1)*c_r;
c(1) = c_l;
t = 0;
for n=2:Nt
    i=1;nres=1;ndc=1;
    cn = c;
    while (nres>tol || ndc>tol)&&i<imax
        r = residual_CdSorption(c, cn, Lim(theta), Lex(theta), cs, rhob);
        nres = norm(N'*r);
        J = jacobian_CdSorption(c, I, Lim(theta), dcsdc, rhob, Grid);
        dc = solve_LBVP(J, -r, B, g, N);
        ndc = norm(N'*dc);
        c = c + dc;
        i = i + 1;
        %fprintf('%d: nres = %3.2e, ndc = %3.2e;\n', i, nres, ndc);
    end
end



figure
plot(Grid.xc,c*1000)
xlabel('x [cm]')
ylabel('Cd [\mug/L]')
title('Numerical')