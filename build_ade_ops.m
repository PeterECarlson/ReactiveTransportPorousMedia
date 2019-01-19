function [Lim,Lex] = build_ade_ops(Grid,I,D,G,q,phi,Disp,dt)
% file: build_ade_ops.m
% author: Peter Carlson
% date: 4/7/2015
% Description:
% Assembles the function handles for the implicit and explicit matrices
% arising from the discretization of the advection-dispersion equation.
% The hydrodynamic dispersion term is always treated implicitly. The
% advection term is treated using the theta-method:
% theta = 1 is the explicit Forward Euler method,
% theta = 0 is the implicit Backward Euler method.
% Input:
% Grid = grid data structure.
% I = Nc by Nc identity matrix.
% D = Nc by Nf divergence matrix.
% G = Nf by Nc gradient matrix.
% q = Nf by 1 flux vector
% phi = Nc by 1 vector of porosities.
% Disp = Hydrodynamic dispersion.
% dt = time step
% Output:
% Lim = @(theta) Implicit N by N matrix from the discretization of the
% advection-diffusion equation. Lim should be an anonymous function
% of theta.
% Lex = @(theta) Explicit N by N matrix from the discretization of the
% advection-diffusion equation. Lex should be an anonymous function
% of theta.


qn = min(q(1:Grid.Nx),0);
qp = max(q(2:Grid.Nfx),0);
A = full(spdiags([qp,qn],[1,0],Grid.Nx,Grid.Nfx))';
Lc = D*(A-Disp*speye(Grid.Nfx)*G);



Lim = @(theta) phi + dt*(1-theta)*Lc;
Lex = @(theta) phi - dt*theta*Lc;
end