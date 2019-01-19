function [u] = solve_LBVP(L,f,B,g,N)
% author: Peter Carlson
% date: 10 Feb 2015
% Description
% Computes the solution $u$ to the linear differential problem given by
%
% $$\mathcal{L}(u)=f \quad x\in\Omega $$
%
% with boundary conditions
%
% $$\mathcal{B}(u)=g \quad x\in\partial\Omega$$.
%
% Input:
% L = matrix representing the discretized linear operator of size N by N,
% where N is the number of degrees of freedom
% f = column vector representing the discretized r.h.s. and contributions
% due non-homogeneous Neumann BC’s of size N by 1
% B = matrix representing the constraints arising from Dirichlet BC’s of
% size Nc by N
% g = column vector representing the non-homogeneous Dirichlet BC’s of size
% Nc by 1.
% N = matrix representing a orthonormal basis for the null-space of B and
% of size N by (N-Nc).
% Output:
% u = column vector of the solution of size N by 1
%
% Example call:
% >> Grid.xmin = 0; Grid.xmax = 1; Grid.Nx = 10;
% >> Grid = build_grid(Grid);
% >> [D,G,I] = build_ops(Grid);
% >> L = -D*G; r = ones(Grid.Nx,1);
% >> B = I(1,:); g = 1; N = spnull(B);
% >> h = solve_LBVP(L,r,B,g,N);

h_p = B'*(B*B'\g);
h_o = N*((N'*L*N)\N'*(f-L*h_p));

u = h_p+h_o;