function [q] = comp_flux(D,Kd,G,h,fs,Grid)

% comp_flux.m
% author: Peter Carlson
% date: 2/11/2015
% Description:
% Computes the mass conservative fluxes across all boundaries from the
% residual of the compatability condition over the boundary cells.
% Input:
% D = Nx by Nx+1 discrete divergence matrix.
% Kd = Nx+1 by Nx+1 conductivity matrix or a scalar K.
% G = Nx+1 by Nx discrete gradient matrix.
% h = Nx by 1 vector of heads in cell centers.
% fs = N by 1 right hand side vector containing only source terms.
% Grid = structure containing grid information.
% Output:
% q = Nx+1 by 1 vector of fluxes across cell faces
% Example call:
% >> Grid.xmin = 0; Grid.xmax = 1; Grid.Nx = 10;
% >> Grid = build_grid(Grid);
% >> [D,G,I] = build_ops(Grid);
% >> L = -D*G; fs = ones(Grid.Nx,1);
% >> dof_dir = 1;
% >> B = I(dof_dir,:); g = 1;
% >> N = I; N(:,dof_dir) = [];
% >> h = solve_lbvp(L,fs,B,g,N);
% >> q = comp_flux(D,1,G,h,fs,Grid);

%Problem 2:1

q = -Kd*G*h;
q(1) = Grid.dx*(D(1,:)*q - fs(1));
q(Grid.Nfx) = -1*Grid.dx*(D(Grid.Nx,:)*q - fs(Grid.Nx));
