function [D,G,I]=build_ops(Grid)
% author: Peter Carlson
% date: 2/3/2015
% Description:
% This function computes the discrete divergence and gradient matrices on a
% regular staggered grid using central difference approximations. The
% discrete gradient assumes homogeneous boundary conditions.
% Input:
% Grid = structure containing all pertinent information about the grid.
% Output:
% D = discrete divergence matrix
% G = discrete gradient matrix
% I = identity matrix
% Example call:
% >> Grid.xmin = 0; Grid.xmax = 1; Grid.Nx = 5;
% >> Grid = build_grid(Grid);
% >> [D,G,I]=build_ops(Grid);

 Xcinv = diag(Grid.xc.^(1-Grid.geom));
 Xf = diag(Grid.xf.^(Grid.geom-1));
 
 D = full(spdiags(repmat([-1,1],[Grid.Nx,1]),[0,1],Grid.Nx,Grid.Nx+1)/Grid.dx);
 G = full(-D');
 G(1,1)=0;
 G(Grid.Nfx,Grid.Nx) = 0;
 D = Xcinv*D*Xf;
 I = speye(Grid.Nx);
end