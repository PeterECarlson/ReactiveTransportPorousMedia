function [A] = comp_mean(a,p,Grid)
% author: Peter Carlson
% date: 2/16/15
% Description:
% Takes coefficient field, a, defined at the cell centers and computes the
% mean specified by the power, p and returns it in a sparse diagonal
% matrix, A.
%
% Input:
% a = Nx by 1 column vector of cell centered values
% p = power of the generalized mean
% 1 (arithmetic mean)
% -1 (harmonic mean)
% Grid = structure containing information about the grid.
%
% Output
% A = Nx+1 by Nx+1 diagonal matrix of power means at the cell faces.
%
% Example call:
% a = @(x) 1+x.^3;
% Grid.xmin = 0; Grid.xmax = 1; Grid.Nx = 10;% Grid = build_grid(Grid);
% A = comp_mean(a(Grid.xc),1,Grid);

a_p = a.^p;
pm_op = full(spdiags(repmat([1,1],[Grid.Nx,1]),[0,1],Grid.Nx,Grid.Nx+1)/2)'; %Nx by Nx+1 grid with two diagonals of 0.5
pm_op(1,1) = 1;
pm_op(end,end) = 1;
A = diag((pm_op*a_p).^(1/p));