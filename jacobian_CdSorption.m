function [J] = jacobian_CdSorption(c,I,Lim,dcsdc,rhob,Grid)
% author: Peter Carlson
% data: 28 Apr 2015
% Input:
% c = current iterate of concentration at n+1
% I = Identity matrix (sparse)
% Lim = Implicit advection diffusion matrix
% dcsdc = Derivative of the isotherm
% rhob = bulk density of porous medium
% Output:
% J = Jacobian matrix evluated at c

dCsdC = spdiags(dcsdc(c),0,Grid.Nx,Grid.Nx);

J = rhob*dCsdC + Lim;
end