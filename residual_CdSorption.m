function [r] = residual_CdSorption(c,cn,Lim,Lex,cs,rhob)
% author: Peter Carlson
% data: 28 Apr 2015
% Input:
% c = current iterate of concentration at n+1
% cn = c at n
% I = Identity matrix (sparse)
% Lim = Implicit advection-diffusion matrix
% Lex = Explicit advection-diffusion matrix
% cs = Anonymous function for isotherm
% rhob = bulk density of porous medium
% Output:
% r = residual, N by 1 vector evaluated at c and cn

r = rhob*cs(c)+Lim*c-rhob*cs(cn)-Lex*cn;