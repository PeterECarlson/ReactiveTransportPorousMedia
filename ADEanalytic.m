function [c] = ADEanalytic(x, t, Pe)
%Author: Peter Carlson
%Date: 03/31/2015
%Description:
%   This evaluates the analytic solution of the advection-diffusion
%   equation.
%   t = nondimensional time
%   x = nondimensional length
%   Pe = Peclay number
%   c = nondimensional concentration
%Output:

xarray = repmat(x, length(t), 1);
tarray = repmat(t', 1, length(x));
c = 1/2*erfc(((Pe./(4*tarray)).^0.5).*(xarray - tarray)) + exp(...
    log(1/2*erfc(((Pe./(4*tarray)).^0.5).*(xarray + tarray))).*Pe.*xarray);
c(:,1) = 1/2*erfc(((Pe./(4*t')).^0.5).*(zeros(length(t),1) - t')) + ...
    1/2*erfc(((Pe./(4*t')).^0.5).*(zeros(length(t),1) + t'));


