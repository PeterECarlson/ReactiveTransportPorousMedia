% pec543_Danube_discharge.m
% author: Peter Carlson
% date: 2/11/2015
% Description:
% Calculates total discharge into both rivers using different methods.
%Problem 2:3


pec543_test_comp_flux


Q_numeric = (-q(1) + q(Grid.Nfx))*Width*b*yr2s %m/s*m*m*s/yr = m^3/yr
Q_precip = qp*Length*Width*yr2s %m/s*m*m*s/yr = m^3/yr

%Total mass is conserved, because Q_numeric is equal to the total amount
%of precipitation

gradL = (gw(2)-gw(1))/(dist_gw(2)-dist_gw(1));
gradR = (gw(end)-gw(end-1))/(dist_gw(end)-dist_gw(end-1));
Q_actual = K*(gradL-gradR)*Width*b*yr2s

%Q_actual is roughly 50% of the modelled versions. I would adjust the
%model by accounting for evapotranspiration.