% plot_Cd_shock.m
% Author: Peter Carlson
% Date: 4/20/2015
% Analytic solution to the Reimann problem of Cd
% Problem 1.3



Length = 45 ; %cm
r = 2.5 ; %cm
Q = 10/60; %ml/s
cmax = 9.5; %ug/g
Kcd = 32; %cm^3/ug
c_r = 1/1000; %ug/cm^3
c_l = 40/1000; %ug/cm^3
phi = 0.4; %loamy sand.
v = Q/(pi*r^2);
rho_s = 2.65; %g/cm3


t = 50*3600;

Kd = (1-phi)/phi*cmax*Kcd;
c_s_r = cmax*Kcd*c_r/(1+Kcd*c_r);
c_s_l = cmax*Kcd*c_l/(1+Kcd*c_l);
dcsdc = (c_s_r-c_s_l)/(c_r-c_l);
R = 1 + (1-phi)/phi*rho_s*dcsdc;
v_c = v/R;
x_shock = (v_c * t);


c = linspace(c_l,c_r);


t_constant = 45/v_c/3600/24 %days

figure
plot([linspace(0,x_shock) linspace(x_shock,Length)],...
    [c_l*ones(1,100)*1000 c_r*ones(1,100)*1000])
xlabel('x, cm')
ylabel('c, \mug/L')

