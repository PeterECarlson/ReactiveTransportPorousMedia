% plot_Cd_rarefaction.m
% Author: Peter Carlson
% Date: 4/20/2015
% Analytic solution to the Reimann problem of Cd
% Problem 1.2



Length = 45 ; %cm
r = 2.5 ; %cm
Q = 10/60; %ml/s
cmax = 9.5; %ug/g
Kcd = 32; %cm^3/ug
c_r = 40/1000; %ug/cm^3
c_l = 1/1000; %ug/cm^3
phi = 0.4; %loamy sand.
v = Q/(pi*r^2);
rho_s = 2.65; %g/cm3


t = 50*3600; 

Kd = (1-phi)/phi*cmax*Kcd;

dcsdc = @(c) cmax*Kcd./(1+(Kcd*c)).^2;
R = @(c) 1 + (1-phi)/phi*rho_s*dcsdc(c);
v_c = @(c) v./R(c);
c_rare = linspace(c_l, c_r);
x_rare = (v_c(c_rare)*t);

v_l = v./R(c_l);
v_r = v./R(c_r);

x_l = v_l*t;
x_r = v_r*t;


t_constant = 45/v_l/3600/24 %days

x = [linspace(0,x_l) x_rare linspace(x_r,Length)];
c = [c_l*ones(1,100) c_rare c_r*ones(1,100)];
figure
plot(x,c*1000) 
xlabel('x, cm')
ylabel('c, \mug/L')
title('Analytical')
