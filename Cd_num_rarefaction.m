% CdSorption_numerical.m
% Author: Peter Carlson
% Date: 4/28/2015
% Driver script
% Question 2

close all
clear all
clc

Length = 45 ; %cm
r = 2.5 ; %cm
Q = 10/60; %ml/s
cmax = 9.5; %ug/g
Kcd = 32; %cm^3/ug
c_r = 40/1000; %ug/cm^3
c_l = 1/1000; %ug/cm^3
phi = 0.4; %loamy sand.
v = Q/(pi*r^2);
rhos = 2.65; %g/cm3;
rhob = phi*rhos+(1-phi);
tau = 2^.5;
d = 0.023;
Dm = 2E-4;
hr = 0;
g = 9.81;
mu = 8.9E-4;

Nx = 1000;
CFL = 100;

cs = @(c) cmax*Kcd*c./(1+Kcd*c);
dcsdc = @(c) cmax*Kcd./(1+Kcd*c).^2;

Cd_Sorption_flow
CdSorption_trans
plot_Cd_rarefaction
