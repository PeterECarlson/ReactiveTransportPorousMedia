close all
clear all
clc

%field
L_f = 50000;
DH_f = 10;
k = 10E-12;
phi = 0.2;
tau = 1/sqrt(2);
mu = 10E-3;
rho = 10E3;
D_m = 10E-9;
r = 0.01;
g = 9.81;
K = k*rho*g/mu
L_l = 0.6;
alphaL_f = 10E2;
alphaL_l = 10E-3;

q_f = K*DH_f/L_f

Dh_f = tau*phi*D_m + alphaL_f*q_f

%Field Pe
Pe_f = L_f*q_f/Dh_f

A = pi*r^2
Q_max = 1000000000000000000000;

%Lab Pe
Pe_l = @(Q) L_l*Q/(A*tau*phi*D_m + alphaL_l*Q); 

Pe_l_max = Pe_l(Q_max)

Q_l = A*tau*phi*D_m/(L_l/Pe_f-alphaL_l)

t_s = L_l*A/Q_l  %sec
t_h = t_s/3600 %hour
t_d = t_h/24 %day