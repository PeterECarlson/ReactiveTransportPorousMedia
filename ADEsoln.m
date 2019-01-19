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

close all
clear all
clc

x = linspace(0,1);
t = [0.1, 0.3, 0.5, 0.7, 0.9];
Pe = 1000;
c = ADEanalytic(x, t, Pe);

figure
plot(x, c(1,:), x, c(2,:), x, c(3,:), x, c(4,:), x, c(5,:))
xlabel('x')
ylabel('c')
legend('t=0.1', 't=0.3', 't=0.5', 't=0.7', 't=0.9')