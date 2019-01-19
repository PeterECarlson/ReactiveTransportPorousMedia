% plot_Cdisotherm.m
% Author: Peter Carlson
% Date: 4/20/2015
% Reads Cd_data.mat and plots

close all
clear all
clc

F = fopen('Cd_data.mat','r')




%{
figure
plot(Grid.xc, c_reverse,...
    Grid.xc, c_forward,...
    Grid.xc, c_cranknicholson,...
    x*Length, c_analytic)
xlabel('x')
ylabel('c')
legend('Reverse',...    
    'Forward',...
    'Crank-Nicholson',...
    'Analytic Numerical Diffusion')

%}