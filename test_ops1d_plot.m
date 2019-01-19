% author: Peter Carlson
% date: 2/3/2015
% Description:
% This tests the divergence and gradient matrices by performing a forward
% test on the double derivative of e^cos(2pix)
close all
clear all
clc


Grid.xmin = 0; Grid.xmax = 3; Grid.Nx = 30;
Grid = build_grid(Grid)
[D,G,I]=build_ops(Grid);

f = @(x) exp(cos(2*pi*x));
ddfdxx = @(x) 4*pi^2*exp(cos(2*pi*x)).*(sin(2*pi*x).*sin(2*pi*x)-cos(2*pi*x));


L = D*G;
ldisc = L*f(Grid.xc');

figure 
plot(Grid.xc,ddfdxx(Grid.xc),'r',Grid.xc,ldisc,'go')
