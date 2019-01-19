% author: Peter Carlson
% date: 2/3/2015
% Description:
% This computes the error of the discrete gradient and Laplacian

close all
clear all
clc

Grid.xmin = 0; Grid.xmax = 3;
fx = @(x) exp(cos(2*pi*x));
dfxdx = @(x) -2*pi*exp(cos(2*pi*x)).*sin(2*pi*x);
ddfxdxx= @(x) 4*pi^2*exp(cos(2*pi*x)).*(sin(2*pi*x).*sin(2*pi*x)-cos(2*pi*x));

EG = zeros(1,11);
EL = zeros(1,11);
n = [10 20 30 50 70 100 200 300 500 700 1000];

for i=1:1:11
    Grid.Nx = n(i);
    Grid = build_grid(Grid);
    [D,G,I]=build_ops(Grid);
    
    f=fx(Grid.xc)';
    full(G);
    gdisc = G*f;
    fprime=dfxdx(Grid.xf)';
    EG(i) = norm(fprime-gdisc)./norm(fprime);
    
    L = D*G;
    ldisc = L*f;
    fdprime=ddfxdxx(Grid.xc)';
    EL(i) = norm(fdprime-ldisc)./norm(fdprime);
end    
 

figure 
loglog(n,EG,n,EL)
xlabel('Nx')
ylabel('Error')
legend('Gradient Error','Lapacian Error')