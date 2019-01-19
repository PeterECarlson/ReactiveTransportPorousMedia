% plot_Cdisotherm.m
% Author: Peter Carlson
% Date: 4/20/2015
% Reads Cd_data.mat and plots
%Question 1.1

close all
clear all
clc

load('Cd_isotherm_data.mat')

cmax = 9.5; %ug/g
Kcd = 3.2e-2; %L/ug
csCd_lang = cmax*Kcd.*cCd./(1+Kcd.*cCd);

cCd_nondim = cCd/Kcd/1000;
csCd_nondim = csCd/cmax;
csCd_lang_nondim = cCd_nondim./(1+cCd_nondim);

figure
plot(cCd, csCd,...
    cCd, csCd_lang)
xlabel('c, \mug/L')
ylabel('c_s, \mug/g')


figure
plot(cCd_nondim, csCd_nondim,...
    cCd_nondim, csCd_lang_nondim)
xlabel('c, []')
ylabel('c_s, []')