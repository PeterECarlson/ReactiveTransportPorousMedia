%% Physical properties
cm2m = 1/100;        % cm to m conversion
yr2s = 365*24*60^2;  % yr to s conversion

Length = 85070;      % Distance between Danube and Tisza rivers [m]
Width = 5430;        % Width of segment considered [m]
K =2e-2*cm2m;        % Hydraulic conductivity [m/s]
qp = 1.5*cm2m/yr2s;  % Average annual precipitation [m3/m2/s]
hL = 90;             % Elevation of Danube river[m]
hR = 80;             % Elevation of Tisza river [m]
b = 100;             % Aquifer thickness [m]

% Data
mm2m = 1.81e3; % mm on map to m in reality
dist_topo = [ 0 14  17  19  21 23  25  28  32  33  35  37  42 43 44 47]*mm2m;
topo = [90 95 100 115 115 110 115 120 115 110 105 100 95 90 85 80];
dist_gw = [ 0 15 18 21 25 28 34 36 39 43 44 47]*mm2m;
gw = [90 95 100 105 110 110 105 100 95 90 85 80];