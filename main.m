%% Clean up
clearvars; close all; clc;

%% Constants

contact_resistance = 10^-3;
T_initial = 20+273;
T_inf = 220;
eps = 0.8;
h = 35;
k_p = 0.25;
k_m = 0.6;
k_b = 0.45;
c_p = 2000;
c_m = 3500;
c_b = 1060;
T_final = 53+273;

%% Initial conditions

dt = 0.1;
dr = 0.3125;

% Define spatial discretization
r = 0:dr:7.5; % radial distance from the center
n = length(r)+2; % number of spatial points

% Initialize temperature array
T = ones(n, 1) * T_initial;