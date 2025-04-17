clc;
clear;
close all;

%% Constants
N= 1;                           % Number of boundary points
n = 2;                          % relative refractive index
p = -5:5;                       % Number of terms to approximate
k = 5;                          % wave number
r = 2;                          % radius of collection
theta = linspace(0, 2*pi, N);   % Boundary data
phi = linspace(0, 2*pi, N);     % Test Data


%% Computing Far Field Pattern

u_infty_1 = findFarField(theta, phi, p, n, k, r)

u_infty_2 = findFarField(-phi, -theta, p, n , k, r)
