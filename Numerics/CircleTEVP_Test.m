clear all;
clc;
close all;

%% Constants

N_I = 20;                   % number of interior points
N_B = 20;                   % Number of boundary points
N_K = 30;                   % number of bessel expansions to run through
k = linspace(0, 10, N_K);   % bessel expansion numbers

%% Set up Boundary Conditions

theta_i = rand(N_I, 1);     % randomly selected interior angles
theta_b = rand(N_B, 1);     % randomly selected boundary angles
r_i = rand(N_I, 1);         % randomly selected interior radii
x_b = cos(theta_b);
y_b = sin(theta_b);
x_i = r_i.*cos(theta_i);
y_i = r_i.*sin(theta_i);

%% Set up A

A = zeros(N_B, N_B+N_I, N_K);

% Loop for boundary
for i = 1:N_B
    for j = 1:N_I
        for l = 1:N_K
            A(i, j, l) = besselj(i, k(l))*cos(i*theta_b(i));
        end
    end
end

% Loop for interior
for i = 1:N_B
    for j = 1:N_I
        for l = 1:N_K
            A(i, j+N_B, l) = besselj(i, k(l)*r_i(j))*cos(i*theta_i(j));
        end
    end
end

%% QR factorization

vect = zeros(N_K, 1);
for i =1:N_K
    [Q,R] = qr(A(:, :, i));
    vect(i) = min(svd(Q)); 
end

%% Plot

plot(vect);