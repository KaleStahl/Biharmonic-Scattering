
%% Constants
N= 100;                         % Number of discretized points
Np = 100;
n = 2;                          % relative refractive index
p = -5:5;                       % Number of terms to approximate
k = 5;                          % wave number
r = 1;                          % radius of object
theta = linspace(0, 2*pi, N);   % Boundary data
phi = linspace(0, 2*pi, N);     % Test Data


%% Testing Reciprocity Relation

u_infty_1 = findFarField(theta, phi, p, n, k, r);

u_infty_2 = findFarField(phi+pi, theta+pi, p, n , k, r);

%% Discretize Far-Field Operator

[Theta, Phi] = meshgrid(theta, phi);
F = 2*pi/N*findFarField(Theta, Phi, p, n, k, r);
g = @(z, t) exp(1i*k*dot(z, [cos(t), sin(t)]));
g_vec = @(z) arrayfun(@(t) g(z, t), theta).';
D = @(z) dot(F*g_vec(z), g_vec(z));
[X, Y] = meshgrid(linspace(-2*r, 2*r, Np));
z = cat(3, X, Y);
Z = zeros(Np);

%%
for i = 1:Np
    for j = 1:Np
        z = [X(i, j), Y(i, j)];
        Z(i, j) = D(z);
    end
end
%%

surf(X, Y, abs(Z))




