clear; 
close all; 
clc

%% Constants
k = linspace(-1, 1, 1000);
n = -5:5;
R = 1;

%% Function

gamma_n = @(n) 1i*k.*(besselh(n-1,1, 1i*k*R)- besselh(n,1, 1i*k*R)*n./(1i.*k*R))./besselh(n,1,  1i*k*R);

%% Plots

fo = figure;
set(0,'defaulttextinterpreter','latex');
set(0,'defaultLineLineWidth',2);
set(0,'DefaultAxesFontSize',14);
set(0,'DefaultLegendInterpreter','latex');


hold on;
for dim = n
    gamma = gamma_n(dim);
    p = plot(k, gamma, 'DisplayName', ['n=', num2str(dim)]);
end
for dim = n
    limit = -1*abs(dim*ones(1, max(size(k))))/R;
    q = plot(k, limit, 'black--', 'DisplayName', '');
end
zero = xline(0, 'black--','DisplayName', '');
xlabel('$k$');
ylabel('$\gamma_n$');
hold off;
