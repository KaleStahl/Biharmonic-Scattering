clear; 
close all;
clc;

%% Constants
k = linspace(-1/100, 1/100, 10000);
n = -5:5;
R = 1;

%% Function

gamma_n = @(n) 1i*k.*(besselh(n-1,1, 1i*k*R)- besselh(n,1, 1i*k*R)*n./(1i.*k*R))./besselh(n,1,  1i*k*R);

%% Plots
set(0,'defaulttextinterpreter','latex');
set(0,'defaultLineLineWidth',2);
set(0,'DefaultAxesFontSize',14);
set(0,'DefaultLegendInterpreter','latex');

for dim = n
    fo = figure;
    hold on;
    gamma = gamma_n(dim);
    p = plot(k, gamma, 'DisplayName', ['n=', num2str(dim)]);
    limit = -1*abs(dim*ones(1, max(size(k))))/R;
    q = plot(k, limit, 'black--', 'DisplayName', '');
    zero = xline(0, 'black--');
    title(append('Graph of $\gamma_n$ for n=', num2str(dim)))
    xlabel('$k$');
    ylabel('$\gamma_n$');
    hold off;
    file_name = append('gamma_n', num2str(dim), '.png');
    saveas(fo, file_name, 'png'); 
end
