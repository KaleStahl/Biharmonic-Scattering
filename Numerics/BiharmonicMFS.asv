clear all;
clc;
close all;

%% Constants

N= 20;                                              % accuracy parameter
np = 2*N;                                           % # of points
k = 1:N;                                            % orders in Bessel expansion

%% Boundary Conditions

te = 1.5*pi*(.5:np-.5)'/np;                         % angles of boundary pts
re = 1./max(abs(sin(te)),abs(cos(te)));             % radii of boundary pts
ti = 1.5*pi*rand(np,1);                             % angles of interior pts
ri = rand(np,1)./max(abs(sin(ti)),abs(cos(ti)));    % radii interior pts
t = [te; ti];                                      % combined boundary and interior angles
r = [re; ri]';                                      % combined boundary and interior radii

%% Set up A and Find Relative Subspace Angles

lambda_vec = linspace(.1, 25, 125);                              % trial values of lambda
S = [];
for rj = r
    sj = [];
    for lambda = lambda_vec
        A = sin(2*t*k/3).*besselj(2*k/3, sqrt(lambda)*rj);
        [Q,R] = qr(A,0);
        s = min(svd(Q(1:np,:)));                            % subspace angle for this lambda

        sj =[sj s];
    end
    S = [S; sj];
end
%% Convert to signed subspace angles

I = 1:length(lambda_vec);               % all lambda points
J = I(2 : end-1);                       % interior points
J = J( S(J)<S(J-1) & S(J)<S(J+1) );     % local minima via smart indexing
J = J+ (S(J-1)>S(J+1));                 % points where sign changes
K = 0*I; 
K(J) = 1;
S = S.*(-1).^cumsum(K);               % introduce sign flips

%% Find Eigenvalues

for j = 1:length(J)
    if (J(j)+9 >125 )
        index = J(j):125;
    else
        index = J(j):J(j)+9;
    end
    index_lambda = lambda_vec(index);
    P = polyfit(S(index)/norm(S(index)), index_lambda, 9);
    lambda_min = polyval(P, 0);
    disp(lambda_min) % display eigenvalue? I think?
end