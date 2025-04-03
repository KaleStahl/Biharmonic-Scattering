
%%
n = 2;
k = 1;
r = 1;
theta = pi/2;

%%
A = [ 
    besselh(n, k*r),                                besselh(n, 1i*k*r),                                     -besselj(n, k^4*sqrt(n)*r),                                                 -besselj(n, 1i*k^4*sqrt(n)*r);
    .5*k*(besselh(n-1, k*r)-besselh(n+1, k*r)),     .5*1i*k*(besselh(n-1, 1i*k*r)-besselh(n+1, 1i*k*r)),    k^4*sqrt(n)*.5*(besselj(n-1, k^4*sqrt(n)*r)-besselj(n+1, k^4*sqrt(n)*r)),   1i*k^4*sqrt(n)*.5*(besselj(n-1, 1i*k^4*sqrt(n)*r)-besselj(n+1, 1i*k^4*sqrt(n)*r));
    -k^2*besselh(n, k*r),                           k^2*besselh(n, 1i*k*r),                                 k^2*sqrt(n)*besselj(n, k^4*sqrt(n)*r),                                      -k^2*sqrt(n)*besselj(n, 1i*k^4*sqrt(n)*r);
    -.5*k^3*(besselh(n-1, k*r)-besselh(n+1, k*r)),  .5*1i*k^3*(besselh(n-1, 1i*k*r)-besselh(n+1, 1i*k*r)),  -k^6*n*.5*(besselj(n-1, k^4*sqrt(n)*r)-besselj(n+1, k^4*sqrt(n)*r)),        1i*k^6*n*.5*(besselj(n-1, 1i*k^4*sqrt(n)*r)-besselj(n+1, 1i*k^4*sqrt(n)*r));
    ]*exp(1i*n*theta);

y = [
    -exp(1i*k*r*theta);
    -1i*k*theta*exp(1i*k*r*theta);
    -1i*k*r*exp(1i*k*r*theta);
    -1i*k*theta*exp(1i*k*r*theta) + k^2*r*theta*exp(1i*k*r*theta);
    ];

x = A\y
