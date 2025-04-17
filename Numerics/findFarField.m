function u_infty = findFarField(theta, phi, p_range, n, k, r)
    u_infty = 0;
    for p = p_range
        x = findCoeffs(n, p ,k, r, theta, phi);
        a_p = x(1);
        u_infty = u_infty + 4/1i*a_p*exp(1i*p*(theta - phi));
    end
end