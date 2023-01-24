function dtheta = KS_fcn_RK(~,theta,Kij,n,k,Aij,omega)

    theta_i = theta * ones(1,n);                            % (nxn) matrix with same theta values in each row
    theta_j = ones(n,1) * theta';                           % (nxn) matrix with same theta values in each col
    theta_ij = theta_i - theta_j;                           % (nxn) matrix of all possible (theta_i - theta_j) combinations
    dtheta = omega + k.*sum(Kij.*sin(theta_ij - Aij),1)';   % Kuramoto Sakaguchi Equation in vector form. (k is k/n)
end