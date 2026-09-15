function [x,rho] = pw_method(F, n, max_iter, tol)
% NONLINEAR POWER METHOD FOR A NONLINEAR, HOMOGENEOUS AND ORDER PRESERVING MAP F
%
% INPUT: 
%        F: map
%        n: dimension of the space
%        max_iter: max number of iterations
%        tol: tolerance on the normalized eigenvalue equation for stopping 
% OUTPUT: 
%        x: dominant eigenvector
%        rho: spectral radius
% -------------------------------------------------------------------------

x=ones(n,1);

for j=1:max_iter
    if j==max_iter
        % IF ALGORITHM DIDN'T CONVERGE LET THE USER KNOW
        fprintf('Power Method did not converge');
        x=NaN;
        rho=NaN;
        break
    end
    % ITERATIVE STEP 
    % i) APPLY MAP 
    x=F(x);
    % ii) NORMALIZE
    x=x/norm(x,1);
    % iii) COMPUTE NORM
    rho=(norm(F(x),1));

    if norm(F(x)-rho*x,1)/rho<tol
       % IF EIGENVALUE EQ. IS SATISFIED STOP
       break
    end
end
