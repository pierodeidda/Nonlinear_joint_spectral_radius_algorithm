function is_ordered = order_check(x,Y,strict)
% ORDER_CHECK: checks whether x is dominated by some column of Y.
%
% INPUT:
%   x      : vector
%   Y      : matrix whose columns are the points to be checked
%   strict : false -> x <= y
%            true  -> x <= y and x ~= y (numerically up to 1e-10)
%
% OUTPUT:
%   is_ordered : 
%               logical value - TRUE if x<=y (+ x~=y if strict=true)
%                             - FALSE otherwise
% -------------------------------------------------------------------------
% Calculate the difference matrix
    D = Y - x;

    % Check componentwise order
    dominated = max(min(D)) > 0;

    % Check for numerical duplicates 
    duplicate = min(vecnorm(D)) < 1e-10;

    if strict
        % Numerical duplicates are not considered as ordered
        is_ordered = any(dominated & ~duplicate);
    else
        % Also regard numerical duplicates as already ordered
        duplicate = min(vecnorm(D)) < 1e-10;
        is_ordered = dominated || duplicate;
    end

end