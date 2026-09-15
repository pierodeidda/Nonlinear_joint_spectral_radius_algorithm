function Y_hull = reorder_hull(Y_hull)
%%  REORDER THE POINTS OF A CONVEX HULL SO THAT THE FIRST IS ZERO
%   INPUT: LIST OF POINTS CONTAINING ZERO
%   OUTPUT: ORDERD LIST OF POINTS
% -------------------------------------------------------------------------

    idx0 = find(Y_hull(1,:) == 0 & Y_hull(2,:) == 0, 1);

    if ~isempty(idx0)
        n = size(Y_hull,2);
        Y_hull = Y_hull(:, [idx0:n, 1:idx0-1]);
    end

end