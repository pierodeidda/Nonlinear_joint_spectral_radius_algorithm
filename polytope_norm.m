function [Y_hull, vertices, info] = polytope_norm(F, x0, rho, Nmax)
%POLYTOPE_METHOD Construct an extremal polytope norm for a family of 
%                linear and order-preserving maps.
%
% INPUT:
%   F    : cell array of function handles
%          e.g. F = {F_1, F_2}
%
%   x0   : dominant eigenvector of spectrum maximizing matrix product
%
%   rho  : joint spectral radius (known)
%
%   Nmax : maximum number of generations
%
% OUTPUT:
%   Y_hull : vertices of the final convex hull
%
%   info : structure containing information about the algorithm
%
%          info.converged
%          info.iterations
%          info.num_vertices
%          info.rho
% -------------------------------------------------------------------------

% Number of maps
    m = length(F);

    % Check dimension
    x0 = x0(:);
    if length(x0) ~= 2
        error('POLYTOPE_METHOD currently requires x0 to be two-dimensional.');
    end

    % NORMALIZATION of the family of matrices
    F_norm = cell(1,m);
    for q = 1:m
        F_norm{q} = @(x) F{q}(x) / rho;
    end

    % INITIALIZATION
    % Initial set of points:
    % (0,0), (x0(1),0), x0, (0,x0(2))
    Y = [[0;0], [x0(1);0], [0;x0(2)], x0];
    vertices=Y;
    % Initial convex hull
    KH = convhull(Y(1,:), Y(2,:));
    Y_hull = Y(:,KH(1:end-1));
    Y_hull = reorder_hull(Y_hull);

    % The points that will be propagated in the next generation.
    % Initially we only propagate x0.
    current_generation = x0;

    % INFORMATION
    info.converged = false;
    info.iterations = 0;
    info.num_vertices = zeros(Nmax,1);
    info.rho = rho;

    for i = 1:Nmax
        % We propagate every point of the current generation.
        next_generation = [];
        for j = 1:size(current_generation,2)
            x = current_generation(:,j);
            for q = 1:m
                % APPLY NORMALIZED MAP
                y = F_norm{q}(x);
                % CHECK IF y IS INSIDE CURRENT CONVEX HULL
                [in,on] = inpolygon(y(1), y(2),Y_hull(1,:), Y_hull(2,:));
                % IF y IS OUTSIDE
                if ~(in || on)
                    % KEEP y AND ITS AXIS PROJECTIONS
                    vertices = [[y(1);0],[0;y(2)], vertices, y];
                    % UPDATE CONVEX HULL
                    KH = convhull(vertices(1,:), vertices(2,:), "Simplify", true);
                    Y_hull = vertices(:,KH(1:end-1));
                    Y_hull = reorder_hull(Y_hull);
                    % SAVE THE LIST OF VERITICES THAT HAVE BEEN PRESERVED IN
                    % THE SAME ORDER THEY APPEARED IN Y
                    hidx = KH(1:end-1);         
                    vertices = vertices(:, sort(hidx));
                    % vertices = vertices(:, 4:end);
                    % Keep only the extremal points
                    % Y = Y_hull;
                    % y belongs to the next generation
                    next_generation = [next_generation, y];
                end
            end
        end
        
        % INFORMATION
        info.iterations = i;
        info.num_vertices(i) = size(Y_hull,2);
        % CHECK CONVERGENCE
        if isempty(next_generation)
            info.converged = true;
            info.iterations = i;
            info.num_vertices = info.num_vertices(1:i);
            fprintf('\nPolytope is stable.\n');
            fprintf('Iterations: %d\n', i);
            fprintf('Number of vertices: %d\n', size(Y_hull,2));
            fprintf('rho = %.16g\n', rho);
            return
        end
        % CHECK NUMERICAL STABILITY OF THE HULL
        % The hull at consecutive generations is compared to check if there
        % is significative difference
        if i > 1
            % Compare the number of vertices first.
            % If the number of vertices is the same, compare their
            % coordinates.
            if size(Y_hull,2) == size(Y_hull_previous,2)
                Diff = max(abs(Y_hull - Y_hull_previous),[],'all');
                % If numerically the two convex hulls are the same stop
                if Diff < 1e-10
                    info.converged = true;
                    info.iterations = i;
                    info.num_vertices = info.num_vertices(1:i);
                    fprintf('\nPolytope is stable.\n');
                    fprintf('Iterations: %d\n', i);
                    fprintf('Number of vertices: %d\n', size(Y_hull,2));
                    fprintf('rho = %.16g\n', rho);
                    return
                end
            end
        end
        % Save current hull for the next iteration
        Y_hull_previous = Y_hull;
        % NEXT GENERATION
        current_generation = next_generation;
     end


    % MAXIMUM NUMBER OF ITERATIONS
    info.converged = false;
    info.iterations = Nmax;
    info.num_vertices = info.num_vertices(1:Nmax);
    fprintf('\nMaximum number of iterations reached.\n');
    fprintf('Iterations: %d\n', Nmax);
    fprintf('Number of vertices: %d\n', size(Y_hull,2));
    fprintf('rho = %.16g\n', rho);
end