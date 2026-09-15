function [vertices, rho, info] = polytope_prenorm(F, x0, rho, Nmax)
%POLYTOPAL_PRENORM Iteratively builds an extremal prenorm for a family
%                  of continuous, homogenous and order-preserving maps.
%
% INPUT:
%   F    : cell array of function handles
%          e.g. F = {F_1, F_2, F_3}
%
%   x0   : initial point (possibly dominant eigenvector of the presumed spectrum 
%          maximizing composition)
%
%   rho  : presumed joint spectral radius
%
%   Nmax : maximum number of iterations for each construction
%
% OUTPUT:
%   Y    : matrix containing all relevant vertices generated in the
%          final construction
%
%   info : structure containing information about the algorithm
%          info.converged
%          info.iterations (iterations required by the last restart)
%          info.restarts (number of times the JSR is updated)
% -------------------------------------------------------------------------

% Number of maps in the family
    m = length(F);
% Number of restarts
    restart = 0;

    while true
        % NORMALIZATION
        % At this stage rho is the proposed JSR. Each individual map
        % is normalized by rho.
        normalization = rho;
        F_norm = cell(1,m);
        for q = 1:m
            F_norm{q} = @(x) F{q}(x) / normalization;
        end

        % INITIALIZATION
        vertices = x0(:);
        new_vertices = x0(:);
        next_vertices = [];

        % Path associated with each point
        % path{j} = [k_1,...,k_l]
        % means
        % Y(:,j) = F_{k_l} o ... o F_{k_1}(x0)
        path = cell(1,1);
        path{1} = [];
        next_path = cell(1,1);

        % Store information
        info.restarts = restart;
        
        % MAIN LOOP
        restart_required = false;
        for i = 1:Nmax
            % Apply all maps to all new vertices
            for j = 1:size(new_vertices,2)
                    xj = new_vertices(:,j);
                    % if next_vertices nonempty and some column dominates xj, skip this xj
                    if ~isempty(next_vertices)
                        dominated_by_next = all(next_vertices - xj >= 0, 1); 
                        if any(dominated_by_next)
                            continue;   % skip to next j (do not enter inner for q = 1:m)
                        end
                    end
                for q = 1:m
                    % Apply q-th normalized map
                    x = F_norm{q}(xj);
                    % CHECK IF x > x0
                    % x > x0 means that x0 is dominated by x and the JSR 
                    % needs to be updated.
                    if order_check(x0, x, true)
                        % The path of maps generating x is
                        % path{j} followed by q
                        current_path = [path{j}, q];
                        % Length of the composition
                        k = length(current_path);
                        % BUILD THE COMPOSITION
                        G = @(z) z;
                        for l = 1:k
                            q_l = current_path(l);
                            G_old = G;
                            G = @(z) F{q_l}(G_old(z));
                        end
                        % POWER METHOD to compute the spectral radius and
                        % dominant eigenvector of G
                        [x_new,rho_map] = pw_method(G,length(x0),1000, 1e-15);
                        % Check whether power method succeeded
                        if any(isnan(x_new))
                            error(['Power method failed for the ', ...
                                   'composition corresponding to path ', ...
                                   mat2str(current_path)]);
                        end
                        % UPDATE THE JSR and restart
                        rho_new = rho_map^(1/k);
                        fprintf('\nRestart required.\n');
                        fprintf('Path: %s\n',mat2str(current_path));
                        fprintf('Length of composition: %d\n',k);
                        fprintf('New presumed JSR: %.16g\n',rho_new);
                        % RESTART
                        x0 = x_new;
                        rho = rho_new;
                        restart = restart + 1;
                        restart_required = true;
                        break
                    % CHECK IF x IS A NEW RELEVANT POINT
                    elseif ~order_check(x,vertices,false)
                            % remove vertices dominated by x
                            dom_vertices = all((x - vertices)>0, 1); 
                            vertices=vertices(:, ~dom_vertices);
                            % mask true for columns dominated by x
                            if ~isempty(next_vertices)
                                dom = all((x - next_vertices) > 0, 1);
                                next_vertices = next_vertices(:, ~dom);
                                next_path = next_path(~dom);
                            end
                            % save x as a new vertex
                            l = size(vertices,2);
                            vertices(:,l+1) = x;
                            % Store x and how it was generated for the next
                            % iteration
                            l = size(next_vertices,2);
                            next_vertices(:,l+1) = x;
                            next_path{l+1} = [path{j},q];
                    end
                end
                if restart_required
                    break
                end
            end

            % If a point x > x0 was found, restart everything
            if restart_required
                break
            end

            % CONVERGENCE
            if size(next_vertices,2) == 0
                info.converged = true;
                info.iterations = i;
                info.restarts = restart;

                fprintf('\nPolytopal prenorm exists.\n');
                fprintf('Iterations: %d\n',i);
                fprintf('Number of vertices: %d\n',size(vertices,2));
                fprintf('rho = %.16g\n',rho);
                fprintf('Number of restarts: %d\n',restart);
                return
            end

            % NEXT GENERATION
            new_vertices=next_vertices;
            next_vertices=[];
            path=next_path;
            next_path=cell(1,1);
        end

        % MAXIMUM NUMBER OF ITERATIONS
        if ~restart_required
            info.converged = false;
            info.iterations = Nmax;
            info.restarts = restart;
            fprintf('\nMaximum number of iterations reached.\n');
            return
        end
        % Otherwise the while loop restarts with the new x0 and rho.
    end
end