% laplace2d
% Solution to the Laplace equation using a square grid of points, 
% we can model irregular shapes and map them out in matrix form. 
% Solves for a point using the average of its neighbors.
% For the insulated boundary case specifically, 
% we must assume boundary condition are unknown rather than 0; 
% often leading to a system of linear equations.
%
% Parameters
% ==========
%    U:	Matrix of points over which solution will be approximated, including boundary conditions
%
% Return Values
% =============
%    U_out: The matrix Uout is the final approximated solution for all the points provided in U

function [U_soln] = laplace2d( U )
    % Error and Warning Checking
    % ==========================
    %
    % Check if argument is a 2D matrix, check if any edges are -infinity

    if ~ismatrix(U) || ~isnumeric(U)
        throw(MException('MATLAB:invalid_argument', 'The argument U must be a numeric 2D matrix'));
    end
    
    % Edges of the 2D Matrix
    edges = [U(1,:), U(end,:), U(:,1).', U(:,end).'];

    % Check if -Inf is present
    if any(edges == -Inf)
        throw(MException('MATLAB:invalid_argument', 'The argument U cannot have boundaries equal to -infinity'));
    end

    % Initialization
    % ==============
    %
    % Create the matrix U_out to store the solution, same size as U

    [n_x, n_y] = size( U );
    U_soln = U;

    u_to_w = zeros( n_x, n_y );
    w_to_u = zeros( 2, n_x * n_y );
    m = 0;

    % Mapping the unknown points to a unique number from 1 to m
    % =========================================================
    %
    % Associate each -inf to an integer from 1 to m, which will later be used in the equations.

    for j = 1:n_x
        for k = 1:n_y
            if U(j, k) == -Inf
                m = m + 1;
                u_to_w(j, k) = m;
                w_to_u(:, m) = [j, k]';
            end
        end
    end

    % Creating and solving a system of linear equations
    % =================================================
    %
    % Create a column vector b with size m, and a matrix M with size m by m.
    % For each unknown point i, take all its neighbors. 
    % If it’s NaN: ignore, 
    % If it’s a Dirichlet condition, subtract 1 from the ith diagonal entry of M, 
    % and subtract the value from the ith entry of b. 
    % If its another unknown, subtract 1 from the diagonal entry of M 
    % and add 1 to the (i, j)th entry of M. 
    % Solve the system w = M \ b.

    % Create the sparse system of linear equations
    M = spalloc( m, m, 5*m );
    b = zeros( m, 1 );

    for k = 1:m
        % Coords of kth unknown
        row = w_to_u(1, k);
        col = w_to_u(2, k);

        % Four neighbors
        neighbors = [row-1, col;
                     row+1, col;
                     row, col-1;
                     row, col+1];
        
        % Loop over all neighbors
        for n = 1:4
            % Coords of neighbors
            r = neighbors(n,1);
            c = neighbors(n,2);

            % Skip NaN points
            if isnan(U(r,c))
                continue;
            end

            % All other cases contribute to dia.
            M(k,k) = M(k,k) - 1;

            if U(r,c) == -Inf
                % Unknown
                j = u_to_w(r,c);
                M(k,j) = M(k,j) + 1;
            else
                % Dirichlet value
                b(k) = b(k) - U(r,c);
            end
        end
    end

    w = M \ b;

    % Substituting the values back into the matrix U_out
    % ===================================================
    %
    % Substitute the values obtained (w) back into the matrix U_soln , this is your solution

    for k = 1:m
        % Copy the value from w into U_soln
        U_soln(w_to_u(1,k), w_to_u(2,k)) = w(k);
    end
end
