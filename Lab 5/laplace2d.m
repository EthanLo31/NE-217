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
    % Argument checking
    if ~ismatrix(U) || ~isnumeric(U)
        throw(MException('MATLAB:invalid_argument', 'The argument U must be a numeric 2D matrix'));
    end
    
    % Edges of the 2D Matrix
    edges = [U(1,:), U(end,:), U(:,1).', U(:,end).'];

    % Check if -inf is present
    if any(edges == -inf)
        throw(MException('MATLAB:invalid_argument', 'The argument U cannot have boundaries equal to -infinity'));
    end

    [n_x, n_y] = size( U );
    U_soln = U;

    % Step 2
    u_to_w = zeros( n_x, n_y );
    w_to_u = zeros( 2, n_x * n_y );
    m = 0;

    for j = 1:n_x
        for k = 1:n_y
            if U(j, k) == -Inf
                m = m + 1;
                u_to_w(j, k) = m;
                w_to_u(:, m) = [j, k]';
            end
        end
    end

    % Create the sparse system of linear equations
    M = spalloc( m, m, 5*m );
    b = zeros( m, 1 );

    for k = 1:m
        % Get the coordinates of the kth point
        % For each of the 4 adjacent points, determine if
        % the point is an insluated boundary point, a Dirichlet
        % boundary point or an unknown value and modify M as appropriate.
    end

    w = M \ b;

    for k = 1:m
        % Copy the value from w into U_soln
    end
end
