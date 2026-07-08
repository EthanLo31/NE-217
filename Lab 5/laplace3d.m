% laplace3d
% Solution to the Laplace equation over a cubic volume of points, 
% Solves for a point using the average of its neighbors.
% For the insulated boundary case specifically, 
% we must assume boundary condition are unknown rather than 0; 
%
% Parameters
% ==========
%    U:	Matrix of points over which solution will be approximated, including boundary conditions
%
% Return Values
% =============
%    U_out: The matrix Uout is the final approximated solution for all the points provided in U

function [U_soln] = laplace3d( U )
    % Error and Warning Checking
    % ==========================
    %
    % Check if argument is a 3D matrix, check if any edges are -infinity

    if ndims(U) ~= 3 || ~isnumeric(U)
        throw(MException('MATLAB:invalid_argument', 'The argument U must be a numeric 3D matrix'));
    end

    % 6 boundary faces of the 3D array
    faces = [ ...
    reshape(U(1,:,:), 1, []), ...   
    reshape(U(end,:,:), 1, []), ... 
    reshape(U(:,1,:), 1, []), ...  
    reshape(U(:,end,:), 1, []), ... 
    reshape(U(:,:,1), 1, []), ...  
    reshape(U(:,:,end), 1, [])     
    ];

    % Check if -Inf is present
    if any(faces == -Inf)
        throw(MException('MATLAB:invalid_argument', 'The argument U cannot have boundaries equal to -infinity'));
    end

    % Initialization
    % ==============
    %
    % Create the matrix U_out to store the solution, same size as U

    [n_x, n_y, n_z] = size( U );
    U_soln = U;

    % Mapping the unknown points to a unique number from 1 to m
    % =========================================================
    %
    % Associate each -inf to an integer from 1 to m, which will later be used in the equations.

    u_to_w = zeros( n_x, n_y, n_z );
    w_to_u = zeros( 3, n_x * n_y * n_z );
    m = 0;

    for i = 1:n_x
        for j = 1:n_y
            for k = 1:n_z
                if U(i, j, k) == -Inf
                    m = m + 1;
                    u_to_w(i, j, k) = m;
                    w_to_u(:, m) = [i, j, k]';
                end
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
    M = spalloc( m, m, 7*m );
    b = zeros( m, 1 );

    % Six neighbor dirs
    dir = [
        -1 0 0;
        1 0 0;
        0 -1 0;
        0 1 0;
        0 0 -1;
        0 0 1
        ];

    for k = 1:m
        
        % Coords of unknown
        x = w_to_u(1,k);
        y = w_to_u(2,k);
        z = w_to_u(3,k);

        % Loop through 6 neighbors
        for n = 1:6

            % Coords of neighbor
            i = x + dir(n,1);
            j = y + dir(n,2);
            l = z + dir(n,3);

            % Skip NaN points
            if isnan(U(r,c))
                continue;
            end

            % All other cases contribute to dia.
            M(k,k) = M(k,k) - 1;

            if U(i,j,l) == -Inf
                % Unknown
                p = u_to_w(i,j,l);
                M(k,p) = M(k,p) + 1;
            else
                % Dirichlet value
                b(k) = b(k) - U(i,j,l);
            end

    end

    w = M \ b;

    % Substituting the values back into the matrix U_out
    % ===================================================
    %
    % Substitute the values obtained (w) back into the matrix U_soln , this is your solution

    for n = 1:m
        % Copy the value from w into U_soln
        U_soln(w_to_u(1,n), w_to_u(2,n), w_to_u(3,n)) = w(n);
    end
end
