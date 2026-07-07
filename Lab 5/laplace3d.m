function [U_soln] = laplace3d( U )
    [n_x, n_y, n_z] = size( U );
    U_soln = U;

    % Step 2
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

    % Create the sparse system of linear equations
    M = spalloc( m, m, 7*m );
    b = zeros( m, 1 );

    for k = 1:m
        % Get the coordinates of the kth point
        % For each of the 6 adjacent points, determine if
        % the point is an insluated boundary point, a Dirichlet
        % boundary point or an unknown value and modify M as appropriate.
    end

    w = M \ b;

    for k = 1:m
        % Copy the value from w into U_soln
    end
end
