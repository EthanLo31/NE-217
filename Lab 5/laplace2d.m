function [U_soln] = laplace2d( U )
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
