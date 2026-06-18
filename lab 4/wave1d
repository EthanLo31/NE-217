% wave1d
% Approximates the one-dimensional wave equation using a finite-difference method.
%
% Parameters
% ==========
% c:        Wave speed
% x_int:    Spatial interval [a, b]
% n_x:      Number of points in the spatial interval
% t_int:    Time interval [t_initial, t_final]
% n_t:      Number of points in the time interval
% u_init:   Function handle for the initial displacement u(x, t_initial)
% du_init:  Function handle for the initial velocity du/dt at t_initial
% u_bndry:  Function handle for the boundary conditions
%           If a boundary is insulated, the function should return NaN
%           for that boundary.
%
% Return Values
% =============
% x_out: Column vector of x values
% t_out: Row vector of t values
% U_out: Matrix of approximations where U_out(i,k) approximates u(x_i,t_k)

function [x_out, t_out, U_out] = wave1d( c, x_int, n_x, t_int, n_t, u_init, du_init, u_bndry )

    % Check if c is scalar
    if ~isscalar(c)
        throw(MException('MATLAB:invalid_argument', ...
            'The argument c is not a scalar'));
    end

    % Check if c is positive
    if c <= 0
        throw(MException('MATLAB:invalid_argument', ...
            'The argument c must be positive'));
    end

    % Check if x_int is a 2-element vector
    if ~isvector(x_int) || length(x_int) ~= 2
        throw(MException('MATLAB:invalid_argument', ...
            'The argument x_int is not a 2-element vector'));
    end

    % Check if the x interval is valid
    if x_int(1) >= x_int(2)
        throw(MException('MATLAB:invalid_argument', ...
            'The entries of x_int must satisfy x_int(1) < x_int(2)'));
    end

    % Check if n_x is scalar
    if ~isscalar(n_x)
        throw(MException('MATLAB:invalid_argument', ...
            'The argument n_x is not a scalar'));
    end

    % Check if n_x is an integer
    if n_x ~= round(n_x)
        throw(MException('MATLAB:invalid_argument', ...
            'The argument n_x is not an integer'));
    end

    % Need enough points for the finite-difference method
    if n_x < 4
        throw(MException('MATLAB:invalid_argument', ...
            'The argument n_x must be at least 4'));
    end

    % Check if t_int is a 2-element vector
    if ~isvector(t_int) || length(t_int) ~= 2
        throw(MException('MATLAB:invalid_argument', ...
            'The argument t_int is not a 2-element vector'));
    end

    % Check if the time interval is valid
    if t_int(1) >= t_int(2)
        throw(MException('MATLAB:invalid_argument', ...
            'The entries of t_int must satisfy t_int(1) < t_int(2)'));
    end

    % Check if n_t is scalar
    if ~isscalar(n_t)
        throw(MException('MATLAB:invalid_argument', ...
            'The argument n_t is not a scalar'));
    end

    % Check if n_t is an integer
    if n_t ~= round(n_t)
        throw(MException('MATLAB:invalid_argument', ...
            'The argument n_t is not an integer'));
    end

    % Need at least two time points
    if n_t < 2
        throw(MException('MATLAB:invalid_argument', ...
            'The argument n_t must be at least 2'));
    end

    % Check if u_init is a function handle
    if ~isa(u_init, 'function_handle')
        throw(MException('MATLAB:invalid_argument', ...
            'The argument u_init is not a function handle'));
    end

    % Check if du_init is a function handle
    if ~isa(du_init, 'function_handle')
        throw(MException('MATLAB:invalid_argument', ...
            'The argument du_init is not a function handle'));
    end

    % Check if u_bndry is a function handle
    if ~isa(u_bndry, 'function_handle')
        throw(MException('MATLAB:invalid_argument', ...
            'The argument u_bndry is not a function handle'));
    end

    % calculate for stability ratio r
    h = (x_int(2) - x_int(1))/(n_x - 1);
    dt = (t_int(2) - t_int(1))/(n_t - 1);

    r = (c*dt/h)^2;

    % if r is greater than 1
    if r >= 1

        n_t_recommended = floor(c*(t_int(2) - t_int(1))/h) + 2;
        throw(MException('MATLAB:invalid_argument', ...
            'The ratio (c*dt/h)^2 = %.4f >= 1, use n_t = %d', ...
            r, n_t_recommended));
    end

    % initialize

    x_out = linspace(x_int(1), x_int(2), n_x)';
    t_out = linspace(t_int(1), t_int(2), n_t);
    U_out = zeros(n_x, n_t);
    initial_values = u_init(x_out);
    U_out(:,1) = initial_values(:);
    bndry = u_bndry(t_out);

    % left & right bound insulation 
    left_known = ~isnan(bndry(1,:));
    U_out(1,left_known) = bndry(1,left_known);
    right_known = ~isnan(bndry(2,:));
    U_out(end,right_known) = bndry(2,right_known);

    if isnan(bndry(1,1))
        U_out(1,1) = (4*U_out(2,1) - U_out(3,1))/3;
    end

    if isnan(bndry(2,1))
        U_out(end,1) = (4*U_out(end-1,1) - U_out(end-2,1))/3;
    end

    % solve for first iteration

    % slides use euler's method
    initial_velocity = du_init(x_out);
    initial_velocity = initial_velocity(:);
    U_out(2:end-1,2) = U_out(2:end-1,1) + dt*initial_velocity(2:end-1);

    % left & right insulation check
    if isnan(bndry(1,2))
        U_out(1,2) = (4*U_out(2,2) - U_out(3,2))/3;
    end

    if isnan(bndry(2,2))
        U_out(end,2) = (4*U_out(end-1,2) - U_out(end-2,2))/3;
    end

    % solve for rest

    for k = 2:n_t-1

        U_out(2:end-1,k+1) = 2*(1 - r)*U_out(2:end-1,k) - U_out(2:end-1,k-1) ...
            + r*(U_out(1:end-2,k) + U_out(3:end,k));

        % another insulation check
        if isnan(bndry(1,k+1))
            U_out(1,k+1) = (4*U_out(2,k+1) - U_out(3,k+1))/3;
        end

        if isnan(bndry(2,k+1))
            U_out(end,k+1) = (4*U_out(end-1,k+1) - U_out(end-2,k+1))/3;
        end
    end
end
