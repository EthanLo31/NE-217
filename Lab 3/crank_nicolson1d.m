% crank_nicolson1d
% Approximates the heat diffusion equation in 1D
%
% Parameters
% ==========
% kappa: Thermal conductivity of the material
% x_rng: Interval in 1D space over which we’re analyzing
% t_rng: Time over which it’s being analyzed
%
% u_init: Initial heat distribution in the rod
% u_bndry: Boundary conditions
%
% nx: Number of points on X over which solution will be approximated
% nt: Number of time points over which solution is approximated
%
% Return Values
% =============
% x_out: List of points in X over which solution was approximated
% t_out: List of time points over which solution was approximated
% U_out: Contains all the solution approximations for all points (x_out, t_out)

function [x_out, t_out, U_out] = crank_nicolson1d(kappa, x_rng, nx, t_rng, nt, u_init, u_bndry)

    % Check if kappa is scalar
    if ~isscalar(kappa)
        throw(MException('MATLAB:invalid_argument', 'The argument kappa is not a scalar'));
    end

    % Check if x_rng is 2-element vector
    if ~isvector(x_rng) || length(x_rng) ~= 2
        throw(MException('MATLAB:invalid_argument', 'The argument x_rng is not a 2-element vector'));
    end

    % Check if nx is scalar
    if ~isscalar(nx)
        throw(MException('MATLAB:invalid_argument', 'The argument nx is not a scalar'));
    end

    % Check if t_rng is 2-element vector
    if ~isvector(t_rng) || length(t_rng) ~= 2
        throw(MException('MATLAB:invalid_argument', 'The argument t_rng is not a 2-element vector'));
    end

    % Check if nt is scalar
    if ~isscalar(nt)
        throw(MException('MATLAB:invalid_argument', 'The argument nt is not a scalar'));
    end

    % Check if u_init is function handle
    if ~isa(u_init, 'function_handle')
        throw(MException('MATLAB:invalid_argument', 'The argument u_init is not a function handle'));
    end

    % Check if u_bndry is function handle
    if ~isa(u_bndry, 'function_handle')
        throw(MException('MATLAB:invalid_argument', 'The argument u_bndry is not a function handle'));
    end

    % Step 1: Variable setup
    % Calculate h
    h = (x_rng(2) - x_rng(1)) / (nx - 1);

    % Calculate dt
    dt = (t_rng(2) - t_rng(1)) / (nt - 1);

    % Calculate ratio
    r = ((kappa * dt) / h ^ 2);

    % Step 2: Initialization
    % Setup the vectors and matrices that will be solved later.
    % Create column vector of x values
    x_out = linspace(x_rng(1), x_rng(2), nx)';

    % Create row vector of t values
    t_out = linspace(t_rng(1), t_rng(2), nt);

    % Create nx by nt matrix
    U_out = zeros(nx, nt);

    % Fill initial conditions
    U_out(:, 1) = u_init(x_out);

    % Bounds
    bndry = u_bndry(t_out);

    % Left
    U_out(1, :) = bndry(1, :);
    % Right
    U_out(end, :) = bndry(2, :);

    % Step 3: Solve
    % Setup tridiagonal matrices
    a_diag = (1 + r) * ones(nx - 2, 1);
    a_off = (-r / 2) * ones(nx - 3, 1);
    A = diag(a_diag) + diag(a_off, 1) + diag(a_off, -1);

    b_diag = (1 - r) * ones(nx - 2, 1);
    b_off = (r / 2) * ones(nx - 3, 1);
    B = diag(b_diag) + diag(b_off, 1) + diag(b_off, -1);

    % Iterate through the matrix, avoid boundary conditions
    for k = 1:nt - 1
        % Extract current state
        U_current = U_out(2:nx - 1, k);

        % Calculate boundary terms
        bndry_term = zeros(nx - 2, 1);
        bndry_term(1) = (r / 2) * (U_out(1, k) + U_out(1, k + 1));
        bndry_term(end) = (r / 2) * (U_out(end, k) + U_out(end, k + 1));

        % Use matrix division to calculate values
        U_out(2:nx - 1, k + 1) = A \ (B * U_current + bndry_term);
    end

end
