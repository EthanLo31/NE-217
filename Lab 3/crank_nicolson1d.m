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

    % Safely pre-fill boundaries only if they are NOT NaN (Dirichlet)

    left_valid = ~isnan(bndry(1, :));

    U_out(1, left_valid) = bndry(1, left_valid);

    right_valid = ~isnan(bndry(2, :));

    U_out(end, right_valid) = bndry(2, right_valid);

    % Step 3: Solve

    % Setup standard tridiagonal matrices

    a_diag = (1 + r) * ones(nx - 2, 1);
    a_off = (-r / 2) * ones(nx - 3, 1);
    A = diag(a_diag) + diag(a_off, 1) + diag(a_off, -1);

    b_diag = (1 - r) * ones(nx - 2, 1);
    b_off = (r / 2) * ones(nx - 3, 1);
    B = diag(b_diag) + diag(b_off, 1) + diag(b_off, -1);

    % Modify Matrices for Insulated Boundaries
    % This forces the matrix to inherently solve the zero-flux condition
    if isnan(bndry(1, 1))

        A(1, 1) = 1 + r / 2;

        B(1, 1) = 1 - r / 2;

    end

    if isnan(bndry(2, 1))

        A(end, end) = 1 + r / 2;

        B(end, end) = 1 - r / 2;

    end

    % Iterate through the matrix, avoid boundary conditions
    for k = 1:nt - 1
        % Extract current state
        U_current = U_out(2:nx - 1, k);

        bndry_term = zeros(nx - 2, 1);

        % Left Boundary RHS Handling

        if isnan(bndry(1, k + 1))

            bndry_term(1) = 0; % Math is absorbed into the modified matrix

        else

            bndry_term(1) = (r / 2) * (U_out(1, k) + U_out(1, k + 1));

        end

        % Right Boundary RHS Handling

        if isnan(bndry(2, k + 1))

            bndry_term(end) = 0; % Math is absorbed into the modified matrix

        else

            bndry_term(end) = (r / 2) * (U_out(end, k) + U_out(end, k + 1));

        end

        % Use matrix division to calculate values
        U_out(2:nx - 1, k + 1) = A \ (B * U_current + bndry_term);
        % Apply zero-flux condition to boundary nodes for accurate plotting
        if isnan(bndry(1, k + 1))
            U_out(1, k + 1) = U_out(2, k + 1);
        end

        if isnan(bndry(2, k + 1))
            U_out(end, k + 1) = U_out(end - 1, k + 1);
        end

    end

end
