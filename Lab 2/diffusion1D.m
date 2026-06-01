% diffusion1D
% Approximates the heat diffusion equation in 1D
% 
% Parameters 
% ========== 
%    kappa: Thermal conductivity of the material
%    x_rng: Interval in 1D space over which we’re analyzing 
%    t_rng: Time over which it’s being analyzed 
% 
%    u_init: Initial heat distribution in the rod
%    u_bndry: Boundary conditions
% 
%    nx: Number of points on X over which solution will be approximated
%    nt: Number of time points over which solution is approximated
% 
% Return Values 
% ============= 
%    x_out: List of points in X over which solution was approximated
%    t_out: List of time points over which solution was approximated
%    U_out: Contains all the solution approximations for all points (x_out, t_out)

 

function [x_out, t_out, U_out] = diffusion1D( kappa, x_rng, nx, t_rng, nt, u_init, u_bndry )

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

    % Step 1: Error checking
    % Ensure ratio is proper for stability
    % Calculate h
    h = (x_rng(2) - x_rng(1))/(nx - 1);

    % Calculate dt
    dt = (t_rng(2) - t_rng(1))/(nt - 1);

    % Calculate ratio
    r = ((kappa*dt)/h^2);

    % If ratio is too large
    if r > 0.5
        % Calculate smallest nt only if necessary (performance overhead)
        s_nt = ceil((2*kappa*(t_rng(2) - t_rng(1)))/(h^2) + 1);
        % Tell user
        throw(MException('MATLAB:invalid_argument', sprintf('The ratio kappa*dt/h^2 = %.4f. The smallest nt value for stability would be %d', r, s_nt)));
    end

    % Step 2: Initialization
    % Setup the vectors and matrices that will be solved later. 
    % Create column vector of x values
    x_out = linspace(x_rng(1), x_rng(2), nx)';

    % Create row vector of t values
    t_out = linspace(t_rng(1), t_rng(2), nt);

    % Create nx by nt matrix
    U_out = zeros(nx, nt);

    % Fill initial conditions
    U_out(:,1) = u_init(x_out);

    % Bounds
    bndry = u_bndry(t_out);

    % Left
    U_out(1,:) = bndry(1,:);
    % Right
    U_out(end,:) = bndry(2,:);

    % Step 3: Solve
    % Iterate through the matrix, avoid boundary conditions
    for k = 1:nt-1
        for i = 2:nx-1
            % Use formula to calculate values
            U_out(i,k+1) = U_out(i,k) + r * (U_out(i-1,k) - 2*U_out(i,k) + U_out(i+1,k));
        end
    end

end
