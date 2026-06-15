%Base Parameters
kappa = 1;
x_rng = [0, pi];
t_rng = [0, 2];

%Base Case
nx_base = 10;
nt_base = 20;
[x3b, t3b, U3b] = crank_nicolson1d(kappa, x_rng, nx_base, t_rng, nt_base, @u3b_init, @u3b_bndry);

% Plot 1
figure;
mesh(t3b, x3b, U3b);
title('Heat Distribution 3.8b: e33lo, ylepage, jteeter, n4du'); 
xlabel('Time');
ylabel('Position');
zlabel('Temperature');

% Plot 2:
figure;
plot(x3b, u3b(x3b, t3b(end)), 'or');
hold on;
plot(x3b, U3b(:,end), 'xb');
title('Final Time Error Comparison: e33lo, ylepage, jteeter, n4du');
legend('Actual Solution', 'Approximation');

% Calculate Base Error
error_base = norm(U3b(:,end) - u3b(x3b, t3b(end))) / length(x3b);
fprintf('Base Error (nx=%d, nt=%d): %f\n', nx_base, nt_base, error_base);

% Test 1: Double points in space (nx)
nx_double = 20; 
[x_nx, t_nx, U_nx] = crank_nicolson1d(kappa, x_rng, nx_double, t_rng, nt_base, @u3b_init, @u3b_bndry);
error_nx = norm(U_nx(:,end) - u3b(x_nx, t_nx(end))) / length(x_nx);
fprintf('Error with double nx (nx=%d, nt=%d): %f\n', nx_double, nt_base, error_nx);

% Test 2: Double points in time (nt)
nt_double = 40; 
[x_nt, t_nt, U_nt] = crank_nicolson1d(kappa, x_rng, nx_base, t_rng, nt_double, @u3b_init, @u3b_bndry);
error_nt = norm(U_nt(:,end) - u3b(x_nt, t_nt(end))) / length(x_nt);
fprintf('Error with double nt (nx=%d, nt=%d): %f\n', nx_base, nt_double, error_nt);

% The boundary conditions
function [u] = u3b_bndry(t)
    u = [0*t;
         0*t];
end
