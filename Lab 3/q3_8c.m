clc

clear

kappa = 1;

x_rng = [0 1];

t_rng = [0 2];

nx = 401;

nt = 4001;

[x3c, t3c, U3c] = crank_nicolson1d(kappa, x_rng, nx, t_rng, nt, @u3c_init, @u3c_bndry);

right_boundary = U3c(end, :);

% find t* when right boundary value is below 0.05.

idx = find(right_boundary <= 0.05, 1, 'first');

t_before = t3c(idx - 1);

t_after = t3c(idx);

u_before = right_boundary(idx - 1);

u_after = right_boundary(idx);

% estimate t*

t1_star = t_before + (0.05 - u_before) .* (t_after - t_before) / (u_after - u_before);

figure;

plot(t3c, right_boundary, 'LineWidth', 1.5);

hold on;

yline(0.05, '--');

xline(t1_star, '--');

xlabel('Time');

ylabel('Temperature at x = 1');

title('n4du e33lo ylepage jteeter');

grid on;

% one side insulated and one side not

function [u] = u3c_init(x)

    u = 0 * x + 1;

end

function [u] = u3c_bndry(t)

    u = [0 * t;

         NaN(size(t))];

end
