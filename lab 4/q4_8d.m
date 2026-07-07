clc
clear

n = 11;
[x_out, t_out, U_out] = wave1d(1, [0, 4], n, [0, 12], 3 * n, @u4d_init, @du4d_init, @u4d_bndry);

% plot(x_out, U_out(:, end), 'r.-')
% title(['n=' num2str(n) ' e33lo, ylepage, jteeter, n4du'])

plot(t_out, max(U_out), 'r.-')
title(['n=' num2str(n) ' e33lo, ylepage, jteeter, n4du'])

% local functions

function [u] = du4d_init(x)

    % initial velocity is zero

    u = 0 .* x;

end

function [u] = u4d_init(x)

    u = zeros(size(x));

    % (0,0) to (1,1)

    left = (x <= 1);

    u(left) = x(left);

    % (1,1) to (4,0)

    right = (x > 1);

    u(right) = (4 - x(right)) / 3;

end

function [u] = u4d_bndry(t)

    % boundary conditions: fixed at both ends
    u = [0 * t;
         0 * t];

end
