% Lab 5 question 5.8c

% 81x81 grid of unknowns
U = -Inf(81,81);

% Boundary conditions
U(:,1) = 100; % Hot wall
U(1,:) = 0; % Bottom wall
U(end, :) = 0; % Top wall
U(:,end) = NaN; % Insulated

U5c_out = laplace2d(U);

% Plot
mesh(U5c_out);
colorbar;
title('ylepage e33lo n4du jteeter');
xlabel('x');
ylabel('y');
zlabel('Temperature (C)');

% Temperature of the middle of the insulated wall
% U5c_out(41,81) remains NaN since it's a bound, but since no heat flow,
% value will be pretty much identical to the one at end-1 (80)
mid_temp = U5c_out(41,80)