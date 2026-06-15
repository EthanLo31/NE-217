% Setup Parameters for the simulation
% Thermal conductivity
kappa = 1.5;  
% Spatial domain
x_range = [0 1];
% Number of spatial points
nx = 6;    
% Time domain
t_range = [0 1];

% Run the simulation with nt = 21
nt1 = 21;
[x3a, t3a, U3a] = crank_nicolson1d(kappa, x_range, nx, t_range, nt1, @u3a_init, @u3a_bndry);

% Plot 1
figure;
mesh(t3a, x3a, U3a); 
title('Heat Distribution: e33lo, ylepage, jteeter, n4du'); 
xlabel('Time');
ylabel('Position');
zlabel('Temperature');

% Run the simulation with nt = 41
nt2 = 41;
[x3a_2, t3a_2, U3a_2] = crank_nicolson1d(kappa, x_range, nx, t_range, nt2, @u3a_init, @u3a_bndry);

% Plot 2
figure;
mesh(t3a_2, x3a_2, U3a_2);
title('Heat Distribution (nt=41): e33lo, ylepage, jteeter, n4du');
xlabel('Time');
ylabel('Position');
zlabel('Temperature');

