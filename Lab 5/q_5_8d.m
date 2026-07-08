% Lab 5 question 5.8d

% 21x41 Grid of unknowns
U5d = -Inf(21,41);

% Boundary conditions
U5d(1,:) = 100; % Top surface
U5d(end,:) = 0; % Bottom surface
U5d(:, 1) = NaN; % Insulated
U5d(:,end) = NaN; % Insulated

% Circle equation is used to create wires
% didn't want to mess around with index, more elegant like this

% First wire (20V)
for i = 1:21
    for j = 1:41
        if sqrt((i-11)^2 + (j-11)^2) <= 5
            U5d(i,j) = 20;
        end
    end
end

% Second wire (40V)
for i = 1:21
    for j = 1:41
        if sqrt((i-11)^2 + (j-31)^2) <= 5
            U5d(i,j) = 40;
        end
    end
end

U5d_out = laplace2d(U5d);

mesh(U5d_out);
colorbar;
title('ylepage e33lo n4du jteeter');
xlabel('x');
ylabel('y');
zlabel('Voltage (V)');
