%  Lab 5 question 5.8e part 2, plotting for largest acceptable n

n = 9;

U = -Inf(51,51,n);

% Ground plane
U(:,:,1) = 0;

% Insulated boundaries
U(1,:,:)   = NaN;
U(end,:,:) = NaN;
U(:,1,:)   = NaN;
U(:,end,:) = NaN;
U(:,:,end) = NaN;

% Conducting strip (21 x 21 at 5 V)
U(16:36,16:36,end) = 5;

U_out = laplace3d(U);

% Plot
figure;
isosurf1(U_out,51);
title('ylepage e33lo n4du jteeter');
