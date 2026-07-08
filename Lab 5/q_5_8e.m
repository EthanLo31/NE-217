%  Lab 5 question 5.8e

n_acc = 0;

for n = 3:20

    % Unknowns
    U = -Inf(51,51,n);

    % Ground plane
    U(:,:,1) = 0;

    % Insulated side walls
    U(1,:,:)   = NaN;
    U(end,:,:) = NaN;
    U(:,1,:)   = NaN;
    U(:,end,:) = NaN;

    % Insulated top
    U(:,:,end) = NaN;

    % Conducting strip
    U(16:36,16:36,end) = 5;

    % Solve
    U_out = laplace3d(U);

    % Maximum voltage on the flare boundary
    max_volt = max( ...
        [ max(U_out(11,:,:),[],'all'), ...
          max(U_out(41,:,:),[],'all'), ...
          max(U_out(:,11,:),[],'all'), ...
          max(U_out(:,41,:),[],'all') ] );

    fprintf('n = %2d   max edge voltage = %.4f\n', n, max_volt);

    if max_volt <= 1
        n_acc = n;
        U_final = U_out;
    end

end

fprintf('\nLargest acceptable n = %d\n', n_acc);

isosurf1(U_final,51);