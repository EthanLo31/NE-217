[x2a, t2a, U2a] = diffusion1D( 0.1, [0,1], 100, [1, 3], 3916, @u2a_init, @u2a_bndry );
%2.8b error msg

U2aa = U2a

figure;
mesh(x2a,t2a,U2aa)
xlabel('x');
ylabel('y');
zlabel('u(t,x)')
title('n4du...')
colorbar;
