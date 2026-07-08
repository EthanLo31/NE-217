% Simple test case from the slides
U = [5 5 5 5 5; NaN -Inf -Inf -Inf 1; NaN -Inf 9 -Inf 1; 1 1 1 1 1]

U_soln = laplace2d(U)


surf(U_soln);
colorbar;
xlabel('x');
ylabel('y');
zlabel('u');
