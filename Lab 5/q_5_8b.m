% Lab 5 question 5.8b
U5b = NaN*ones( 6, 9 );
U5b(1:3,1:3) = 5;
U5b(3:6,8:9) = -5;
U5b(2:3,2:3) = -Inf;
U5b(3:4,3:4) = -Inf;
U5b(4:5,4:8) = -Inf;

U5b_out = laplace2d(U5b)