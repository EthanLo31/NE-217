% Lab 5 question 5.8a
U5a = zeros( 6, 9 );
U5a(1:3,1:3) = 5;
U5a(3:6,8:9) = -5;
U5a(2:3,2:3) = -Inf;
U5a(3:4,3:4) = -Inf;
U5a(4:5,4:8) = -Inf;

U5a_out = laplace2d(U5a)