clc 
clear 
format long 

% n should be sufficiently large at 10^2 

n = 100; 

% approximate the period and nt 

T = 2*pi; 
n_t = 3*n + 1; 

[x4a, t4a, U4a] = wave1d( 1, [0, pi], 10, [0, 3*T], n_t, @u4a_init, @du4a_init, @u4a_bndry ); 

norm0toT = norm( U4a(:, 1) - U4a(:, n + 1) ) 

normTto2T = norm( U4a(:, n + 1) - U4a(:, 2*n + 1) ) 

norm2Tto3T = norm( U4a(:, 2*n + 1) - U4a(:, 3*n + 1) ) 

function [u] = u4a_bndry(t) 

   u = [0*t; 

        0*t]; 

end 

function [u] = u4a_init(x) 

   u = sin(x); 

end 

function [u] = du4a_init(x) 

   u = 0*x; 

end 
