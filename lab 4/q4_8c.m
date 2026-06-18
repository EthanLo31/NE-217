clc 
clear 

% n should be sufficiently large at 10^2 

n = 100; 

% approximate the period and nt 

T = 2*pi; 
n_t = 3*n + 1; 

[x4a, t4a, U4a] = wave1d( 1, [0, pi], 10, [0, 3*T], n_t, @u4a_init, @du4a_init, @u4c_bndry ); 
 

mesh( t4a, x4a, U4a ) 

title( 'n4du e33lo ylepage jteeter' ) 

xlabel( 't' ) 

ylabel( 'x' ) 

zlabel( 'u(x,t)' )

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

function [u] = u4c_bndry(t)
    u = [0*t;
         NaN*t];
end
