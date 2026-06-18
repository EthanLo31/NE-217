clc
clear

[x4a, t4a, U4a] = wave1d( 1, [0, pi], 20, [0, 10], 62, @u4a_init, @du4a_init, @u4a_bndry );

mesh( t4a, x4a, U4a )

title( 'n4du e33lo ylepage jteeter' );

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
