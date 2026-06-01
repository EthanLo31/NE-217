function [u] = u2a_bndry(t)
    u = [0*t + 1.7;
         0*t + 4.1];
end

function [u] = u2a_init(x)
    u = 0*x + 0.9;
end

[x1 t1 U1] = diffusion1D(0.1, [0,1], 100, [1,3], 4000, @u2a_init, @u2a_bndry);
mesh(t1, x1, U1)
title( 'ylepage and e33lo and n4du and jteeter' );