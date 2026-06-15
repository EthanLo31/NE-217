function [u] = u3a_bndry(t) 

    u = [0*t + -1; 

         0*t +  2]; 

end 

 

function [u] = u3a_init(x) 

    u = 0*x + 1; 

end 

[x3a, t3a, U3a] = crank_nicolson1d( 1.5, [0 1], 6, [0 1], 76, @u3a_init, @u3a_bndry ); 

mesh( t3a, x3a, U3a )  
title( 'e33lo, ylepage, jteeter, n4du' ); 