clc 
clear 

 

% (Main Script goes here) 

 

% local functions  

function [u] = du4d_init(x) 

% initial velocity is zero 

u = 0*x; 

end 

function [u] = u4d_init(x) 

u = zeros(size(x)); 

% (0,0) to (1,1) 

left = (x <= 1); 

u(left) = x(left); 

% (1,1) to (4,0) 

right = (x > 1); 

u(right) = (4 - x(right))/3; 

end 
