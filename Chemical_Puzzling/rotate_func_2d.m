function [sse]=rotate_func_2d(x,recreate,true_pos_trans)

theta=x;

R=[cos(theta) -sin(theta); sin(theta) cos(theta)];

diff=recreate*R-true_pos_trans;    
sse=norm(diff);