function [sse]=rotate_func_2d(theta,rec_pos_trans,true_pos_trans)

%This function rotates the reconstructed position and gives the sum squared error
%between this rotation and the true position. It works for 2 dimensional
%inputs.

%The inputs are angle (theta), centered reconstructed positions (rec_pos_trans), and
%centered true positions (true_pos_trans). The output is the sum squared
%error (sse)

R=[cos(theta) -sin(theta); sin(theta) cos(theta)];

diff=rec_pos_trans*R-true_pos_trans;    
sse=norm(diff);