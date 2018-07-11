function [h,H]= compute_lm_model(lm)

global XX

dx= lm(1) - XX(1);
dy= lm(2) - XX(2);
d2= dx^2 + dy^2;
d= sqrt(d2);

% calculate h
h= [d; 
    pi_to_pi( atan2(dy,dx) - XX(3) )];

% calculate H
H = [-dx/d, -dy/d,  0;
      dy/d2, -dx/d2, -1];