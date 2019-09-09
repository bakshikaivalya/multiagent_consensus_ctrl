% Author: Kaivalya Bakshi; Date: 27 Mar 2018; Function for computing
% initial and final boundary conditions correspoding to position and 
% momentum ODEs of example in figure 1 of Quadratic MFG paper by Ullmo et.
% al (https://arxiv.org/pdf/1708.07730.pdf), required by MATLAB bvp4c,
% to compute the solution to the BVP of position, momentum ODE dynamics. We
% extend the published example by admitting non trivial dynamics in the SDE
% dynamics
function res = PosMombvp4cbc(ya,yb)
res = [ ya(1,1) + 3/2; yb(2,1) + yb(1,1) - 7/2 ];
end