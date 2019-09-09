% CONFIDENTIAL (C) Mitsubishi Electric Research Labs (MERL) 2018
% Author: Kaivalya Bakshi; Date: 2 Apr 2018; Function for computing
% initial and final boundary conditions correspoding to soliton ODE 
% dynamics of example in Quadratic MFG paper by Ullmo et. al
% (https://arxiv.org/pdf/1708.07730.pdf), required by MATLAB bvp4c,
% to compute the solution to the BVP of position, momentum ODE dynamics. We
% extend the published example by admitting non trivial dynamics in the SDE
% dynamics
function res = Solitonbvp4cbc(ya,yb)

global mu sigma X0 Sigma0

res = [ ya(1,1) - X0; yb(2,1) + yb(1,1) - 7/2; 
        ya(3,1) - Sigma0;  yb(4,1) - mu*sigma^2 + 2*yb(3,1)^2*1];
end