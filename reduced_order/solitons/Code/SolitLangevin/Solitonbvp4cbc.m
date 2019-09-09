% CONFIDENTIAL (C) Mitsubishi Electric Research Labs (MERL) 2018
% Author: Kaivalya Bakshi; Date: 12 Apr 2018; Function for computing
% initial and final boundary conditions correspoding to soliton ODE 
% dynamics of Hamiltonian system described in MERL report based on 
% Quadratic MFG paper by Ullmo et. al (https://arxiv.org/pdf/1708.07730.pdf), 
% required by MATLAB bvp4c, to compute the solution to the BVP of position,
%  momentum ODE dynamics. We extend the published example by admitting 
% coupled and non trivial passive drift in the SDE dynamics
function res = Solitonbvp4cbc(ya,yb)

global mu sigma q10 q20 Sigmastar

pterm = 2*yb(3,1)*yb(4,1) - mu*sigma^2 + 2*yb(3,1)^2*1;
% pterm = yb(3,1) - mu*sigma^2;

res = [ ya(1,1) - q10; yb(2,1) + yb(1,1) - 7/2; 
        ya(3,1) - q20;  pterm];
end