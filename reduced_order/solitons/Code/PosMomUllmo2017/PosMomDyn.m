% Author: Kaivalya Bakshi; Date: 27 Mar 2018; Function for computing
% dynamics of soliton position and momentum ODEs of example in figure 1 of 
% Quadratic MFG paper by Ullmo et. al(https://arxiv.org/pdf/1708.07730.pdf) 
% We extend the published example by admitting non trivial dynamics in the 
% SDE dynamics
function dXdt = PosMomDyn(t,X)

global mu par

% dXdt = [X(2,1)/mu; X(1,1)^3 + X(1,1)]; % integrator system
% dXdt = [X(2,1)/mu + par.alpha*(-X(1,1)^3 + X(1,1)); X(1,1)^3 + X(1,1)]; % linear passive drift
dXdt = [X(2,1)/mu + par.alpha*sin(X(1,1)); X(1,1)^3 + X(1,1)]; % trigonometric drift
% dXdt = [X(2,1)/mu + par.alpha*(-X(1,1)^3 + X(1,1)); X(1,1)^3 + X(1,1)]; % bistable potential passive drift
end