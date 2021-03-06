% Author: Kaivalya Bakshi; Date: 24 Apr 2018; Function for computing
% dynamics of soliton position and momentum ODEs for Hamiltonian system  
% described in MERL report based on Quadratic MFG paper by Ullmo et. al
% (https://arxiv.org/pdf/1708.07730.pdf). We extend the published example 
% by admitting coupled and non trivial passive drift in the SDE dynamics
function dXdt = SolitonDyn(t,X)

global mu sigma par c % c = 0/1 for uncoupled/coupled pos, var dynamics

q1 = X(1,1); p1 = X(2,1); q2 = X(3,1); p2 = X(4,1); g = par.g; alpha = par.alpha;
k = g/(1+alpha)/sqrt(1 + alpha)/(2*pi)^(alpha/2);

% 1D SDE dynamics with U0 = -x^4/4 - x^2/2
dXdt = [p1/mu; ...
        q1^3 + (1 + c*3*q2^2)*q1; ...
        p2/mu; ...
        -mu*sigma^4/(4*q2^3) + k*alpha/q2^(alpha+1) ...
        - c*q2*(-3*(q1^2 + q2^2) - 1)]; % integrator system
    
% dXdt = [X(2,1)/mu + par.alpha*(-X(1,1)^3 + X(1,1)); X(1,1)^3 + X(1,1)]; % linear passive drift
% dXdt = [X(2,1)/mu + par.alpha*sin(X(1,1)); X(1,1)^3 + X(1,1)]; % trigonometric drift
% dXdt = [X(2,1)/mu + par.alpha*(-X(1,1)^3 + X(1,1)); X(1,1)^3 + X(1,1)]; % bistable potential passive drift
end