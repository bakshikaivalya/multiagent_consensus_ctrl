% CONFIDENTIAL (C) Mitsubishi Electric Research Labs (MERL) 2018
% Author: Kaivalya Bakshi; Date: 27 Mar 2018; Function for computing
% dynamics of soliton position and momentum ODEs of example in figure 1 of 
% Quadratic MFG paper by Ullmo et. al(https://arxiv.org/pdf/1708.07730.pdf) 
% We extend the published example by admitting non trivial dynamics in the 
% SDE dynamics
function dXdt = SolitonDyn(t,X)

global mu sigma par

Xt = X(1,1); Pt = X(2,1); Sigmat = X(3,1); Lambdat = X(4,1); g = par.g; alpha = par.alpha;
c = 1; % c = 0/1 for uncoupled/coupled pos, var dynamics
dXdt = [Pt/mu; ...
        Xt^3 + (1 + c*3*Sigmat^2)*Xt; ...
        Lambdat/2/mu/Sigmat; ...
        2*g*alpha/(1+alpha)*1/sqrt(1 + alpha)/(2*pi)^(alpha/2)*(1/Sigmat)^(alpha) ...
        + (Lambdat^2 - mu^2*sigma^4)/(2*mu*Sigmat^2) ...
        + c*( Sigmat^2*(3*(Xt^2 + Sigmat^2) + 1) )]; % integrator system
    
% dXdt = [X(2,1)/mu + par.alpha*(-X(1,1)^3 + X(1,1)); X(1,1)^3 + X(1,1)]; % linear passive drift
% dXdt = [X(2,1)/mu + par.alpha*sin(X(1,1)); X(1,1)^3 + X(1,1)]; % trigonometric drift
% dXdt = [X(2,1)/mu + par.alpha*(-X(1,1)^3 + X(1,1)); X(1,1)^3 + X(1,1)]; % bistable potential passive drift
end