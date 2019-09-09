% Author: Kaivalya Bakshi; Date: 2 Aug 2017; Function for computing 

function [dX,dV] = ssMFGCtrlSCSFullMicroDynamics(par,V)
%% CS dynamics for micro system represented by aggregated [X; V] where X/V =
% [x_i]/[v_i] wherein x_i, v_i are in R

%% renaming parameters for easier reading of code
d = par.d;
N = par.N;
sigma = par.sigma;
mu = par.mu;
dt = par.dt;
dW = randn(N*d,1)*sqrt(dt); % Wiener noise

%% Dynamics
dX = V*dt;
dV = -(V - ones(N*d,1)*mu)*dt + sigma*dW;

end