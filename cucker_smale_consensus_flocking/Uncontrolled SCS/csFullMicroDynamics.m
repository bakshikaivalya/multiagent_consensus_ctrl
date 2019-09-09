% Author: Kaivalya Bakshi; Date: 18 May 2017; v1 Date 18 May 2017: v mean for 2 particles plotted
% v3 Date 24 May: case for dimensionality of position or velocity, d > 1 natural number included 

function [dX,dV] = csFullMicroDynamics(cr_case,par,X,V)
% CS dynamics for micro system represented by aggregated [X; V] where X/V =
% [x_i]/[v_i] wherein x_i, v_i are in R^d

% renaming parameters for easier reading of code
d = par.d;
N = par.N;
lambda = par.lambda;
beta = par.beta;
D = par.D;
dt = par.dt;
dW = randn(N*d,1)*sqrt(dt); % Wiener noise

Fv = zeros(N*d,1); % initializing velocity drift vector
cr_sum(:,1) = zeros(N,1); %
for i = 1:1:N
    for j = 1:1:N
        cr_sum(i,1) = cr_sum(i,1) + cr(beta,X(((j - 1)*d + 1):1:j*d),X(((i - 1)*d + 1):1:i*d)); % generating sum of comm rates of particle i w.r.t. all other particles
    end
    switch cr_case % cr_case is 0 or 1 for unnormalized and normalized cases respectively
        case 0
            for j = 1:1:N
        %        Fv(i) = Fv(i) + lambda/N*cr(beta,X(j),X(i))*(V(j) - V(i)) - sqrt(D)/N*dW(j)/dt; % Euler Murayama SDE discretization for d = 1 % insert ((k - 1)*d + 1):1:i*d instead of k (which is either i or j) for multidimensional case
               Fv(((i - 1)*d + 1):1:i*d) = Fv(((i - 1)*d + 1):1:i*d) ...
                                            + lambda/N * cr(beta,X(((j - 1)*d + 1):1:j*d),X(((i - 1)*d + 1):1:i*d)) * (V(((j - 1)*d + 1):1:j*d) - V(((i - 1)*d + 1):1:i*d)) ...
                                            - sqrt(D)/N*dW(((j - 1)*d + 1):1:j*d)/dt; % Euler Murayama SDE discretization
            end
        case 1
            for j = 1:1:N
        %        Fv(i) = Fv(i) + lambda/N*cr(beta,X(j),X(i))*(V(j) - V(i)) - sqrt(D)/N*dW(j)/dt; % Euler Murayama SDE discretization for d = 1 % insert ((k - 1)*d + 1):1:i*d instead of k (which is either i or j) for multidimensional case
               Fv(((i - 1)*d + 1):1:i*d) = Fv(((i - 1)*d + 1):1:i*d) ...
                                            + lambda/cr_sum(i,1) * cr(beta,X(((j - 1)*d + 1):1:j*d),X(((i - 1)*d + 1):1:i*d)) * (V(((j - 1)*d + 1):1:j*d) - V(((i - 1)*d + 1):1:i*d)) ...
                                            - sqrt(D)/N*dW(((j - 1)*d + 1):1:j*d)/dt; % Euler Murayama SDE discretization
            end
    end
end

dX = V*dt;
dV = Fv*dt + sqrt(D)*dW; 
% dZ = [dX;dV]; % state increment as a function of [X; V]

end