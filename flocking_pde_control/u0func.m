% Author: Kaivalya Bakshi; v1 Date 29 Jun 2017: script for assigning
% initial condition to nodal points from mesh on the domain
% i/p is the region or mesh nodal points and output is the initial
% condition for all PDEs or dependent variables

function u0 = u0func(region)
%% MFG flocking model parameters

global alpha mu sigma beta r N

nr = length(region.x);
u0 = zeros(N,nr); % initializing initial conditions for both PDE systems or both dependent variables as 0 on the entire domain

%%
s = sqrt(sigma^2*sqrt(r)/2)/alpha;
s1 = sqrt(sigma^2*sqrt(r)/2);

% setting an initial condition perdiodic in the position variable 
% for both density and value function
u0(1,:) = sin(2.*region.x)./sqrt(2).*((region.y - mu).^2./s1.^2 - 1);
u0(2,:) = sin(2.*region.x)./sqrt(2).*((region.y - mu).^2./s1.^2 - 1);

end