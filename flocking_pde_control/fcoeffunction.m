% Author: Kaivalya Bakshi; v1 Date 27 Jun 2017: function defining "f"
% coefficient which is a vector depending on x,y and partial derivatives 
% of u1 and u2 for the MFG flocking PDE system. Further f2 depends on 
% a functional of u1, namely the CC perturbation. This function returns the
%  "f" function matrix computed on the input 2D domain given by "region"
%  with i/p region and state which contain the unstructured grid node
%  points and u1, u2 solutions at the node points at the present iteration
%  of pdesolve

function [fmatrix] = fcoeffunction(region,state)

%% MFG flocking model parameters
global alpha mu sigma beta r N

nr = numel(region.x);
fmatrix = zeros(N,nr);

%% Assigning f coefficient components
fmatrix(1,:) = -region.y.*state.ux(1,:) - (region.y - mu).*alpha./sqrt(r).*state.uy(1,:) - (region.y - mu).*alpha./(sigma.^2.*r.*sqrt(r)).*state.uy(2,:);
fmatrix(2,:) = region.y.*state.ux(2,:) - (region.y - mu).*alpha./sqrt(r).*state.uy(2,:) + ccpertfunc(region,state,nr,beta,sigma,r,mu,alpha);

end