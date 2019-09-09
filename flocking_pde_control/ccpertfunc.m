% Author: Kaivalya Bakshi; v1 Date 28 Jun 2017: function defining the CC
% perturbation which is a scalar depending on x,y and the function u1 on
% the domain defined by region. This function returns the CC perturbation
% vector giving the CC perturbation values at all node points in "region"
% with i/p region and state which contain the unstructured grid node
% points and u1, u2 solutions at the node points at the present iteration
% of pdesolve

function [ccpertmatrix] = ccpertfunc(region,state,nr,beta,sigma,r,mu,alpha)

global lengthx widthv

s1 = sqrt(sigma^2*sqrt(r)/2);
ccpertmatrix = zeros(1,nr);

if numel(state.u(1,:)) ~= 1
    X = region.x; Y = region.y; P = [X; Y]; % extract the pdesolve mesh specified by user
    V = state.u(1,:);
%         size(P) % DEBUGGING
%         size(V) % DEBUGGING
    FV = scatteredInterpolant(P',V','linear'); % create interpolant using scattered data or data over unstructured mesh P
    dx = 0.05; dy = 0.05; [Xq,Yq] = meshgrid(-lengthx:dx:+lengthx,-widthv:dy:+widthv); Q = [Xq(:),Yq(:)]; % creating mesh Q used for intergration
%         size(Q) % DEBUGGING
    Vq = FV(Q); % compute interpolant on the required mesh Q
%         zeros(numel(Xq(:)),1) % DEBUGGING
    I = zeros(numel(Xq(:)),1);
%         size(I); % DEBUGGING
    for k = 1:1:numel(Xq(:)) % compute CC perturbation integral at all mesh grid points
%         k; % DEBUGGING
        for l = 1:1:numel(Xq(:)) % sum over all mesh grid points to compute CC perturbation integral at mesh gridpoint indexed k
%               Vq(l); % DEBUGGING
            I(k,1) = I(k,1) + ( 1 + (Q(l,1) - Q(k,1))^2 )^(-beta) * exp((Q(l,2) - mu)^2/2/s1^2)/sqrt(2*pi*s1^2) * ( Q(l,2) - Q(k,2) )* Vq(l);
        end
        I(k,1) = 2*(Q(k,2) - mu)/alpha*I(k,1)*dx*dy;        
    end      
    FI = scatteredInterpolant(Q,I,'linear'); % create interpolant for CC perturbation integral using data over required mesh Q
    ccpertmatrix(1,:) = FI(P'); % compute CC perturbation interpolant over the mesh P
%         size(state.u) % DEBUGGING       
end

end