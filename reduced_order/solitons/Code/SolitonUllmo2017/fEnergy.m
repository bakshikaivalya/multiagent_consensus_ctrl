% CONFIDENTIAL (C) Mitsubishi Electric Research Labs (MERL) 2018
% Author: Kaivalya Bakshi; Date: 12 Apr 2018; Function for computing
% energy of solitons introduced in Quadratic MFG paper by Ullmo et. al
% (https://arxiv.org/pdf/1708.07730.pdf). We extend the published example 
% by admitting coupled and further non trivial passive drift in the SDE 
% dynamics

function E = fEnergy(X)

global mu sigma par

Xt = X(1,:); Pt = X(2,:); Sigmat = X(3,:); Lambdat = X(4,:); g = par.g; alpha = par.alpha;
k = par.g/(1+par.alpha)/sqrt(1 + par.alpha)/(2*pi)^(par.alpha/2);
c = 1; % c = 0/1 for uncoupled/coupled pos, var dynamics

% Total energy of soliton with 1D SDE dynamics and U0 = -x^4/4 - x^2/2
E = Pt.^2./(2*mu) + c*(Lambdat.^2 - mu^2*sigma^4)./(8.*mu.*Sigmat.^2) ...
    + c*k./Sigmat.^alpha ...
    - (Xt.^4 + 6*Sigmat.^2.*Xt.^2 + 3.*Sigmat.^4)./4 - (Xt.^2 + Sigmat.^2)./2; 

end