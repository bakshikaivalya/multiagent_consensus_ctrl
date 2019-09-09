% CONFIDENTIAL (C) Mitsubishi Electric Research Labs (MERL) 2018
% Author: Kaivalya Bakshi; Date: 24 Apr 2018; Function for computing
% total energy of solitons with Hamiltonian system described in MERL report
% based on Quadratic MFG paper by Ullmo et. al (https://arxiv.org/pdf/1708.07730.pdf). 
% We extend the published example by admitting coupled and further non 
% trivial passive drift in the SDE dynamics

function E = fEnergy(X)

global mu sigma par Sigmastar c % c = 0/1 for uncoupled/coupled pos, var dynamics

q1 = X(1,:); p1 = X(2,:); q2 = X(3,:); p2 = X(4,:); g = par.g; alpha = par.alpha;
k = g/(1+alpha)/sqrt(1 + alpha)/(2*pi)^(alpha/2);

% Total energy of soliton with 1D SDE dynamics and U0 = -x^4/4 - x^2/2
if c == 1
    E = p1.^2./(2*mu) + p2.^2./(2*mu) - mu*sigma^4./(8.*q2.^2) ...
        + k./(q2.^alpha) - (q1.^4 + 6*q2.^2.*q1.^2 + 3.*q2.^4)./4 ...
        - (q1.^2 + q2.^2)./2;
else
   E(1,:) = p1.^2./(2*mu) + (-q1.^4./4 - q1.^2./2);
   E(2,:) = p2.^2./(2*mu) ...
            - mu*sigma^4./(8.*q2.^2) + k./(q2.^alpha);
end

end