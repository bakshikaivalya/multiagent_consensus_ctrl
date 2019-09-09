clc
clear all
close all

c = 1; % c = 0/1 for uncoupled/coupled pos, var dynamics
par.alpha = 3;
par.g = 1;
mu = 1;
sigma = 1;

g = par.g;
alpha = par.alpha;

k = g/(1 + alpha)/sqrt(1 + alpha)/(2*pi)^(alpha/2);

% Computing fixed point and eigenproperties in case c=0/1 and alpha <= 2
if c == 0
    fprintf('UNCOUPLED SOLITON DYNAMICS')
else
    fprintf('COUPLED SOLITON DYNAMICS')
end
if c == 0
    qstaranalyticc0 = (mu*sigma^4/(4*k*alpha))^(1/(2-alpha)) % for c = 0, display fixed point value
    qstar = qstaranalyticc0;
end
syms qstar
p = -1/4*mu*sigma^4*qstar^(alpha + 1) ...
                        + k*alpha*qstar^3 + c*( qstar^(alpha + 5) ...
                        + 3*qstar^(alpha + 7) );
qstar = roots(sym2poly(p));
qstar = qstar(imag(qstar) == 0);
qstar = qstar(qstar>0) % display all real, positive fixed points
if alpha <= 2
    fprintf('minimal equilibrium value')
    qstar = min(qstar) % display minimal equilibrium value
    Jacobian = [0, 1/mu, 0, 0; ...
        3*qstar + 1, 0, 0, 0; ...
        0, 0, 0, 1/mu; ...
        0, 0, ...
        +3/4*mu*sigma^4/qstar^4 ...
        - k*alpha*(1 + alpha)/qstar^(2 + alpha) ...
        + c*((3*qstar^2 + 1) + 6*qstar^2), 0];
    if alpha == 1
       Jacobian 
    end
    [v,e] = eig(Jacobian)
end
% Computing multiple fixed points and eigenproperties for case c=1, alpha > 2
if c == 1 && alpha > 2 
    % unstable fixed point
    fprintf('UNSTABLE EQUILIBRIUM POINT')
    qstar1 = min(qstar)
    Jacobian1 = [0, 1/mu, 0, 0; ...
        3*qstar1 + 1, 0, 0, 0; ...
        0, 0, 0, 1/mu; ...
        0, 0, ...
        +3/4*mu*sigma^4/qstar1^4 ...
        - k*alpha*(1 + alpha)/qstar1^(2 + alpha) ...
        + c*((3*qstar1^2 + 1) + 6*qstar1^2), 0];
    [v1,e1] = eig(Jacobian1)
    % stable fixed point    
    fprintf('STABLE EQUILIBRIUM POINT')
    qstar2 = max(qstar)
    Jacobian2 = [0, 1/mu, 0, 0; ...
        3*qstar2 + 1, 0, 0, 0; ...
        0, 0, 0, 1/mu; ...
        0, 0, ...
        +3/4*mu*sigma^4/qstar2^4 ...
        - k*alpha*(1 + alpha)/qstar2^(2 + alpha) ...
        + c*((3*qstar2^2 + 1) + 6*qstar2^2), 0];
    [v2,e2] = eig(Jacobian2)
end

% Plotting zero momentum energy contours
q1 = linspace(-2,2,20);
q2 = linspace(0,2,30);
[q1,q2] = meshgrid(q1,q2);
Z = - mu*sigma^4./(8.*q1.^2) ...
        + k./q2.^alpha - (q1.^4 + 6*q1.^2.*q2.^2 + 3.*q2.^4)./4 ...
        - (q1.^2 + q2.^2)./2;
figure
grid on
plot(zeros(numel(qstar)),qstar,'r*')
hold on
contour(q1,q2,Z,'ShowText','on')
xlabel('q_1')
ylabel('q_2')