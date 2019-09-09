% Author: Kaivalya Bakshi; Date: 24 Apr 2018; Script for computing two
% point boundary value problem (BVP) solution for soliton position and 
% momentum ODEs of example in figure 1 of Quadratic MFG paper by Ullmo et. 
% al (https://arxiv.org/pdf/1708.07730.pdf) with external potential 
% U0(x) = -x^4/4 - x^2/2. We extend the published example by admitting 
% coupled and non trivial passive drift in the SDE dynamics

clc
clear all
close all

global par mu sigma q10 q20 c

c = 1 % c = 0/1 for uncoupled/coupled pos, var dynamics
par.alpha = 3;
par.g = 1;
mu = 1;
sigma = 1;

g = par.g;
alpha = par.alpha

OPTIONS=optimset('fsolve');
q10 = -3/2;
p10 = 0; % initial guess for momentum at initial time
q20 = 0.5;
p20 = mu*sigma^2; % initial guess for momentum of variance at initial time
t0=0;
dt = 0.001;

%% Computing equilibrium/fixed point of coupled ODEs (q1*,p2*,q2*,p2*) 
% Equilibrium point (0,0,q2*,0) for U0(x) = -x^4/4 - x^2/2
k = g/(1 + alpha)/sqrt(1 + alpha)/(2*pi)^(alpha/2);
% Computing fixed point and eigenproperties in case c=0/1 and alpha <= 2
if c == 0
    fprintf('UNCOUPLED SOLITON DYNAMICS')
else
    fprintf('COUPLED SOLITON DYNAMICS')
end
if c == 0
    if alpha ~= 2
        qstaranalyticc0 = (mu*sigma^4/(4*k*alpha))^(1/(2-alpha)) % for c = 0, display fixed point value
        qstar = qstaranalyticc0;
    else 
        qstaranalyticc0 = 0 % for c = 0, display fixed point value
    end
end
syms q2star
p = -1/4*mu*sigma^4*q2star^(alpha + 1) ...
                        + k*alpha*q2star^3 + c*( q2star^(alpha + 5) ...
                        + 3*q2star^(alpha + 7) );
q2star = roots(sym2poly(p));
q2star = q2star(imag(q2star) == 0);
q2star = q2star(q2star>0) % display all real, positive fixed points
if alpha <= 2
    fprintf('minimal equilibrium value')
    q2star = min(q2star) % display minimal equilibrium value
    Jacobian = [0, 1/mu, 0, 0; ...
        3*q2star + 1, 0, 0, 0; ...
        0, 0, 0, 1/mu; ...
        0, 0, ...
        +3/4*mu*sigma^4/q2star^4 ...
        - k*alpha*(1 + alpha)/q2star^(2 + alpha) ...
        + c*((3*q2star^2 + 1) + 6*q2star^2), 0];
    if alpha == 1
       Jacobian 
    end
    [v,e] = eig(Jacobian)
end
% Computing multiple fixed points and eigenproperties for case c=1, alpha > 2
if c == 1 && alpha > 2 
    % unstable fixed point
    fprintf('EQUILIBRIUM VARIANCE 1')
    q2star1 = min(q2star)
    Jacobian1 = [0, 1/mu, 0, 0; ...
        3*q2star1 + 1, 0, 0, 0; ...
        0, 0, 0, 1/mu; ...
        0, 0, ...
        +3/4*mu*sigma^4/q2star1^4 ...
        - k*alpha*(1 + alpha)/q2star1^(2 + alpha) ...
        + c*((3*q2star1^2 + 1) + 6*q2star1^2), 0];
    [v1,e1] = eig(Jacobian1)
    % stable fixed point    
    fprintf('EQUILIBRIUM VARIANCE 2')
    q2star2 = max(q2star)
    Jacobian2 = [0, 1/mu, 0, 0; ...
        3*q2star2 + 1, 0, 0, 0; ...
        0, 0, 0, 1/mu; ...
        0, 0, ...
        +3/4*mu*sigma^4/q2star2^4 ...
        - k*alpha*(1 + alpha)/q2star2^(2 + alpha) ...
        + c*((3*q2star2^2 + 1) + 6*q2star2^2), 0];
    [v2,e2] = eig(Jacobian2)
end

% Plotting zero momentum energy contours
q1 = linspace(-2,2,20);
q2 = linspace(0,2,30);
[q1,q2] = meshgrid(q1,q2);
Z = - mu*sigma^4./(8.*q1.^2) ...
        + k./q2.^alpha - (q1.^4 + 6*q1.^2.*q2.^2 + 3.*q2.^4)./4 ...
        - (q1.^2 + q2.^2)./2;
figure(1)
grid on
plot(zeros(numel(q2star)),q2star,'r*')
hold on
contour(q1,q2,Z,'ShowText','on')
xlabel('q_1')
ylabel('q_2')
title(['Zero momentum energy contours, c = ', num2str(c), ', \alpha = ', num2str(alpha)])

%% Solving soliton dynamics BVPs using bvp4c for various terminal time values
i = 1; % index for storing boundary condition residuals
for tf = 0.01:0.5:10 % several terminal time values
t = 0:dt:tf;
solinit = bvpinit(t,[0;0;3;0]);
sol = bvp4c(@SolitonDyn,@Solitonbvp4cbc,solinit);

% Plots
T = sol.x;
X = sol.y;
% Plotting time and phase trajectories of soliton quantities
figure(2)
plot(T,X(1,:),T,X(2,:))
xlabel('Time')
title(['Pos and Mom time trajectories, c = ', num2str(c), ', \alpha = ', num2str(alpha)])
grid on
hold on

figure(3)
plot(T,X(3,:),T,X(4,:))
xlabel('Time')
title(['Var and VarMom time trajectories, c = ', num2str(c), ', \alpha = ', num2str(alpha)])
grid on
hold on

figure(4)
plot(X(1,:),X(2,:))
xlabel('q_1')
ylabel('p_1')
title(['Pos and Mom phase, c = ', num2str(c), ', \alpha = ', num2str(alpha)])
hold on
x = min(X(1,:)):0.1:max(X(1,:));
plot(x,7/2 - x,'--r')
plot(0,0,'*k','linewidth',5)
grid on

figure(5)
plot(X(3,:),X(4,:))
xlabel('q_2')
ylabel('p_2')
title(['Canonical Var and VarMom phase, c = ', num2str(c), ', \alpha = ', num2str(alpha)])
hold on
x = min(X(3,:)):0.01:max(X(3,:));
plot(x,mu*sigma^2./2./x - x*1,'--r') 
plot(q2star,zeros(numel(q2star)),'*k','linewidth',5)
grid on

figure(6)
grid on
hold on
E = fEnergy(X);
if c== 1
    plot(T,E)
else
    plot(T,E(1:2,:))
end
xlabel('Time')
ylabel('Energy')
title(['Energy trajectories for various terminal time values, c = ', num2str(c), ', \alpha = ', num2str(alpha)])

figure(7)
grid on
hold on
xlabel('Time')
ylabel('Residuals')
title(['Residuals of boundary conditions for various terminal time values, c = ', num2str(c), ', \alpha = ', num2str(alpha)])
ya =  [q10;  X(2,1); q20 ; X(4,1)];
yb = [X(1,end); X(2,end); X(3,end); X(4,end)];
res(i,:) = Solitonbvp4cbc(ya,yb); % computing boundary condition residuals
plot(T(end),res(i,:),'r.')
i = 1+1;
end