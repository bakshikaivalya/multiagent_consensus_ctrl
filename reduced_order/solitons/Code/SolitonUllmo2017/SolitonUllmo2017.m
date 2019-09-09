% CONFIDENTIAL (C) Mitsubishi Electric Research Labs (MERL) 2018
% Author: Kaivalya Bakshi; Date: 2 Apr 2018; Script for computing two
% point boundary value problem (BVP) solution for soliton position and 
% momentum ODEs of example in figure 1 of Quadratic MFG paper by Ullmo et. 
% al (https://arxiv.org/pdf/1708.07730.pdf. We extend the published example
% by admitting non trivial passive drift in the SDE dynamics

clc
clear all
close all

global par mu sigma X0 Sigma0

par.alpha = 1;
par.g = 1;
mu = 1;
sigma = 1;

OPTIONS=optimset('fsolve');
X0 = -3/2;
P0 = 0; % initial guess for momentum at initial time
Sigma0 = 0.4;
Lambda0 = mu*sigma^2; % initial guess for momentum of variance at initial time
t0=0;
dt = 0.001;

%% Computing equilibrium/fixed point of coupled ODEs (Xt*,Pt*,Sigmat*,Lambdat*) 
% Equilibrium point (Xstar,0,Sigmat*,0) for U0(x) = -x^4/4 - x^2/2
k = par.g/(1+par.alpha)/sqrt(1 + par.alpha)/(2*pi)^(par.alpha/2);
Xstar = 0;
fsigfp = @(Sigmastar) Sigmastar^(6 + par.alpha) + Sigmastar^(2 + par.alpha) ...
                    + 2*par.g*par.alpha/(1+par.alpha)...
                        *1/sqrt(1 + par.alpha)/(2*pi)^(par.alpha/2)*Sigmastar^2 ...
                    -(mu^2*sigma^4/2/mu)*Sigmastar^par.alpha;
Sigmastar = fsolve(@(Sigmastar) fsigfp(Sigmastar),Sigma0);
% for alpha = 1
r = roots([1 0 0 1 0 ...
    2*par.g*par.alpha/(1+par.alpha)*1/sqrt(1 + par.alpha)/(2*pi)^(par.alpha/2)...
    -(mu^2*sigma^4/2/mu) 0]);
r = r(r == real(r));
% Sigmastar = sqrt(pi)*mu*sigma^4/par.g; % simpler express par.alpha = 1, uncoupled dynamics case
Sigmastar = r(r > 0);
J = [0 1/mu 0 0; ...
    -3*Sigmastar - 1 0 0 0; ...
    0 0 0 1/2/mu/Sigmastar; ...
    0 0 mu*sigma^4/Sigmastar^3 - k*par.alpha/Sigmastar^(1 + par.alpha) - 2*Sigmastar*(3*Sigmastar^2 + 1) - 6*Sigmastar^3 0];
[v,e] = eig(J);
e = eig(J)

%% For single terminal time
% tf=5;
% P0=fsolve('bvpsm',P0,OPTIONS,'fPosMom','psiPosMom',t0,tf,X0);
% [T,X]=ode45('fPosMom',[0:0.01:tf],[X0;P0]);
% 
% %% Plots
% % Plotting time trajectories of position and momentum
% figure(1)
% plot(T,X)
% title('Pos and Mom time trajectories')
% grid on
% hold on
% 
% figure(2)
% plot(X(:,1),X(:,2))
% title('Pos and Mom phase')
% hold on
% x = min(X(:,1)):0.01:max(X(:,1));
% plot(x,7/2 - x,'--r')
% grid on

%% For multiple terminal times
for tf = 0.01:0.1:15 % for several terminal times
% %% Multiple shooting method for BVP with fsolve. Initial guess for P0 for 
% %tf^(k+1) is the solution P0 for tf^k except for k=1 where provided
% % initial guess is used
% % fPosMomDyn = @ PosMomDyn;
% % ferrtermPosMom = @ errtermPosMom;
% P0=fsolve(@(P0) errODEpropPosMom(t0,tf,X0,P0,@PosMomDyn,@errtermPosMom),P0,OPTIONS);
% t = 0:dt:tf;
% [T,X]=ode45(@PosMomDyn,t,[X0;P0]);
% % Plots
% % Plotting time and phase trajectories of position and momentum
% figure(1)
% plot(T,X)
% title('Pos and Mom time trajectories')
% grid on
% hold on
% 
% figure(2)
% plot(X(:,1),X(:,2))
% title('Pos and Mom phase')
% hold on
% x = min(X(:,1)):0.1:max(X(:,1));
% plot(x,7/2 - x,'--r') % plotting terminal time constraint on 
% grid on
%% Using bvp4c
t = 0:dt:tf;
solinit = bvpinit(t,[X0;P0;Sigma0;Lambda0]);
sol = bvp4c(@SolitonDyn,@Solitonbvp4cbc,solinit);
% Plots
T = sol.x;
X = sol.y;
% Plotting time and phase trajectories of soliton quantities
figure(1)
plot(T,X(1,:),T,X(2,:))
xlabel('Time')
title('Pos and Mom time trajectories')
grid on
hold on

figure(2)
plot(T,X(3,:),T,X(4,:))
xlabel('Time')
title('Var and VarMom time trajectories')
grid on
hold on

figure(3)
plot(X(1,:),X(2,:))
xlabel('X_t')
ylabel('P_t')
title('Pos and Mom phase')
hold on
x = min(X(1,:)):0.1:max(X(1,:));
plot(x,7/2 - x,'--r')
plot(Xstar,0,'*k','linewidth',5)
grid on

figure(4)
plot(X(3,:),X(4,:))
xlabel('\Sigma_t')
ylabel('\Gamma_t')
title('Var and VarMom phase')
hold on
x = min(X(3,:)):0.01:max(X(3,:));
plot(x,mu*sigma^2 - 2*x.^2*1,'--r') 
plot(Sigmastar,0,'*k','linewidth',5)
plot(r,zeros(numel(r),1),'or','linewidth',2)
grid on

figure(5)
options = odeset('RelTol',1e-10);
[tspan,X] = ode45(@SolitonDyn,t,X(:,1),options);
size(X)
E = fEnergy(X');
plot(tspan,E)
xlabel('Time')
ylabel('Energy')
title('Energy trajectories')
grid on
hold on
end