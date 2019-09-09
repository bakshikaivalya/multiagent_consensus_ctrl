% Author: Kaivalya Bakshi; Date: 27 Mar 2018; Script for computing two
% point boundary value problem (BVP) solution for soliton position and 
% momentum ODEs of example in figure 1 of Quadratic MFG paper by Ullmo et. 
% al (https://arxiv.org/pdf/1708.07730.pdf. We extend the published example
% by admitting non trivial passive drift in the SDE dynamics

clc
clear all
close all

global mu par

par.alpha = 1;
mu = 1;

OPTIONS=optimset('fsolve');
X0 = -3/2;
P0 = 0; % initial guess for momentum at initial time
t0=0;
dt = 0.001;
% tf=5; % for a single terminal time
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
for tf = 0.01:0.1:5 % for several terminal times
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
solinit = bvpinit(t,[X0;P0]);
sol = bvp4c(@PosMomDyn,@PosMombvp4cbc,solinit);
% Plots
T = sol.x;
X = sol.y;
% Plotting time and phase trajectories of position and momentum
figure(1)
plot(T,X)
xlabel('Time')
ylabel('Position')
title('Pos and Mom time trajectories')
grid on
hold on

figure(2)
plot(X(1,:),X(2,:))
xlabel('Position')
ylabel('Momentum')
title('Pos and Mom phase')
hold on
x = min(X(1,:)):0.1:max(X(1,:));
plot(x,7/2 - x,'--r')
grid on
end