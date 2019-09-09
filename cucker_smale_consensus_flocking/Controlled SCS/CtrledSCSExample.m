% Author: Kaivalya Bakshi; Date: 2 Aug 2017; Script for plotting steady state
% feedback controlled (x,v) in R^2 trajectories at equilibrium/fixed point
% of Stochastic Cucker Smale MFG PDE system by Nourian, Caines, Malhame as
% presented in IFAC 2011 paper

clc
clear all
close all

par.d = 1; % no. of position/velocity dimensions
par.N = 100; % no. of particles should be > 1
par.K = 250; % no. of realizations per particle
par.sigma = 2; % Wiener noise covariance
par.mu = 0; % initial mean which is the desired state consensus
par.T = 5; % total simulation time
par.dt =  0.01; % time discretization interval length

% renaming parameters for easier reading of code
d = par.d;
N = par.N;
K = par.K; % no. of realizations per particle
sigma = par.sigma;
mu = par.mu;
dt = par.dt;
dW = randn(N*d,1)*sqrt(dt);

%% Initial distribution of particle positions and velocities
% with (component of x_i/v_i) ~ uni[-50,+50] for 1 <= i <= N - 1
% and the x_N/v_N satisfying sum_i (component of x_i) = 0 or sum_i (component of v_i) = 0
    % case d = 1
% X0 = -50*ones(par.N-1,1) + 100*rand(par.N-1,1);
% X0(par.N,1) = -sum(X0(1:1:par.N-1));
% V0 = -50*ones(par.N-1,1) + 100*rand(par.N-1,1);
% V0(par.N,1) = -sum(V0(1:1:par.N-1));
    % case d is a natural number
X0(1:1:(par.N-1)*par.d,1) = -50*ones((par.N-1)*par.d,1) + 100*rand((par.N-1)*par.d,1);
X0(((par.N-1)*par.d+1):1:par.N*par.d,1) = zeros(par.d,1);
for j = 1:1:(par.N-1)
    X0(((par.N-1)*par.d+1):1:par.N*par.d,1) = X0(((par.N-1)*par.d+1):1:par.N*par.d,1) - X0(((j-1)*par.d+1):1:j*par.d);
end
V0(1:1:(par.N-1)*par.d,1) = -50*ones((par.N-1)*par.d,1) + 100*rand((par.N-1)*par.d,1);
V0(((par.N-1)*par.d+1):1:par.N*par.d,1) = zeros(par.d,1);
for j = 1:1:(par.N-1)
    V0(((par.N-1)*par.d+1):1:par.N*par.d,1) = V0(((par.N-1)*par.d+1):1:par.N*par.d,1) - V0(((j-1)*par.d+1):1:par.N*par.d:1:j*par.d);
end

%% Creating trajectory samples for all particles
X_traj(:,1,:) = zeros(par.N*par.d, K); % initializing means of each particle V_mean(i,t) = zeros(1,par.T/par.dt)
V_traj(:,1,:) = zeros(par.N*par.d, K); % initializing means of each particle V_mean(i,t) = zeros(1,par.T/par.dt)
V_mean(:,:) = zeros(par.N*par.d, par.T/par.dt+1); % initializing means of each particle V_mean(i,t) = zeros(1,par.T/par.dt)
X_mean(:,:) = zeros(par.N*par.d, par.T/par.dt+1);
for k = 1:1:K % iterating over realizations
    X_traj(:,1,k) = X0;
    V_traj(:,1,k) = V0;
    for t = 1:1:par.T/par.dt % time step iteration
        [dX,dV] = ssMFGCtrlSCSFullMicroDynamics(par,V_traj(:,t,k)); % one time step update of (x_i,v_i) for all particles i
        X_traj(:,t+1,k) = X_traj(:,t,k) + dX;
        V_traj(:,t+1,k) = V_traj(:,t,k) + dV;
    end
    V_mean(:,:) = V_mean(:,:) + 1/par.K*V_traj(:,:,k); %  compute mean of sampled velocity trajectories for all particles
    X_mean(:,:) = X_mean(:,:) + 1/par.K*X_traj(:,:,k); %  compute mean of sampled position trajectories for all particles
end

%% Plotting sample trajectories of component 1 <= c <= d of position and velocity for particles # 1 and # 2
% and respective means over all samples of first component for all particles except particle # N, 
% since its initial position and velocity are not visually appealing
c = 1; % component number 1 <= c <= d of position and velocity for particles # 1 and # 2
% for which all sample trajectories are plotted

figure(1) % position trajectories
plot((0:par.dt:par.T),X_mean(1:par.d:((par.N-2)*par.d+1),:),'linewidth',1.2) % plot mean first component of state trajectories of all particles except # N
grid on
hold on
for k = 1:1:par.K
    plot((0:par.dt:par.T),X_traj(c:par.d:c+par.d,:,k)) % plot sample k of component c of position trajectories of particles # 1 and # 2
end
plot((0:par.dt:par.T),X_mean(c:par.d:c+par.d,:),'black','linewidth',1.2)
xlabel('Time')
ylabel('position')
title(['Position trajectories, \mu = ' num2str(mu)])

figure(2) % velocity trajectories
plot((0:par.dt:par.T),V_mean(1:par.d:((par.N-2)*par.d-1),:),'linewidth',1.2) % plot mean first component of |velocity| trajectories of all particles except # N
grid on
hold on
for k = 1:1:par.K
    plot((0:par.dt:par.T),V_traj(c:par.d:c+par.d,:,k)) % plot sample k of component c of velocity trajectories for particles # 1 and # 2
end
plot((0:par.dt:par.T),V_mean(c:par.d:c+par.d,:),'black','linewidth',1.2)
xlabel('Time')
ylabel('velocity')
title(['Velocity trajectories, \mu = ' num2str(mu)'])