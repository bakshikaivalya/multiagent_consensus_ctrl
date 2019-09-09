% Author: Kaivalya Bakshi; v1 Date 18 May 2017: v mean for 2 particles plotted
% v2 Date 23 May: x,v means computed in sample generator loop; x,v means
% for all particles plotted
% v3 Date 24 May: case for dimensionality of position or velocity, d > 1 natural number included 

clc
clear all
close all
tic; % runtime starts

par.d = 1; % no. of position/velocity dimensions
par.N = 100; % no. of particles should be > 1
par.K = 250; % no. of realizations per particle
par.lambda = 10; % coupling strength SCS
par.beta = 0.8; % coupling constant in comm rate
par.D = 10; % Wiener noise covariance
par.T = 4; % total simulation time
par.dt =  0.01; % time discretization interval length
cr_case = 1; % case 0 and 1 represent unnormalized and normalized comm rate cases

% initial distribution of particle positions and velocities
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

% creating trajectory samples for all particles
V_mean(:,:) = zeros(par.N*par.d, par.T/par.dt+1); % initializing means of each particle V_mean(i,t) = zeros(1,par.T/par.dt)
X_mean(:,:) = zeros(par.N*par.d, par.T/par.dt+1);
for k = 1:1:par.K % iterating over realizations
    X_traj(:,1,k) = X0;
    V_traj(:,1,k) = V0;
    for t = 1:1:par.T/par.dt % time step iteration
        [dX,dV] = csFullMicroDynamics(cr_case,par,X_traj(:,t,k),V_traj(:,t,k)); % one time step update of (x_i,v_i) for all particles i
        X_traj(:,t+1,k) = X_traj(:,t,k) + dX;
        V_traj(:,t+1,k) = V_traj(:,t,k) + dV;
    end
    X_mean(:,:) = X_mean(:,:) + 1/par.K*X_traj(:,:,k); %  compute mean of sampled position trajectories for all particles
    V_mean(:,:) = V_mean(:,:) + 1/par.K*V_traj(:,:,k); %  compute mean of sampled velocity trajectories for all particles
end

% plotting sample trajectories of component 1 <= c <= d of position and velocity for particles # 1 and # 2
% and respective means over all samples of first component for all particles except particle # N, 
% since its initial position and velocity are not visually appealing
c = 1; % component number 1 <= c <= d of position and velocity for particles # 1 and # 2
% for which all sample trajectories are plotted

figure(1) % position trajectories
plot((0:par.dt:par.T),X_mean(1:par.d:((par.N-2)*par.d+1),:),'black','linewidth',1.2) % plot mean first component of state trajectories of all particles except # N
grid on
hold on
for k = 1:1:par.K
    plot((0:par.dt:par.T),X_traj(c:par.d:c+par.d,:,k)) % plot sample k of component c of position trajectories of particles # 1 and # 2
end
plot((0:par.dt:par.T),X_mean(c:par.d:c+par.d,:),'black','linewidth',1.2)
xlabel('Time')
ylabel('position')
title(['Position trajectories, \beta = ' num2str(par.beta)])

figure(2) % velocity trajectories
plot((0:par.dt:par.T),V_mean(1:par.d:((par.N-2)*par.d-1),:),'black','linewidth',1.2) % plot mean first component of |velocity| trajectories of all particles except # N
grid on
hold on
% plot((0:par.dt:par.T),V_mean(1:par.N/2 - 1:par.N/2,:),'black','linewidth',1.2)
for k = 1:1:par.K
    plot((0:par.dt:par.T),V_traj(c:par.d:c+par.d,:,k)) % plot sample k of component c of velocity trajectories for particles # 1 and # 2
end
plot((0:par.dt:par.T),V_mean(c:par.d:c+par.d,:,:),'black','linewidth',1.2)
xlabel('Time')
ylabel('velocity')
title(['Velocity trajectories, \beta = ' num2str(par.beta)])

% plotting cosensus sample trajectories for particles # 1 and # 1 when d = 2
if par.d == 2
    figure(3) % velocity consensus sample trajectories for particle # 1 and # 2 when d = 2
    grid on
    hold on
    plot3(0,0,par.T,'ro','linewidth',2) % consensus (v^x,v^y) = (0,0) at simulation time par.T, which is supposed to be large time
    plot3(V_traj(1,1,1),V_traj(2,1,1),0,'g*','linewidth',2)
    plot3(V_traj(3,1,1),V_traj(4,1,1),0,'g*','linewidth',2)
    for k = 1:1:par.K
    plot3(V_traj(1,:,k),V_traj(2,:,k),(0:par.dt:par.T)) % plot sample k of (v^x,v^y) trajectories of particle # 1
    plot3(V_traj(3,:,k),V_traj(4,:,k),(0:par.dt:par.T)) % plot sample k of (v^x,v^y) trajectories of particle # 2
    end

    figure(4) % state consensus sample trajectories for particle # 1 and # 2 when d = 2
    grid on
    hold on
    plot3(0,0,par.T,'ro','linewidth',2) % consensus (v^x,v^y) = (0,0) at simulation time par.T, which is supposed to be large time
    plot3(X_traj(1,1,1),X_traj(2,1,1),0,'g*','linewidth',2)
    plot3(X_traj(3,1,1),X_traj(4,1,1),0,'g*','linewidth',2)
    for k = 1:1:par.K
    plot3(X_traj(1,:,k),X_traj(2,:,k),(0:par.dt:par.T)) % plot sample k of (v^x,v^y) trajectories of particle # 1
    plot3(X_traj(3,:,k),X_traj(4,:,k),(0:par.dt:par.T)) % plot sample k of (v^x,v^y) trajectories of particle # 2
    end
    
else
    
end

time_min = toc/60 % runtime ends and displayed in min