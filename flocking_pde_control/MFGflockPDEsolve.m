% Author: Kaivalya Bakshi; v1 Date 27 Jun 2017: script for computing
% solution to MFG flocking model perturbation coupled PDE system with
% normalized cost coupling

% The PDE system consists of forward FP and backward HJB PDEs which are
% coupled by the cost coupling (CC) running cost. The PDEs are not directly
% available in the divergence form expected by the MATLAB PDE toolbox
% pdesolve and are massaged to extract that form before coding up here. We
% are therefore creating a model and solver for two time dependent PDE over
% two independent variables. Therefore the domain over which we solve the
% PDEs is the same for both dependent variables or both PDEs

%%
clc
clear all
close all

tic
%% MFG flocking model parameters
par.mu = 0; % velocity consensus mean
par.sigma = 0.2; % Wiener noise variance
par.r = 1; % control cost
par.beta = 1.2; % communication rate constant for CS dynamics
par.alpha = integral(@(x) 1./(1 + x.^2).^par.beta,-Inf,Inf); % state cost normalizing term

par.N = 2; % number of PDEs

global alpha mu sigma beta r N lengthx widthv

alpha = par.alpha;
mu = par.mu;
sigma = par.sigma;
beta = par.beta;
r = par.r;

N = par.N;

%% PDE Model with a two dependent variables
numberOfPDE = 2;
model = createpde(numberOfPDE);

%% Geometry 
% Defining a rectangular geometry for domain of solution of PDE system
lengthx = pi/2;
widthv = 0.5;

% Defining the rectangle by giving the 4 x-locations followed by the 4
% y-locations of the corners using construct solid geometry command
gdm = [3 4 -lengthx +lengthx +lengthx -lengthx -widthv -widthv widthv widthv]';
g = decsg(gdm, 'S1', ('S1')');

% Converting the DECSG geometry into a geometry object
% on doing so it is appended to the PDEModel
geometryFromEdges(model,g);

% Plot the geometry and display the edge labels for use in the boundary
% condition definition.
figure; 
pdegplot(model,'EdgeLabels','on'); 
axis([-lengthx-.1 +lengthx+.1 -widthv-.1 +widthv+.1]);
title 'Geometry With Edge Labels Displayed';

figure; 
pdegplot(model,'FaceLabels','on')
axis([-lengthx-.1 +lengthx+.1 -widthv-.1 +widthv+.1]);
title 'Geometry With Edge Labels Displayed';

%% Mesh
% Create the triangular mesh on the square with approximately
% ten elements in each direction.
hmax = .1; % element size
msh = generateMesh(model,'Hmax',hmax);
figure; 
pdeplot(model); 
axis equal
title 'Rectangle With Triangular Element Mesh'
xlabel 'X-coordinate, meters'
ylabel 'Y-coordinate, meters'
 
% %% Steady State Solution %%
% 
% %% Definition of PDE Coefficients
% % The expressions for the coefficients required by PDE Toolbox can easily
% % be identified by comparing the equation above with the scalar parabolic
% % equation in the PDE Toolbox documentation.
% 
% % The "m", "a" coefficients are zero, "d" is an identity matrix and the
% % "c" coefficient is a constant matrix in the PDE system considered and 
% % "f" coefficient is a vector depending on partial derivatives of u1 and 
% % u2. Further f2 depends on a functional of u1, namely the CC perturbation.
% % We define f by MATLAB function fcoeffunction whichin the MATLAB function 
% % ccpertfunc will define the CC perturbation function. In the steady state
% % case "d" is zero as well.
% m = zeros(4,1); % zeros(2,2); % 
% d = eye(2); %zeros(4,1); % 
% c = [0;0;0;sigma^2/2; 0;0;0;0; 0;0;0;1/2/r; 0;0;0;sigma^2/2];
% a = [0;0;0;0]; % zeros(2,2); % 
% f = @fcoeffunction;
% 
% specifyCoefficients(model,'d',0,'m',0,'c',c,'a',0,'f',f);
% 
% %% Initial and Boundary Conditions
% % Boundary conditions are assumed to be both the density and value function
% % to be zero on the boundary since both these functions are assumed to tend
% % to zero as the radius tends to infinity. Thus we only need to set the
% % Dirichilet boundary condition in this case and not Neumann conditions
% 
% % applyBoundaryCondition(model,'dirichlet','Edge',1,'u',[0; 0]);
% % applyBoundaryCondition(model,'dirichlet','Edge',2,'u',[0; 0]);
% % applyBoundaryCondition(model,'dirichlet','Edge',3,'u',[0; 0]);
% % applyBoundaryCondition(model,'dirichlet','Edge',4,'u',[0; 0]);
% applyBoundaryCondition(model,'dirichlet','edge',1:4,'h',[1,0;0,1],'r',[0;0]);
% % bc2 = applyBoundaryCondition(model,'neumann','Edge',1:4,'g',[0;0], 'q',zeros(2,2));
% 
% % Initial conditions for PDE system. Initial perturbations should satisfy
% % certain properties, namely belonging to the class L2(f^infinity) and the
% % density perturbation should satisfy <density,1>_L2(f^infinity) = 0
% u0 = @u0func;
% % u0 = [1;0];
% setInitialConditions(model,u0);
% 
% %% solvepde steady state case
% % |solvepde| automatically picks the parabolic solver to obtain the
% % solution
% R = solvepde(model);
% u = R.NodalSolution;
% 
% figure; 
% pdeplot(model,'XYData',u(:,1),'Contour','on','ColorMap','jet');
% title 'Perturbation in density, Steady State Solution'
% xlabel 'X-coordinate, meters'
% ylabel 'Y-coordinate, meters'
% axis equal
% figure; 
% pdeplot(model,'XYData',u(:,2),'Contour','on','ColorMap','jet');
% title 'Perturbation in value function, Steady State Solution'
% xlabel 'X-coordinate, meters'
% ylabel 'Y-coordinate, meters'
% axis equal

%% Transient Solution %%

%% Definition of PDE Coefficients
% The expressions for the coefficients required by PDE Toolbox can easily
% be identified by comparing the equation above with the scalar parabolic
% equation in the PDE Toolbox documentation.

% The "m", "a" coefficients are zero, "d" is an identity matrix and the
% "c" coefficient is a constant matrix in the PDE system considered and 
% "f" coefficient is a vector depending on partial derivatives of u1 and 
% u2. Further f2 depends on a functional of u1, namely the CC perturbation.
% We define f by MATLAB function fcoeffunction whichin the MATLAB function 
% ccpertfunc will define the CC perturbation function.
m = zeros(4,1); % zeros(2,2); % 
d = eye(2); %zeros(4,1); % 
c = [0;0;0;sigma^2/2; 0;0;0;0; 0;0;0;1/2/r; 0;0;0;sigma^2/2];
a = [0;0;0;0]; % zeros(2,2); % 
f = @fcoeffunction;
specifyCoefficients(model,'m',0,'d',1,'c',c,'a',0,'f',f);
endTime = 4;
tlist = 0:.1:endTime;
p = msh.Nodes;
numNodes = size(p,2);

%% Initial and Boundary Conditions
% Boundary conditions are assumed to be both the density and value function
% to be zero on the boundary since both these functions are assumed to tend
% to zero as the radius tends to infinity. Thus we only need to set the
% Dirichilet boundary condition in this case and not Neumann conditions

% applyBoundaryCondition(model,'dirichlet','Edge',1,'u',[0; 0]);
% applyBoundaryCondition(model,'dirichlet','Edge',2,'u',[0; 0]);
% applyBoundaryCondition(model,'dirichlet','Edge',3,'u',[0; 0]);
% applyBoundaryCondition(model,'dirichlet','Edge',4,'u',[0; 0]);
applyBoundaryCondition(model,'dirichlet','edge',1:4,'h',[1,0;0,1],'r',[0;0]);
% bc2 = applyBoundaryCondition(model,'neumann','Edge',1:4,'g',[0;0], 'q',zeros(2,2));

% Initial conditions for PDE system. Initial perturbations should satisfy
% certain properties, namely belonging to the class L2(f^infinity) and the
% density perturbation should satisfy <density,1>_L2(f^infinity) = 0
u0 = @u0func;
% u0 = [1;0];
setInitialConditions(model,u0);
% setInitialConditions(model,[0.01;0.01]);

%% Set solver options
model.SolverOptions.RelativeTolerance = 1.0e-3; 
model.SolverOptions.AbsoluteTolerance = 1.0e-4;

%% solvepde transient case
% |solvepde| automatically picks the parabolic solver to obtain the
% solution
R = solvepde(model,tlist);
u = R.NodalSolution;

%% plotting results
figure;
pdeplot(model,'XYData',u(:,1,end),'Contour','on','ColorMap','jet');
title(sprintf('Perturbation in density, Transient Solution( %d seconds)\n', ...
  tlist(1,end)));
xlabel 'X-coordinate, meters'
ylabel 'Y-coordinate, meters'
axis equal;
figure;
pdeplot(model,'XYData',u(:,2,end),'Contour','on','ColorMap','jet');
title(sprintf('Perturbation in value function, Transient Solution ( %d seconds)\n', ...
  tlist(1,end)));
xlabel 'X-coordinate, meters'
ylabel 'Y-coordinate, meters'
axis equal;

%% 
toc