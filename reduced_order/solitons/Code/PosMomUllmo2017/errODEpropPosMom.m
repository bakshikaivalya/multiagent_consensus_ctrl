% CONFIDENTIAL (C) Mitsubishi Electric Research Labs (MERL) 2018
% Author: Kaivalya Bakshi; Date: 27 Mar 2018; Function for computing
% terminal time error related to terminal boundary condition of position,
%  momentum ODEs of example in figure 1 of Quadratic MFG paper by Ullmo 
% et. al (https://arxiv.org/pdf/1708.07730.pdf), by forward propagating 
% position, momentum ODEs given an initial time guess for the 2 point
% boundary value problem (BVP). We extend the published example by 
% admitting non trivial dynamics in the SDE dynamics

function  error = errODEpropPosMom(t0,tf,X0,P0,PosMomDyn,errtermPosMom)
XP0=[X0;P0];
[T,X]=feval('ode45',PosMomDyn,[t0 tf],XP0);
XPf=X(length(X),:)'; % extracing ODE solution at terminal time from ODE solution propagated using ode45
error=feval(errtermPosMom,XPf);
end