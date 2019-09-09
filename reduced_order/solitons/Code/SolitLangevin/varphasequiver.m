% CONFIDENTIAL (C) Mitsubishi Electric Research Labs (MERL) 2018
% Author: Kaivalya Bakshi; Date: 25 Apr 2018; Script for plotting soliton  
% dynamics from Quadratic MFG paper by Ullmo et. al 
% (https://arxiv.org/pdf/1708.07730.pdf) with external potential 
% U0(x) = -x^4/4 - x^2/2. We extend the published example by admitting 
% coupled and non trivial passive drift in the SDE dynamics

clc
clear all
close all

par.alpha = 3;
par.g = 1;
mu = 1;
sigma = 1;

g = par.g;
alpha = par.alpha;

k = g/(1 + alpha)/sqrt(1 + alpha)/(2*pi)^(alpha/2);

qstaranalyticc0 = (mu*sigma^4/(4*k*alpha))^(1/(2-alpha)) % for c = 0, alpha ~= 2

[q2,p2] = meshgrid(.05:0.1:0.5,-1:0.08:1);

q2dot = p2./mu;
% for uncoupled soliton variance dynamics
p2dot = -mu.*sigma.^4./(4.*q2.^3) + k.*alpha./q2.^(alpha+1); % integrator system

figure(1)
x = min(min(q2)):0.01:max(max(q2));
plot(x,mu*sigma^2./2./x - x*1,'--r')
grid on
hold on
scale = 1;
LineSpec = 'filled';
quiver(q2,p2,q2dot,p2dot,scale)
plot(qstaranalyticc0,0,'r*')