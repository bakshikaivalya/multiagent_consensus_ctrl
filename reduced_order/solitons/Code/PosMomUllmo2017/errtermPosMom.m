% Author: Kaivalya Bakshi; Date: 27 Mar 2018; Function for computing
% terminal time error related to terminal boundary condition of position,
%  momentum ODEs of example in figure 1 of Quadratic MFG paper by Ullmo 
% et. al (https://arxiv.org/pdf/1708.07730.pdf). We extend the published 
% example by admitting non trivial dynamics in the SDE dynamics

function e = errtermPosMom(X)
e = X(2,1) + X(1,1) - 7/2 ;
end