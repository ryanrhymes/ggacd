% Example script ACD_Simul and ACD_Fit for ACD(1,1)
% It will first simulate a ACD model and then fit it using the fitting
% function
% Modified by: Liang Wang (University of Helsinki, Finland

clc; clear all; close all;

import acd_garch.*;

%%%%%%%%%%%%%%%%%%%%%%%% OPTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nr=1000;    % How many observations in the simulation

% Choose your distribution (just comment the one you dont want)
% dist='exp';
dist='ggamma';

% Choose your flavor (parameters)

Coeff.w=.1;     % constant in expected duration (psi)
Coeff.q=.2;     % Coeff at duration in t-1 (alpha)
Coeff.p=.7;     % Coeff at expected duration in t-1 (beta)
Coeff.y=.8;     % for weibull and ggamma dist
Coeff.z=1.2;    % just for ggamma dist

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

q=size(Coeff.q,2);
p=size(Coeff.p,2);

simulDur=ACD_Simul(nr,Coeff,q,p,dist);  % Simulation