% Example script ACD_Simul and ACD_Fit for ACD(1,1)
% It will first simulate a ACD model and then fit it using the fitting
% function
% Modified by: Liang Wang (University of Helsinki, Finland

clc; clear all; close all;

import acd_garch.*;

%%%%%%%%%%%%%%%%%%%%%%%% OPTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nr=10000;    % How many observations in the simulation

% Choose your distribution (just comment the one you dont want)

% dist='exp';
dist='ggamma';

% Choose your flavor (parameters)

Coeff.w=.1;     % constant in expected duration (psi)
Coeff.q=.03;     % Coeff at duration in t-1 (alpha)
Coeff.p=.95;     % Coeff at expected duration in t-1 (beta)
Coeff.y=.4;     % for weibull and ggamma dist
Coeff.z=3.2;    % just for ggamma dist

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

q=size(Coeff.q,2);
p=size(Coeff.p,2);

simulDur=ACD_Simul(nr,Coeff,q,p,dist);  % Simulation

[specOut]=ACD_Fit(simulDur,dist,q,p);   % Fitting

plot([specOut.h simulDur]);
title('Duration Simulation and Modelling');
legend('Fitted Duration', 'Real Duration');
xlabel('Observations');
ylabel('Durations');