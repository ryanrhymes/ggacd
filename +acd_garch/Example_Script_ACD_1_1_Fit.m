% Example script ACD_Simul and ACD_Fit for ACD(1,1)
% It will first load some duration data and then fit a ACD(1,1) by using the fitting
% function
% Modified by: Liang Wang (University of Helsinki, Finland

clc; clear all; close all;

import acd_garch.*;

load +acd_garch/Example_Data.mat;

%%%%%%%%%%%%%%%%%%%%%%%% OPTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Choose your distribution (just comment the one you dont want)

%%% dist='weibull'; or dist='exp';
dist='ggamma';

stdMethod=1;   % method for standard error calculation (see ACD_Fit.m)

% Choose your flavor (parameters)

q=1;
p=1;
x=dur;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[specOut]=ACD_Fit(x,dist,q,p,stdMethod);    % Fitting

plot([specOut.h x]);
title('Duration Simulation and Modelling');
legend('Fitted Duration', 'Real Duration');
xlabel('Observations');
ylabel('Durations');
