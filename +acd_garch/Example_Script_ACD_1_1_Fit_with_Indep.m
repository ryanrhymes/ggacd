% Example script ACD_Simul and ACD_Fit for ACD(1,1) with independent
% variables
% It will first load some duration data and then fit a ACD(1,1) by using the fitting
% function

% Please notes that the indep is in the same time as the duration (x). If
% you want to pass lagged variables you'll have to lag the the independent matrix 
% (eg. indep=indep(2:end,:) for first lag)
%
% You should also notes that when you use a indep matrix, the likelihood
% evaluation is much slower agasint the case of a plain ACD(p,q). The reason is
% that the later uses the vectorized function filter() for calculation of
% conditional duration, making it a lot faster.
% Modified by: Liang Wang (University of Helsinki, Finland

clear;
load Example_Data.mat;
addpath('mFiles_ACD');

%%%%%%%%%%%%%%%%%%%%%%%% OPTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Choose your distribution (just comment the one you dont want)

dist='exp';
% dist='weibull';

stdMethod=1;   % method for standard error calculation (see ACD_Fit.m)

% Choose your flavor (parameters)

q=1;
p=1;
x=dur;
indep=indep_example(:,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[specOut]=ACD_Fit(x,dist,q,p,stdMethod,indep);    % Fitting

plot([specOut.h x]);
title('Duration Simulation and Modelling');
legend('Fitted Duration', 'Real Duration');
xlabel('Observations');
ylabel('Durations');

rmpath('m_Files_ACD');