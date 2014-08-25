% Example script ACD_Simul and 
% This script will fisrt simulate a ACD (1,1) with exponential distribution
% and then plot the likelihood of the model (the figure showed at matlab
% exchange site)
% Modified by: Liang Wang (University of Helsinki, Finland

% clear;

addpath('mFiles_ACD');

%%%%%%%%%%%%%%%%%%%%%%%% OPTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nr=1000;            % How many observations in the simulation of the duration series
n_simulations=5000; % How many simulations at likelihood plot

% Choose duration series parameters (only alpha and beta);

Coeff.q=.2; % Coeff at duration in t-1 (alpha)
Coeff.p=.7; % Coeff at expected duration in t-1(beta)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Coeff.w=.1; % This isn't an option (dont change it)

simulDur=ACD_Simul(nr,Coeff,1,1,'exp');  % Simulation

plotLikelihoodFunct(simulDur,n_simulations); % This is the function for plotting the log likelihoods

rmpath('mFiles_ACD');
