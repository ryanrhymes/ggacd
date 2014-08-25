% Function to simulate a ACD(q,p) model to data
%
% USAGE: [simulDur]=ACD_Simul(nr,Coeff,q,p,dist)
%
% INPUT:
%       nr - Number of simulated observations
% 
%       Coeff - a structure with the simulated coefficients (more details
%       at example script)
% 
%       dist - distribution assumption (so far, only 'exp' and 'weibull'
%       and 'ggamma' accepted
% 
%       q - maximum lag at alpha coefficients
% 
%       p - maximum lag at beta coefficients 
% 
% OUTPUT:
%       specOut - A structure with the fitted coefficients and the fitted
%       duration (more details at example script).
% 
% Author: Marcelo Perlin (PhD student at ICMA/Uk)
% Modified by: Liang Wang (University of Helsinki, Finland

function [simulDur]=ACD_Simul(nr,Coeff,q,p,dist)

    import acd_garch.*;

    if nr<(max(p,q)+1)
        error('The number of observations shouldbe be higher than max(q,p)')
    end
    
    
    if p<1||q<1
        error('The input q and p should be higher or equal than 1');
    end
    
    if strcmp(dist,'exp')==0&&strcmp(dist,'weibull')==0&&strcmp(dist,'ggamma')==0
        error('The input dist should be either ''exp'' or ''weibull'' or ''ggamma''');
    end
    

    % Prealocation of large matrices and some precalculations
    
    simulDur=zeros(nr,1);
    durExp=zeros(nr,1);

    first_idx=max(q,p)+1;

    durExp(1,1)=0;
    simulDur(1:first_idx,1)=Coeff.w;
 
    % Matlab does not have ggamma random numbers
    if strcmp(dist,'ggamma')
        system(['python +tool/pyfun_ggamma_rvs.py ' num2str(Coeff.y) ' ' num2str(Coeff.z) ' ' int2str(nr) ' > ggrvs.tmp']);
        load ggrvs.tmp; delete ggrvs.tmp;
    end

    for i=first_idx:nr

        % this is the psi equation
        
        durExp(i,1)=Coeff.w+flipdim(Coeff.q,2)*simulDur(i-q:i-1,1)+flipdim(Coeff.p,2)*durExp(i-p:i-1,1);
        
        switch dist
            
            case 'exp'
                simulDur(i,1)=durExp(i,1)*exprnd(1);            % simulated duration at observation i
            case 'weibull'
                simulDur(i,1)=durExp(i,1)*wblrnd(1,Coeff.y);    % simulated duration at observation i
            case 'ggamma'
                simulDur(i,1)=durExp(i,1)*ggrvs(i);             % This is python fix for GGAMMA.

        end

    end

    plot(simulDur);
    title('Simulation of a Duration Series');
    xlabel('Observations');
    ylabel('Simulated Durations');
