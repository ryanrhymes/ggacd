% Function for calculating standard errors of MS_Regress_Fit
%
% The calculation (aproximation) of the first and second derivative of the likelihood
% function is done by one side finite differences method
%
% Four methods for the calculation of the covariance matrix were implemented here:
% (1) Using second partial derivatives
% (2) Using first partial derivatives (outer product Matrix)
% (3) Using white covariance matrix 
% (4) Using Newey and West Covariance matrix
% 
% The standard erros from the last two methods are also called robust
% standard errors. 
%
% References:
%
% HAMILTON, J., D. (1994). Time Series Analysis. Princeton University Press.
%
% NEWEY, B., WEST, D (1987) �A Simple, Positive semidefinite, Heteroskedasticity
% and Autocorrelation Consistent Covariance Matrix� Econometrica, ed. 55, p. 347-370
%
% WHITE, H. (1984) �Asymptotic Theory for Econometricians� New York:
% Academic Press.
% Modified by: Liang Wang (University of Helsinki, Finland

function [V]=getvarMatrix_ACD(method,param,x,q,p,dist,indep)

import acd_garch.*;

fprintf(1,'\nCalculating Standard Error Vector...');

% A simple (and perhaps brutal) fix for the cases where the parameters are very small. If that
% happens, the derivatives of the lik function goes to zero, making the
% standard errors = NAN. The next line should take care of that by
% replacing any value lower than the machine precision

param(abs(param)<eps)=.0001;

% Definition of a small change for each parameter.

myDelta=1e-5*abs(param);

% First Derivative calculation


[sumlik1,Output,logLikVec1]=ACD_Lik(param,x,q,p,dist,0,indep);
n=length(logLikVec1); % number of observations
s=zeros(n,numel(param));
for i=1:numel(param)

    m=param;
    m(i)=param(i)+myDelta(i);
    [sumlik2,Output,logLikVec2]=ACD_Lik(m,x,q,p,dist,0,indep);

    s(:,i)=(logLikVec2-logLikVec1)/(myDelta(i));

end

sum_s_Matrix=zeros(numel(param));
for idx=2:n
    s_Matrix=s(idx,:)'*s(idx,:);
    sum_s_Matrix=sum_s_Matrix+s_Matrix;
end

OP_Matrix=sum_s_Matrix/n;   % outer product matrix

% Matrix for newey_west

switch method

    case 4
        
        S_param=zeros(numel(param));

        s_sum_part=zeros(length(param));
        for j=2:n

            sum_der=zeros(length(param));
            for i=j+1:n

                der=s(i,:)'*s(i-j,:);
                sum_der=sum_der+der;

            end

            myGamma=sum_der/n;

            s_sum_part=s_sum_part+(1-j/(n+1))*(myGamma+myGamma');

        end

        S_param=OP_Matrix+s_sum_part;   % Matrix for newey_west

end

% Hessian (second partial derivatives Matrix) Calculation

sde=zeros(numel(param));
for i=1:numel(param)

    for j=1:numel(param)

        m1=param;
        m2=param;

        m1(i)=param(i)-myDelta(i);   %x-
        m11=m1;
        m12=m1;
        m11(j)=m1(j)-myDelta(j);     %x-,y-
        m12(j)=m1(j)+myDelta(j);     %x-,y+
        
        m2(i)=param(i)+myDelta(i);   %x+
        m22=m2;
        m21=m2;
        
        m21(j)=m2(j)-myDelta(j);     %x+,y-
        m22(j)=m2(j)+myDelta(j);     %x+,y+  
        
        [sumlik11]=ACD_Lik(m11,x,q,p,dist,0,indep);
        [sumlik12]=ACD_Lik(m12,x,q,p,dist,0,indep);
        [sumlik21]=ACD_Lik(m21,x,q,p,dist,0,indep);
        [sumlik22]=ACD_Lik(m22,x,q,p,dist,0,indep);
        
        % the likeilhood funct provides the negative of log likelihood, so
        % I need to take the negative again
        
        sumlik11=-sumlik11;
        sumlik12=-sumlik12;        
        sumlik21=-sumlik21;        
        sumlik22=-sumlik22;                
        
        % second derivative by one side finite difference
        
        sde(i,j)=(sumlik22-sumlik21-sumlik12+sumlik11)/(4*myDelta(i)*myDelta(j));

    end
end

H=-sde/n; % the Hessian (or the information Matrix) See Hamilton eq 5.8.2, page 143

switch method
    case 1  % Using second partial derivatives
        V=inv(H*n);
    case 2  % Using first partial derivatives (outer product Matrix)
        V=inv(OP_Matrix*n);
    case 3  % Using White's Covariance Matrix (see eq 5.8.7 at Hamilton page 145)
        V=1/n*inv((H*inv(OP_Matrix)*H));
    case 4  % Using Newey and West Covariance Matrix (see eq 5.8.7 at Hamilton page 145)
        V=1/n*inv((H*inv(S_param)*H));
end