% Equality and inequality constraints 
% Modified by: Liang Wang (University of Helsinki, Finland

function [c,ceq]=confuneq(param)

    global global_p;
    global global_q;

    p=global_p;
    q=global_q;
    
    Coeff.w=param(1);
    Coeff.q=param(2:q+1);
    Coeff.p=param(q+2:q+1+p);

    c=sum(Coeff.q)+sum(Coeff.p)-1; % the sum of alpha and beta should be lower than 1
    ceq=0;