% Modified by: Liang Wang (University of Helsinki, Finland

function plotLikelihoodFunct(x,n_simulations)

if n_simulations<0
    error('The input n_simulations should be higher than zero');
end

w=.1;   % this is the constant of the model

sumloglik=zeros(n_simulations,1);

a_points=rand(n_simulations,1); % simulating alpha
b_points=rand(n_simulations,1); % simulating beta

for j1=1:n_simulations

    h1=w+filter(a_points(j1),1,x);          % Vectorized Filter
    h=filter([0 1],[1 -b_points(j1)],h1);   % Vectorized Filter
    
    h(1,1)=w;   % setting first 

    loglik(:,1)=log(1./h(:,1).*exp(-x(:,1)./h(:,1)));   % Exponential log likelihood vector

    sumloglik(j1)=sum(loglik);  %sum of log likelihoods

    disp(['Simulation #' num2str(j1) ' of ' num2str(n_simulations)]);

end

% plot it!

figure('position',[50 80 1150 650]);

plot3(b_points,a_points,sumloglik,'.','MarkerSize',15);
title('Log Likelihoods for ACD Model');

xlabel('Betas');
ylabel('Alphas');
zlabel('LogLikelihoods');
grid on;

[idx1 ]=find(max(max((sumloglik)))==sumloglik);

hold on;

plot3(b_points(idx1),a_points(idx1),sumloglik(idx1),'.','color','r','MarkerSize',35);
legend('Log Likelihoods for different \alpha and \beta',['Maximum Log Lik at \alpha=' num2str(a_points(idx1)), 'and \beta=' num2str(b_points(idx1))]);

hold off;