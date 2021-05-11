function [setting]= init_cma_par1(setting, EA_2, n, n2)
%% So, mean(EA_2) does not violate the initialization condition in the competition
setting.xmean = mean(EA_2);
setting.xmean=setting.xmean';
setting.insigma=0.01;
setting.sigma = setting.insigma;


setting.sigma = max(setting.insigma);              % overall standard deviation
setting.pc = zeros(n,1); setting.ps = zeros(n,1);  % evolution paths for setting.C and setting.sigma

if length(setting.insigma) == 1
    setting.insigma = setting.insigma * ones(n,1) ;
end
setting.diagD = setting.insigma/max(setting.insigma);      % diagonal matrix D defines the scaling
setting.diagC = setting.diagD.^2;
setting.B = eye(n,n);                      % setting.B defines the coordinate system
setting.BD = setting.B.*repmat(setting.diagD',n,1);        % setting.B*D for speed up only
setting.C = diag(setting.diagC);                   % covariance matrix == setting.BD*(setting.BD)'
setting.D = ones(n,1); 
setting.chiN=n^0.5*(1-1/(4*n)+1/(21*n^2));  % expectation of
setting.mu = ceil(n2/2);               % number of parents/points for recombination
% setting.mu = n2;
setting.weights = log(max(setting.mu, n/2) + 1/2)-log(1:setting.mu)'; % muXone array for weighted recombination setting.mu = floor(setting.mu);
setting.mueff=sum(setting.weights)^2/sum(setting.weights.^2); % variance-effective size of setting.mu
setting.weights = setting.weights/sum(setting.weights);     % normalize recombination setting.weights array

% Strategy parameter setting: Adaptation
setting.cc = (4 + setting.mueff/n) / (n+4 + 2*setting.mueff/n); % time constant for cumulation for setting.C
setting.cs = (setting.mueff+2) / (n+setting.mueff+3);  % t-const for cumulation for setting.sigma control
setting.ccov1 = 2 / ((n+1.3)^2+setting.mueff);    % learning rate for rank-one update of setting.C
setting.ccovmu = 2 * (setting.mueff-2+1/setting.mueff) / ((n+2)^2+setting.mueff);  % and for rank-setting.mu update
setting.damps = 0.5 + 0.5*min(1, (0.27*n2/setting.mueff-1)^2) + 2*max(0,sqrt((setting.mueff-1)/(n+1))-1) + setting.cs; % damping for setting.sigma

setting.xold =  setting.xmean;
