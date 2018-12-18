function [ N, t_breakup, tau_icefree, T, L, x_crit, x_max, sigma,...
    sigma_hat, tau_mate, tau_den, k, lambda, sigma_0, sigma_1]= params_F( ~ )

% model parameters
N = 500;     % number of compartments between x_min and x_max 
t_breakup = 108; % number of days before breakup (= July 17 - April 1 + 1)
tau_icefree = 83; % number of ice free days during which bears must fast
T = 24;      % max number of years from reproductive maturity until death at age 30

% bear parameters
L = 1.96;    % asymptotic body length of adult female bear
x_crit = 0;  % below this, bear dies (KJ)
x_max = 26.14*59.76*L^3-390.53*L^3; % bear is satiated, (MJ) (tab: Energy 2 Mass conversions)
sigma = 0.996^(1/365); % daily probability of survival
sigma_hat = sigma^(365-t_breakup); % overwinter probability of survival 
tau_mate = 17; % length of pairing during mating
tau_den = 134; % number of days in maternity den
k = 1.15;      % expected size of a recruited litter

% patch specific parameters
lambda = [1/3.5 1/2.5];   % food encounter probability in patch [L, P]

val0 = 0.651^(1/365); 
scaleby = 3.4539; % multiple of mortality in active ice 
sigma_0 = [val0 1-(1-val0)*scaleby]; % daily probability of COY litter survival [L, P]
val1 = 0.86^(1/365);
scaleby = 1.1;
sigma_1 = [val1 1-(1-val1)*(scaleby)]; % daily probability of yearling litter survival [L, P]

end




