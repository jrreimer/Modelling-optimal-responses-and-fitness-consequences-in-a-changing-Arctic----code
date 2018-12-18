% This is the dynamic programming model for bear foraging habitat. 
% Plots the best decisions at each time for a given physiological state. 

% Patch 1 is the landfast ice, and Patch 2 the active ice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set working directory here if needed 

close all
clear

% initialize arrays
f = [];   % F(x,t,tyear,eta):  current fitness function = max(vL,vP)
vL = [];  % V1(x,t,tyear,eta): value function for patch 1, landfast ice
vP = [];  % V2(x,t,tyear,eta): value function for patch 2, pack ice
dpatch = [];   % dpatch(x,t,tyear,eta):  optimal patch decision
dpreg = [];    % dpreg(x,tyear): abort or keep pregnancy overwinter decision
dCOY = [];     % dCOY(x,tyear): abandon or keep COY litter overwinter decision

% assign parameters 
[ N, t_breakup, tau_icefree, T, L, x_crit, x_max, sigma, sigma_hat, ...
    tau_mate, tau_den, k, lambda, sigma_0, sigma_1] = params_F();

cd functions % set working directory to sub-folder "functions"

f = terminal_cond_F( t_breakup, f, N, T); % set end condition

% solve and output the dynamic programming equation
[ dpatch, dpreg, dCOY, f, vL, vP ] = dpe_solve_F_mat( dpatch, dpreg,...
    dCOY, f, vL, vP, N, t_breakup, tau_icefree, T, L, x_crit, x_max, sigma,...
    sigma_hat, tau_mate, tau_den, k, lambda, sigma_0, sigma_1);


cd ..   % change working directory back to parent directory
save('ResultWorkspace.mat') % save workspace; includes results and parameters

