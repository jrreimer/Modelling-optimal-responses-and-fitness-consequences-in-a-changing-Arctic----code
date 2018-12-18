close all
clear 

% set working directory here, if needed
load('ResultWorkspace.mat')
x_init = 2000; % approximately the mean of the initial conditions vector
eta_init = 1; % begin single

cd 'Monte Carlo functions'

[MC_x, MC_eta, f_MC, dpatch_MC, dpreg_MC, dCOY_MC ] = monte_carlo_F_mat(...
    f, x_init, eta_init,N, t_breakup, tau_icefree, T, L, x_crit, x_max,...
    sigma, sigma_hat, tau_mate, tau_den, k, lambda, sigma_0, sigma_1);

%%%%% PLOT MONTE CARLO SIMULATION %%%%%
cd ..
cd 'Monte Carlo functions'
Monte_carlo_plot( T, x_crit, x_max, t_breakup, tau_icefree,...
     MC_x, MC_eta, L, tau_den)
set(gcf,'Position',[10, -80, 1000, 850])
% we further edited this image using the Matlab plot editing GUI

%%%%% Save plot, if desired 
% SET WORKING DIRECTORY HERE 
cd ..
% saveas(gcf,'Monte_carlo.png')
% saveas(gcf,'Monte_carlo.eps','epsc')
% saveas(gcf,'Monte_carlo.fig')





