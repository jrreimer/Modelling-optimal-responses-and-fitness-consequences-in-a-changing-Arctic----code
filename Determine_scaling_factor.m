
%%% runs many (1000) Monte Carlo simulations and saves the proportion of 
% spring days that the simulated individual spends in the active ice when 
% she has a litter of cubs of the year. The user must then use matlab’s 
% curve fitting GUI to fit an exponential curve. The parameters of this 
% curve are then entered into the bottom part of the script to plot the 
% scaling factor. Creates Figure 1 in the text. 

close all
clear

% NOTE: This script takes a long time to run,
% however. If you wish just to create the plot, using the existing data
% points, skip down to Section B below. 

%%%%% Section A %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Monte carlo simulations (grey dots in figure 1). 

% set working directory, if needed
%cd 'C:\Users\Jody\Google Drive\PhD\Optimal Predation Project\Matlab code\Adult females with variable senescence\Used Files'
num.sims = 100; % number of simulated model runs
eta_init = 1; % start in reproductive state 1

% reasonable bounds on scaling values
min_scale = 2;  
max_scale = 5;

Init.x = csvread('Initial_x.csv'); % initial conditions of body condition data
x_init = datasample(Init.x, num.sims); % randomly sample from distribution of initial values (with replacement)

perpack = zeros(num.sims, 1); % each row is a different simulation; each column is percent of time spent in pack ice
scl = unifrnd(min_scale, max_scale, num.sims, 1); % vector of random scalings to try

parfor c = 1:num.sims
    % run SDP model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp(c)   % simulation number counter

    f = [];   % F(x,t,tyear,eta):  current fitness function = max(v1,v2)
    vL = [];  % V1(x,t,tyear,eta): value function for patch 1
    vP = [];  % V2(x,t,tyear,eta): value function for patch 2
    dpatch = [];   % dpatch(x,t,tyear,eta):  optimal patch decision
    dpreg = [];    % dpreg(x,tyear): abort or keep pregnancy overwinter decision
    dCOY = [];     % dCOY(x,tyear): abandon or keep COY litter overwinter decision

    % assign parameters 
    [ N, t_breakup, tau_icefree, T, L, x_crit, x_max, sigma, sigma_hat, ...
        tau_mate, tau_den, k, lambda, sigma_0, sigma_1] = params_F();
    sigma_0(2) = max(1-(1-sigma_0(1))*scl(c),0); % survival probability for coys in pack ice
    sigma_1(2) = max(1-(1-sigma_1(1))*1.1, 0); % survival probability for yearlings in pack ice
    
    % set terminal condition
    cd functions
    f = terminal_cond_F( t_breakup, f, N, T);

    % solve and output the dynamic programming equation
   [ dpatch, dpreg, dCOY, f, vL, vP ] = dpe_solve_F_mat( dpatch, dpreg,...
        dCOY, f, vL, vP, N, t_breakup, tau_icefree, T, L, x_crit, x_max, sigma,...
        sigma_hat, tau_mate, tau_den, k, lambda, sigma_0, sigma_1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cd ..
    cd 'Monte Carlo functions'
    [MC_x, MC_eta, f_MC, dpatch_MC, dpreg_MC, dCOY_MC ] =...
            monte_carlo_F_mat(  f, x_init(c), eta_init, ...
            N, t_breakup, tau_icefree, T, L, x_crit, x_max, sigma, sigma_hat, ...
            tau_mate, tau_den, k, lambda, sigma_0, sigma_1);
      
    withcubs = (MC_eta(1,1:t_breakup-1,1:T) == 3); 
    patches = dpatch_MC(withcubs);
    disp('percent time in pack ice:')
    disp(sum(patches==2)/length(patches))    
    perpack(c,1) = sum(patches==2)/length(patches);  

end % end simulations
% save csv outputs for later (they take long to run)
csvwrite('scl.csv',scl)
csvwrite('perpack.csv',perpack)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% plot  %%%%%%%%%%%%%%%%%%
% if ran previously, or on a server, can read in csv outputs here:
cd ..
scl = csvread('scl.csv');
perpack = csvread('perpack.csv');

%call curve fitting tool:
cftool
% get coefficients from curve fitting tool
% use x data as scl (scale) and y data as perpack

tar = 0.37; % target percent time in pack ice
%%% exponential %%%
a = 0.8884 ;
b = -0.2536;
scalefac = (1/b)*log(tar/a);
x = min_scale:0.01:max_scale;
perpackfit = a*exp(b.*x);

disp('scaling factor that results in 37% of time spent in the pack ice is:') 
disp(scalefac)

%%%% make figure
figure
scatter(scl, perpack(:,1),10,rgb('lightgrey'),'filled')
hold on
plot(x, perpackfit, 'k','linewidth',2)
line([min_scale scalefac], [tar tar],'Color','black','LineStyle','--');
hold on
line([scalefac scalefac], [0 tar],'Color','black','LineStyle','--'); 
hold on
plot(scalefac, tar,'s', ...
    'MarkerSize',12,...
    'MarkerEdgeColor','black',...
    'MarkerFaceColor','black')
xlabel('scaling factor','FontSize',12)
ylabel('% of days in active ice ','FontSize',12)
legend('simulations','fit exponential curve','Location','NE');
legend('boxoff')
set(gca,'fontsize',12)
xlim([min_scale max_scale])
ylim([0 0.83])
set(gcf,'color','white','Position',[100, 100, 375, 300])

%%% save plots %%%
% saveas(gcf,'Determine_scaling_factor.png')
% saveas(gcf,'Determine_scaling_factor','epsc')
% saveas(gcf,'Determine_scaling_factor.fig')



