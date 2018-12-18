
% Runs 100 Monte Carlo simulations for varying lengths of spring. Used to
% create Figure 4 in text. 
% Note: the simulations take a long time to run. If just wish to create
% plot, just run Part B below. 

close all
clear

%%% PART A %%%
cd ..
Init_x = csvread('Initial_x.csv'); % remember to make sure whole csv is NUMERIC
num_sims = 100; % number of simulated model runs
eta_init = 1;
breakup_range = 87:7:108; % 107 and 1,2,3 weeks earlier

perpack_COY = zeros(num_sims, length(breakup_range)); % each row is a different simulation
perpack_yrl = zeros(num_sims, length(breakup_range)); % each row is a different simulation

for c = 1:4
    t_breakup = breakup_range(c);
    % run SDP model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp('t_breakup = ')
    disp(t_breakup)    
    tau_icefree = 281 - (t_breakup+91); % tau_breakup increases with changes in t_breakup

    f = [];   % F(x,t,tyear,eta):  current fitness function = max(v1,v2)
    vL = [];  % V1(x,t,tyear,eta): value function for patch 1
    vP = [];  % V2(x,t,tyear,eta): value function for patch 2
    dpatch = [];   % dpatch(x,t,tyear,eta):  optimal patch decision
    dpreg = [];    % dpreg(x,tyear): abort or keep pregnancy overwinter decision
    dCOY = [];     % dCOY(x,tyear): abandon or keep COY litter overwinter decision

    % assign parameters 
    [ N, ~, ~, T, L, x_crit, x_max, sigma, sigma_hat, ...
        tau_mate, tau_den, k, lambda, sigma_0, sigma_1] = params_F();

    % set end condition
    cd 'functions'
    f = terminal_cond_F( t_breakup, f, N, T, k);
    % make f, vL, and VP into symbolic matrices

    % solve and output the dynamic programming equation
    [ dpatch, dpreg, dCOY, f, vL, vP ] = dpe_solve_F_mat( dpatch, dpreg,...
        dCOY, f, vL, vP, N, t_breakup, tau_icefree, T, L, x_crit, x_max, sigma,...
        sigma_hat, tau_mate, tau_den, k, lambda, sigma_0, sigma_1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    initvals = datasample(Init_x, num_sims); % randomly sample from distribution of initial values (with replacement)

    parfor j = 1:num_sims
        x_init = initvals(j);
        cd 'C:\Users\Jody\Google Drive\PhD\Optimal Predation Project\Matlab code\Adult females with variable senescence\Used Files\Monte Carlo functions'
        [MC_x, MC_eta, f_MC, dpatch_MC, dpreg_MC, dCOY_MC ] =...
            monte_carlo_F_mat(  f, x_init, eta_init, ...
            N, t_breakup, tau_icefree, T, L, x_crit, x_max, sigma, sigma_hat, ...
        tau_mate, tau_den, k, lambda, sigma_0, sigma_1);
        
        %withcubs = (MC_eta(1,1:t_breakup-1,1:T) == 3) | (MC_eta(1,1:t_breakup-1,1:T) == 4); % indices when eta = 3 or 4
        withcubs_COY = (MC_eta(1,1:t_breakup-1,1:T) == 3); % indices when eta = 3 or 4
        withcubs_yrl = (MC_eta(1,1:t_breakup-1,1:T) == 4);
        patches_COY = dpatch_MC(withcubs_COY);
        patches_yrl = dpatch_MC(withcubs_yrl);
        perpack_COY(j,c) = sum(patches_COY==2)/length(patches_COY); 
        perpack_yrl(j,c) = sum(patches_yrl==2)/length(patches_yrl); 
    end  
end % end t_breakup in breakup_range
cd ..
cd 'Plotting Code'
csvwrite('perpack_COY.csv',perpack_COY)
csvwrite('perpack_yrl.csv',perpack_yrl)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% PART B - PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
perpack_COY = csvread('perpack_COY.csv');
perpack_yrl = csvread('perpack_yrl.csv');
perpack_COY(perpack_COY == 0) = NaN;
perpack_yrl(perpack_yrl == 0) = NaN;
figure
subplot(1,2,1)
boxplot(perpack_COY, 'Colors','k','Symbol','ok','OutlierSize',4)
title('with cubs')
xlabels = {'JN 26','JL 3', 'JL 10', 'JL 17'};
xticklabels(xlabels)
xlabel('t_{breakup}','FontSize',12)
ylabel('% days in active ice ','FontSize',12)
ylim([0 1])
set(gca,'fontsize',12)
subplot(1,2,2)
boxplot(perpack_yrl, 'Colors','k','Symbol','ok','OutlierSize',4)
title('with yearlings')
xlabels = {'June 26','July 3', 'July 10', 'July 17'};
xticklabels(xlabels)
xlabel('t_{breakup}','FontSize',12)
ylim([0 1])
set(gca,'fontsize',12)
set(gcf,'color','white','Position',[100, 100, 725, 300])

%saveas(gcf,'Percent_time_in_patch.png')
%saveas(gcf,'Percent_time_in_patch.eps','epsc')
%saveas(gcf,'Percent_time_in_patch.fig')

for ct = 1:4
    col1 = perpack_COY(:,ct);
    col1(isnan(col1)) = [];
    median(col1)
end
for ct = 1:4
    col1 = perpack_yrl(:,ct);
    col1(isnan(col1)) = [];
    median(col1)
end

