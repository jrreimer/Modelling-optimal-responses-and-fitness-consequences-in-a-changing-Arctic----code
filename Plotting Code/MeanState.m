%%% 10000 Monte Carlo simulations of a 10 year old bear’s energy stores 
% throughout spring. Creates Figure 3 in text. 

close all
clear

%%% Note: this script takes a long time to run; if just want to create the
%%% figure, can skip to Part B below %%%%%%

%%%%%%%%%% PART A %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 10000 Mote Carlo simulations (i.e., the grey lines in Figure 3)

cd ..
load('ResultWorkspace.mat') % first load workspace results
Initx = csvread('Initial_x.csv'); % note: whole csv must be numeric

numsims = 10000; % number of simulated model runs
age = 10;
year = age-4;

eta_init = 1;
x_all = zeros(numsims, t_breakup, T); % each row is a different simulation
eta_all = zeros(numsims, t_breakup, T);
init_vals = datasample(Initx, numsims); % randomly sample from distribution of initial values (with replacement)

cd 'Monte Carlo functions'
parfor j = 1:numsims
    disp('Simulation Number:')
    disp(j)
    
    x_init = init_vals(j);
    cd 'C:\Users\Jody\Google Drive\PhD\Optimal Predation Project\Matlab code\Adult females with variable senescence\Used Files\Monte Carlo functions'
    [MC_x, MC_eta, f_MC, dpatch_MC, dpreg_MC, dCOY_MC ] = ...
        monte_carlo_F_mat(  f, x_init, eta_init,...
        N, t_breakup, tau_icefree, T, L, x_crit, x_max, sigma, sigma_hat, ...
            tau_mate, tau_den, k, lambda, sigma_0, sigma_1);
      
    x_all(j,:,:) = MC_x;   
    eta_all(j,:,:) = MC_eta;
end

mean_x = zeros(1,t_breakup,T);
mean_x_1 = zeros(1,t_breakup,T);
mean_x_2 = zeros(1,t_breakup,T);
mean_x_3 = zeros(1,t_breakup,T);
mean_x_4 = zeros(1,t_breakup,T);

x_all_1 = zeros(size(x_all)); x_all_1(eta_all == 1) = x_all(eta_all == 1);
x_all_2 = zeros(size(x_all)); x_all_2(eta_all == 2) = x_all(eta_all == 2);
x_all_3 = zeros(size(x_all)); x_all_3(eta_all == 3) = x_all(eta_all == 3);
x_all_4 = zeros(size(x_all)); x_all_4(eta_all == 4) = x_all(eta_all == 4);

for k = 1:T
    for l = 1:t_breakup
        live.vals = x_all((x_all(:,l,k) > x_crit),l,k);
        mean_x(1,l,k) = mean(live.vals);
        
        live.vals1 = x_all_1((x_all_1(:,l,k) > x_crit),l,k);
        live.vals2 = x_all_2((x_all_2(:,l,k) > x_crit),l,k);
        live.vals3 = x_all_3((x_all_3(:,l,k) > x_crit),l,k);
        live.vals4 = x_all_4((x_all_4(:,l,k) > x_crit),l,k);
        mean_x_1(1,l,k) = mean(live.vals1);
        mean_x_2(1,l,k) = mean(live.vals2);
        mean_x_3(1,l,k) = mean(live.vals3);
        mean_x_4(1,l,k) = mean(live.vals4);          
    end
end
%%%% STOPPED RUNNING IT HERE 
cd ..
cd 'Plotting Code'
% save outputs (they take a long time; this way you can call them later)
csvwrite('mean_x_1.csv',mean_x_1(:,:,year))
csvwrite('mean_x_2.csv',mean_x_2(:,:,year))
csvwrite('mean_x_3.csv',mean_x_3(:,:,year))
csvwrite('mean_x_4.csv',mean_x_4(:,:,year))
csvwrite('mean_x.csv',mean_x(:,:,year))
csvwrite('x_all.csv',x_all(:,:,year))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PART B %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% now for a given year, mid-age (say age 10), plot lines in grey, and
% average in black against date
x_all = csvread('x_all.csv');
mean_x_1 = csvread('mean_x_1.csv');
mean_x_2 = csvread('mean_x_2.csv');
mean_x_3 = csvread('mean_x_3.csv');
mean_x_4 = csvread('mean_x_4.csv');
mean_x = csvread('mean_x.csv');

% for dates on x axis over spring
springDayz = datetime(2000,4,01) + caldays(0:(t_breakup-1));
figure
for j = 1:numsims
    %plot(1:t_breakup,x_all(j,:,year),'LineWidth',0.5,'Color',rgb('lightgrey'))
    plot(springDayz, x_all(j,:,year),'LineWidth',0.5,'Color',rgb('lightgrey')); hold on
end
% calculate average mass of alive simulated bears

hold on
plot(springDayz, mean_x(1,:,year),'k','LineWidth',3)

%%%%% can plot each reproductive class as its own colour %%%%%
%plot(springDayz, mean_x_1(1,:),'LineWidth',3,'Color',rgb('yellow')); hold on
%plot(springDayz, mean_x_2(1,:),'LineWidth',3,'Color',rgb('green')); hold on
%plot(springDayz, mean_x_3(1,:),'LineWidth',3,'Color',rgb('blue')); hold on
%plot(springDayz, mean_x_4(1,:),'LineWidth',3,'Color',rgb('pink')); hold on

xtickformat('MMM d')
xlabel('day in spring, t','FontSize',12)
ylabel('Energy reserves (MJ)','FontSize',12)
%title(['Energetic state through spring for a bear age ',num2str(age)])
set(gca,'fontsize',12)
xlim([springDayz(1) springDayz(length(springDayz))])
%ylim([x_crit x_max])
ylim([x_crit 7000])
set(gcf,'color','white','Position',[100, 100, 375, 300])

% cd 'C:\Users\Jody\Google Drive\PhD\Optimal Predation Project\Result plots' % save all plots here
% saveas(gcf,'MeanState.png')
% saveas(gcf,'MeanState.eps','epsc')
% saveas(gcf,'MeanState.fig')