
% Calculates and plots changes in the reproductive energy thresholds for 
% varying tbreakup dates, as well as concurrent changes in the female’s 
% expected lifetime fitness. Figure 6 in text.

%note: should make N high in params_F.m for good resolution
close all
clear 
%tic

cd ..
[ N, ~, ~, T, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~] = params_F();
    
breakup_range = 87:1:108; % 
xvec = [87 94 101 108]; 
xlabels = {'June 26','July 3', 'July 10', 'July 17'};
preg.vals = zeros(length(breakup_range),T-1);
COY.vals = zeros(length(breakup_range),T-1);
R0 = zeros(N+1,length(breakup_range));

%%%%% PART A %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c = 1;
for t_breakup = breakup_range
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
    cd functions
    f = terminal_cond_F( t_breakup, f, N, T, k);

    % solve and output the dynamic programming equation
    [ dpatch, dpreg, dCOY, f, vL, vP ] = dpe_solve_F_mat( dpatch, dpreg,...
        dCOY, f, vL, vP, N, t_breakup, tau_icefree, T, L, x_crit, x_max, sigma,...
        sigma_hat, tau_mate, tau_den, k, lambda, sigma_0, sigma_1);
    
    X = dpreg == 2;
    Y = cumsum(X,1)==1 & X;
    [jval,year] = find(Y == 1);
    xvals = x_crit + ((x_max - x_crit)/N).*jval;
    preg.vals(c,year) = xvals;
    
    X = dCOY == 2;
    Y = cumsum(X,1)==1 & X;
    [jval,year] = find(Y == 1);
    xvals = x_crit + ((x_max - x_crit)/N).*jval;
    COY.vals(c,year) = xvals;
    
    R0(:,c) = f(:,1,1,1);
    
    c = c+1;
end

cd ..
%csvwrite('pregvals.csv',preg.vals)
%csvwrite('COYvals.csv',COY.vals)
%csvwrite('R0vals.csv',R0)

%%%%%% PLOT REPRODUCTIVE THRESHOLDS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

preg.vals = csvread('pregvals.csv');
COY.vals = csvread('COYvals.csv');
R0 = csvread('R0vals.csv');

figure
subplot(2,1,1)
age = 10;
plot(breakup_range, preg.vals(:,age-4),'k:','LineWidth',3)
hold on
plot(breakup_range, COY.vals(:,age-4),'k','LineWidth',3)
xlim([breakup_range(1),breakup_range(length(breakup_range))])
ylim([2600,3700])
%xlabel('t_{breakup}','FontSize',12)
ylabel('reproductive threshold (MJ)','FontSize',12)
%title(['Reproductive thresholds for an individual age ',num2str(age)])
xticks(xvec)
xticklabels(xlabels)
legend('abort pregnancy','cease lactation')
legend('boxoff')
text(87.5,3625,'(a)')
set(gca,'fontsize', 12);
%saveas(gcf,'Reproductive_thresholds_against_icefree.png')
%saveas(gcf,'Reproductive_thresholds_against_icefree.eps','epsc')
%saveas(gcf,'Reproductive_thresholds_against_icefree.fig')


%%%%%% PLOT CHANGES IN R_0 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%figure
subplot(2,1,2)
cd .. 
Init.x = csvread('Initial_x.csv'); % remember to make sure whole csv is NUMERIC

x = mean(Init.x);
%put real values into computer values ranging from 0 to N
xc = (x-x_crit)*N/(x_max - x_crit);
j = floor(xc);
change_x = xc - j;

d = 1;
R0x = zeros(1,length(breakup_range)); 
for t_breakup = breakup_range
    R0x(1,d) = (1-change_x)*R0(j+1,d) + change_x*R0(j+2,d); % linear interpolated f values
    d = d+1;
end

plot(breakup_range, R0x,'k','LineWidth',3)
xlabel('t_{breakup}','FontSize',12)
ylabel({'expected lifetime','reproductive success'},'FontSize',12)
xlim([breakup_range(1),breakup_range(length(breakup_range))])
xticks(xvec)
xticklabels(xlabels)
text(87.5,3.8,'(b)')

%title('R_0 as a function of the length of spring','FontSize',12)
set(gca,'fontsize',12)
set(gcf,'color','white','Position',[100, 100, 375, 600])

%saveas(gcf,'R0_and_thresholds_against_icefree.png')
%saveas(gcf,'R0_and_thresholds_against_icefree.eps','epsc')
%saveas(gcf,'R0_and_thresholds_against_icefree.fig')


%toc
