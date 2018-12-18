% to plot results generated from pb_foraging_dyn_prog_model_F_mat.m
close all
clear 
cd ..
load 'ResultWorkspace.mat'


%cd 'C:\Users\Jody\Google Drive\PhD\Optimal Predation Project\Result plots' % save all plots here
darkval = 0.35;
medval = 0.6;
lightval = 0.85;
mapIce = [darkval darkval darkval
    medval medval medval
    lightval lightval lightval];

mapCOY = [darkval darkval darkval
    medval medval medval
    lightval lightval lightval];

fontsz = 12;
yvec_xvals = [2500 5000 7500];
yvec_indexvals = (yvec_xvals-x_crit).*N./(x_max - x_crit);
ylabels = {'x_{max}','7500','5000','2500','x_{crit}'};
xvec = [18 18+30 18+30+31];
xlabels = {'AP', 'MA', 'JN'};
%%%% OPTIMAL DECISIONS - EACH REPRODUCTIVE STATE %%%%%%%%%%%%%%%%%%%%%%%%%
for eta = 1:4
    figure
    
    for j = 1:T % plot all T years
        subplot(4,6,j)
        image(flip(dpatch(:,:,j,eta),1))
        colormap(mapIce)
        set(gca,'xtick',[])
        set(gca,'ytick',[])
        if j >= 19 
            xlabel('date','FontSize',fontsz)
            xticks(xvec)
            xticklabels(xlabels)
        end
        if j == 1 || j == 7 || j == 13 || j == 19
            %[1 yvec_indexvals N]
            yticks([1 1.42 2.83 4.25 5])
            %yticks([1 yvec_indexvals N])
            yticklabels(ylabels)
        end
        if j == 7
            ylim=get(gca,'YLim');
            ht = text(-50,1035,'Energy reserves (MJ)');
            set(ht,'Rotation',90)
            set(ht,'FontSize',fontsz)   
        end       
        
        title(['Age ',num2str(j+4)],'FontSize',fontsz)
                
        if j == 1
            hold on
            for ii = 1:size(mapIce,1)
                p(ii) = patch(NaN, NaN, mapIce(ii,:));
            end
            lbl =  {'fast ice', 'active ice','either'};
            legend(p, lbl,'Orientation','horizontal','Position',[0.1 0.008 0.3 0.02]);
        end
        set(gca,'fontsize',fontsz) 

    end   
    set(gcf, 'Position', [10, 50, 1400, 700])
    
%     cd 'C:\Users\Jody\Google Drive\PhD\Optimal Predation Project\Result plots'   
%     if eta == 1
%         saveas(gcf,'Result_plots_eta_1.png')
%         saveas(gcf,'Result_plots_eta_1.eps','epsc')
%         saveas(gcf,'Result_plots_eta_1.fig')
%     elseif eta == 2
%         saveas(gcf,'Result_plots_eta_2.png')
%         saveas(gcf,'Result_plots_eta_2.eps','epsc')
%         saveas(gcf,'Result_plots_eta_2.fig')
%     elseif eta == 3
%         saveas(gcf,'Result_plots_eta_3.png')
%         saveas(gcf,'Result_plots_eta_3.eps','epsc')
%         saveas(gcf,'Result_plots_eta_3.fig')
%     else
%         saveas(gcf,'Result_plots_eta_4.png')
%         saveas(gcf,'Result_plots_eta_4.eps','epsc')
%     end
end


%%%%%%%%%% PLOT 4 OF THE ABOVE FIGURES in 1 PLOT %%%%%%%%%%%
figure
for eta = 1:4
    j = 6; % take 6th model year; bear = 10 years old
    subplot(2,2,eta)
    image(flip(dpatch(:,:,j,eta),1))

    colormap(mapIce)
    if eta == 1
        hold on
        for ii = 1:size(mapIce,1)
            p(ii) = patch(NaN, NaN, mapIce(ii,:));
        end
        lbl =  {'fast ice', 'active ice','either'};
        legend(p, lbl,'Location','northwest');
    end
    
    %colorbar
    xlabel('date','FontSize',12)
    if eta == 1
        title('single', 'FontSize',12)
        ylabel('Energy reserves (MJ)','FontSize',12)
    elseif eta == 2
        title('pregnant', 'FontSize',12)
    elseif eta == 3
        title('with cubs','FontSize',12)
        ylabel('Energy reserves (MJ)','FontSize',12)
    else
        title('with yearlings','FontSize',12)
    end
             
    yticks([1 yvec_indexvals N])
    yticklabels(ylabels)
    
    xticks(xvec)
    xticklabels(xlabels)
    set(gca,'fontsize',fontsz) 

end

set(gcf, 'Position', [10, 50, 650, 650])
% saveas(gcf,'Result_plots_all_eta.png')
% saveas(gcf,'Result_plots_all_eta.eps','epsc')
% saveas(gcf,'Result_plots_all_eta.fig')


%%%%%% PLOT OPTIMAL OVERWINTER DECISIONS %%%%%%%
figure
subplot(1,2,1)
image(flip(dpreg,1))
colormap(mapCOY)
hold on
for ii = 1:size(mapCOY,1)
    p(ii) = patch(NaN, NaN, mapCOY(ii,:));
end
lbl =  {'abort pregnancy', 'continue pregnancy','either'};
legend(p, lbl,'Location','northwest'); hold off

xlabel('age (years)','FontSize',fontsz)
ylabel('Energy reserves (MJ)','FontSize',fontsz)
title('pregnant')  
yticks([1 yvec_indexvals N])
xticks(1:3:26)
yticklabels(ylabels)
xticklabels(5:3:30)
set(gca,'fontsize',fontsz) 

subplot(1,2,2)
image(flip(dCOY,1))
colormap(mapCOY)
hold on
for ii = 1:size(mapCOY,1)
    p(ii) = patch(NaN, NaN, mapCOY(ii,:));
end
lbl =  {'abandon cubs', 'keep cubs','either'};
legend(p, lbl,'Location','northwest'); hold off

xlabel('age (years)','FontSize',fontsz)
xticks(1:3:26)
xticklabels(5:3:30)
title('with cubs')             

yticks([1 yvec_indexvals N])
yticklabels(ylabels)
set(gca,'fontsize',fontsz) 

set(gcf, 'Position', [50, 50, 650, 300])
% saveas(gcf,'Optimal_repro.png')
% saveas(gcf,'Optimal_repro.eps','epsc')
% saveas(gcf,'Optimal_repro.fig')
cd ..
