function [] = Monte_carlo_plot( T, x_crit, x_max, t_breakup, tau_icefree,...
    MC_x, MC_eta, L, tau_den)

% go to appropriate working directory
cd ..
cd functions

map = [ 0.961 0.933 0.451
    0.314 0.663 0.506
    0.494 0.329 0.647
    0.961 0.588 0.451];

if T ~= 24
    disp('Not correct number of years, fix subplots')
end
    
figure
hold on
set(gcf,'color','white','Position',[10, -80, 1000, 850])

yvec_xvals = [2500 5000 7500];
ylabels = {'x_{crit}', '2500', '5000', '7500', 'x_{max}'};

%dim = [0.82 0.12 .3 .3];
dim = [0.50 0.12 .3 .3];
str = {'Reproductive',' senescence'};
%str = {'Reproductive senescence'};
annotation('textbox',dim,'String',str,'FitBoxToText','on',...
    'EdgeColor','none','VerticalAlignment','bottom',...
    'HorizontalAlignment','left','FontSize',12);

x1 = .13;
x2 = .127; 

y1 = .887;
y2 = (.887/4)*3;
y3 = (.887/4)*2;
y4 = (.887/4)*1;

%%% Add age labels over each year %%%
for age = 1:5
    dim = [(x1 + (age-1)*0.1290) y1 .3 .3];
    str = ['Age ' num2str(age+4)];
    annotation('textbox',dim,'String',str,'FitBoxToText','on',...
        'EdgeColor','none','VerticalAlignment','bottom',...
        'HorizontalAlignment','left','FontSize',12);    
end
age = 6;
dim = [0.7720 y1 .3 .3];
str = ['Age ' num2str(age+4)];
    annotation('textbox',dim,'String',str,'FitBoxToText','on',...
        'EdgeColor','none','VerticalAlignment','bottom',...
        'HorizontalAlignment','left','FontSize',12);
for age = 1:6
    dim = [(x2 + (age-1)*0.1290) y2 .3 .3];
    str = ['Age ' num2str(age+10)];
    annotation('textbox',dim,'String',str,'FitBoxToText','on',...
        'EdgeColor','none','VerticalAlignment','bottom',...
        'HorizontalAlignment','left','FontSize',12);    
end
for age = 1:6
    dim = [(x2 + (age-1)*0.1290) y3 .3 .3];
    str = ['Age ' num2str(age+16)];
    annotation('textbox',dim,'String',str,'FitBoxToText','on',...
        'EdgeColor','none','VerticalAlignment','bottom',...
        'HorizontalAlignment','left','FontSize',12);    
end
for age = 1:6
    dim = [(x2 + (age-1)*0.1290) y4 .3 .3];
    str = ['Age ' num2str(age+22)];
    annotation('textbox',dim,'String',str,'FitBoxToText','on',...
        'EdgeColor','none','VerticalAlignment','bottom',...
        'HorizontalAlignment','left','FontSize',12);    
end


for c = 1:4
    %%%%% set up each subplot %%%%%
    subplot(4,1,c)
    plot(1:6*365,ones(1,6*365),'w')
    axis([1,6*365+1,x_crit,x_max+1])
    yticks([x_crit yvec_xvals x_max])
    yticklabels(ylabels)
    if c == 2
        ht = text(-175,-9000,'Energy reserves (MJ)');
        set(ht,'Rotation',90)
        set(ht,'FontSize',12)   
    end
    if c == 4
        xlabel('Date','FontSize',12)
    end
    t_springs = 1:365:(5*365+1);
    t_breakups = t_springs +  t_breakup - 1;
    t_freezeups = t_breakups + tau_icefree;
    all_xticks = [t_springs t_breakups t_freezeups];
    xtx = sort(all_xticks);
    xlb = {'t_{spring}','t_{breakup}','t_{freezeup}'};
    xlbx = repmat(xlb,1,6);        
    if c ~= 4        
        xticks(xtx);        
        xticklabels(xlbx)
    elseif c == 4
        xticks(xtx(1:length(xtx)-1))
        xticklabels(xlbx(1:length(xlbx)-1)) 
    end
    xtickangle(-20)
    box off
    
    %%%%% For each year: %%%%%    
    for tyear = (1:6) + (c-1)*6 
        yr = mod(tyear,6);
        if yr == 0
            yr = 6;
        end
        Time = ((yr-1)*365 +1):((yr-1)*365 + t_breakup); % spring days
         
         for j = 1:t_breakup-1
            if MC_eta(1,j,tyear) == 1
                line(Time(j:j+1),MC_x(:,j:j+1,tyear),'Color',map(1,:),'LineWidth',3)
            elseif MC_eta(1,j,tyear) == 2
                line(Time(j:j+1),MC_x(:,j:j+1,tyear),'Color',map(2,:),'LineWidth',3)
            elseif MC_eta(1,j,tyear) == 3
                line(Time(j:j+1),MC_x(:,j:j+1,tyear),'Color',map(3,:),'LineWidth',3)
            elseif MC_eta(1,j,tyear) == 4
                line(Time(j:j+1),MC_x(:,j:j+1,tyear),'Color',map(4,:),'LineWidth',3)
            end            
         end
         
         mksz = 9;
         for j = 1:t_breakup-1
            if MC_eta(1,j,tyear) == 1
                if MC_eta(1,j+1,tyear) == 2
                    hold on
                    plot(Time(j+1),MC_x(1,j+1,tyear),'ko','MarkerSize',mksz,'MarkerFaceColor','k') % became pregnant
                end
            elseif MC_eta(1,j,tyear) == 2
            elseif MC_eta(1,j,tyear) == 3
                if MC_eta(1,j+1,tyear) == 1 
                    hold on
                    plot(Time(j+1),MC_x(1,j+1,tyear),'kx','MarkerSize',mksz,'LineWidth',2,'MarkerFaceColor','k') % lost COYs
                end
            elseif MC_eta(1,j,tyear) == 4
                if MC_eta(1,j+1,tyear) == 1 
                    hold on
                    plot(Time(j+1),MC_x(1,j+1,tyear),'kx','MarkerSize',mksz,'LineWidth',2,'MarkerFaceColor','k') % lost COYs
                end
            end            
         end
         
         %%% SEE WHAT HAPPENED OVER WINTER
         if tyear ~= T % don't need to check in final year
            eta1 = MC_eta(1,t_breakup,tyear);
            eta2 = MC_eta(1,1,tyear+1);
            %%%%% finished spring single %%%%%
            if eta1 == 1 
                xsum = zeros(1,tau_icefree+1); % store summer states
                xwin = ones(1, 365 - t_breakup - tau_icefree + 2);
                xsum(1) = MC_x(1,t_breakup,tyear);
                for d = 2:tau_icefree+1
                    xsum(d) = xsum(d-1) - 0.392*(mass_mat(xsum(d-1),L).^0.813);
                end
                xwin = xwin*w_1_mat( MC_x(1,t_breakup,tyear), tau_icefree, L );
                line(Time(j+1):Time(j+1)+tau_icefree,xsum,'Color',map(1,:),'LineWidth',3,'LineStyle',':');
                line(Time(j+1)+tau_icefree: Time(j+1)+tau_icefree + 365 - t_breakup - tau_icefree+1, xwin, 'Color', map(1,:), 'LineWidth',3,'LineStyle',':');

                
            %%%%% finished spring pregnant %%%%%
            elseif eta1 == 2 % finished spring pregnant
                xsum = zeros(1,tau_icefree+1); % store summer states
                xsum(1) = MC_x(1,t_breakup,tyear);
                for d = 2:tau_icefree+1
                    xsum(d) = xsum(d-1) - 0.392*(mass_mat(xsum(d-1),L).^0.813);
                end
                
                if eta2 == 1 % abort litter
                    xwin = ones(1, 365 - t_breakup - tau_icefree + 2);
                    xwin = xwin*w_1_mat( MC_x(1,t_breakup,tyear), tau_icefree, L );
                    line(Time(j+1):Time(j+1)+tau_icefree,xsum,'Color',map(2,:),'LineWidth',3,'LineStyle',':');
                    line(Time(j+1)+tau_icefree: Time(j+1)+tau_icefree + 365 - t_breakup - tau_icefree+1, xwin, 'Color', map(1,:), 'LineWidth',3,'LineStyle',':');
                    
                elseif eta2 == 3 % continued litter
                    xwin = ones(1, 365 - t_breakup - tau_icefree - tau_den + 2);
                    xden = ones(1, tau_den+1);
                    xden(1) = xsum(tau_icefree+1);
                    for f = 2:tau_den + 1
                        xden(f) = xden(f-1) - 0.02*mass_mat(xden(f-1),L).^1.09;
                    end
                    xwin = xwin*w_2_mat( MC_x(1,t_breakup,tyear), tau_icefree, L, tau_den);
                    line(Time(j+1):Time(j+1)+tau_icefree,xsum,'Color',map(2,:),'LineWidth',3,'LineStyle',':');
                    line(Time(j+1)+tau_icefree:Time(j+1)+tau_icefree+tau_den, xden, 'Color', map(2,:), 'LineWidth',3,'LineStyle',':');
                    line(Time(j+1)+tau_icefree+tau_den:Time(j+1)+tau_icefree+tau_den + 365 - t_breakup - tau_icefree - tau_den +1, xwin, 'Color', map(3,:), 'LineWidth',3,'LineStyle',':');
                end
                
            %%%%% finished spring with COYs %%%%%
            elseif eta1 == 3                 
                xwin = ones(1, 365 - t_breakup - tau_icefree + 2);
                
                if eta2 == 1 % abandoned COYs
                    % tau_icefree - round(tau_icefree/2) + 1
                    xsumA = zeros(1,round(tau_icefree/2)+1); % store summer states with COYS
                    xsumB = zeros(1,tau_icefree - round(tau_icefree/2) + 1); % store summer states after losing COYs
                    xsumA(1) = MC_x(1,t_breakup,tyear);

                    for h = 2:(round(tau_icefree/2)+1)
                        xsumA(h) = xsumA(h-1) - 0.392*mass_mat(xsumA(h-1),L).^0.813 - 0.17*mass_mat(xsumA(h-1),L).^0.75;
                    end 
                    xsumB(1) = xsumA(h);
                    for k = 2:(tau_icefree - round(tau_icefree/2) + 1)
                        xsumB(k) = xsumB(k-1) - 0.392*mass_mat(xsumB(k-1),L).^0.813;
                    end
                    xwin = xwin*w_3_loss_mat( MC_x(1,t_breakup,tyear), tau_icefree, L);
                    
                    line(Time(j+1):Time(j+1)+round(tau_icefree/2),xsumA,'Color',map(3,:),'LineWidth',3,'LineStyle',':');
                    line(Time(j+1)+round(tau_icefree/2):Time(j+1)+tau_icefree,xsumB,'Color',map(1,:),'LineWidth',3,'LineStyle',':');
                    line(Time(j+1)+tau_icefree: Time(j+1)+tau_icefree + 365 - t_breakup - tau_icefree+1, xwin, 'Color', map(1,:), 'LineWidth',3,'LineStyle',':');
                    
                elseif eta2 == 4 % stayed with COYs
                    xsum = zeros(1, tau_icefree+1);
                    xsum(1) = MC_x(1,t_breakup,tyear);
                    for l = 2:tau_icefree+1
                        xsum(l) = xsum(l-1) - 0.392*mass_mat(xsum(l-1),L).^0.813 - 0.17*mass_mat(xsum(l-1),L).^0.75;
                    end
                    xwin = xwin*w_3_mat( MC_x(1,t_breakup,tyear), tau_icefree, L);
                    line(Time(j+1):Time(j+1) + tau_icefree,xsum,'Color',map(3,:),'LineWidth',3,'LineStyle',':');
                    line(Time(j+1)+tau_icefree: Time(j+1)+tau_icefree + 365 - t_breakup - tau_icefree+1, xwin, 'Color', map(3,:), 'LineWidth',3,'LineStyle',':');
                 end
                
                
            %%%%% finished spring with yearlings %%%%%
            else 
                xsum = zeros(1,tau_icefree+1); % store summer states
                xwin = ones(1, 365 - t_breakup - tau_icefree + 2);
                xsum(1) = MC_x(1,t_breakup,tyear);
                for d = 2:tau_icefree+1
                    xsum(d) = xsum(d-1) - 0.392*(mass_mat(xsum(d-1),L).^0.813);
                end
                xwin = xwin*w_1_mat( MC_x(1,t_breakup,tyear), tau_icefree, L );
                line(Time(j+1):Time(j+1)+tau_icefree,xsum,'Color',map(4,:),'LineWidth',3,'LineStyle',':');
                line(Time(j+1)+tau_icefree: Time(j+1)+tau_icefree + 365 - t_breakup - tau_icefree+1, xwin, 'Color', map(4,:), 'LineWidth',3,'LineStyle',':');
             end
         end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for tyear = (1:6) + (c-1)*6
        yr = mod(tyear,6);
        if yr == 0
            yr = 6;
        end
        Time = ((yr-1)*365 +1):((yr-1)*365 + t_breakup); % spring days
        %xlabel(['Age ',num2str(tyear+4)], 'FontSize',12)
         
         %%%%%%%%%%%% MAKE MARKERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         if tyear ~= T % don't need to check in final year
            eta1 = MC_eta(1,t_breakup,tyear);
            eta2 = MC_eta(1,1,tyear+1);
            %%%%% finished spring pregnant %%%%%
            if eta1 == 2 % finished spring pregnant
                xsum = zeros(1,tau_icefree+1); % store summer states
                xsum(1) = MC_x(1,t_breakup,tyear);
                for d = 2:tau_icefree+1
                    xsum(d) = xsum(d-1) - 0.392*(mass_mat(xsum(d-1),L).^0.813);
                end
                
                if eta2 == 1 % abort litter
                    hold on
                    plot(Time(j+1)+tau_icefree,xsum(length(xsum)),'kx','MarkerSize',mksz,'LineWidth',2,'MarkerFaceColor','k') % abort litter

                elseif eta2 == 3 % continued litter
                    xden = ones(1, tau_den+1);
                    xden(1) = xsum(tau_icefree+1);
                    for f = 2:tau_den + 1
                        xden(f) = xden(f-1) - 0.02*mass_mat(xden(f-1),L).^1.09;
                    end
                    hold on
                    plot(Time(j+1)+tau_icefree+tau_den,xden(length(xden)),'ks','MarkerSize',mksz,'MarkerFaceColor','k') % abort litter
                end
                
            %%%%% finished spring with COYs %%%%%
            elseif eta1 == 3                 
                xwin = ones(1, 365 - t_breakup - tau_icefree + 2);
                
                if eta2 == 1 % abandoned COYs
                    % tau_icefree - round(tau_icefree/2) + 1
                    xsumA = zeros(1,round(tau_icefree/2)+1); % store summer states with COYS
                    xsumA(1) = MC_x(1,t_breakup,tyear);

                    for h = 2:(round(tau_icefree/2)+1)
                        xsumA(h) = xsumA(h-1) - 0.392*mass_mat(xsumA(h-1),L).^0.813 - 0.17*mass_mat(xsumA(h-1),L).^0.75;
                    end 
                    hold on
                    plot(Time(j+1)+round(tau_icefree/2),xsumA(length(xsumA)),'kx','MarkerSize',mksz,'LineWidth',2,'MarkerFaceColor','k') % abort litter
         
                elseif eta2 == 4 % stayed with COYs
                    xsum = zeros(1, tau_icefree+1);
                    xsum(1) = MC_x(1,t_breakup,tyear);
                    for l = 2:tau_icefree+1
                        xsum(l) = xsum(l-1) - 0.392*mass_mat(xsum(l-1),L).^0.813 - 0.17*mass_mat(xsum(l-1),L).^0.75;
                    end
                    xwin = xwin*w_3_mat( MC_x(1,t_breakup,tyear), tau_icefree, L);
                    hold on
                    plot(Time(j+1)+tau_icefree + 365 - t_breakup - tau_icefree+1,xwin(length(xwin)),'kd','MarkerSize',mksz,'MarkerFaceColor','k') % abort litter
                end           
                
            %%%%% finished spring with yearlings %%%%%
            elseif eta1 == 4
                xwin = ones(1, 365 - t_breakup - tau_icefree + 2);
                xwin = xwin*w_1_mat( MC_x(1,t_breakup,tyear), tau_icefree, L );
                hold on
                plot(Time(j+1)+tau_icefree + 365 - t_breakup - tau_icefree+1,xwin(length(xwin)),'k*','LineWidth',2,'MarkerSize',mksz,'MarkerFaceColor','k') % abort litter
            end
         end
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         if tyear == T % don't need to check in final year
            eta1 = MC_eta(1,t_breakup,tyear);            
            %%%%% finished spring pregnant or with COYS %%%%%
            if (eta1 == 2) || (eta1 == 3) 
                hold on
                plot(Time(j+1),MC_x(1,j+1,tyear),'kx','MarkerSize',mksz,'LineWidth',2,'MarkerFaceColor','k') % lose litter due to senescence
            %%%%% finished spring with yearlings %%%%%
            elseif eta1 == 4
                hold on
                plot(Time(j+1),MC_x(1,j+1,tyear),'k*','LineWidth',2,'MarkerSize',mksz,'MarkerFaceColor','k') % recruit litter
            end
         end
         
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% Plot Legend %%%%%
    for ii = 1:size(map,1)
        p(ii) = line(NaN, NaN, 'Color',map(ii,:),'LineWidth',3);
    end
    lbl =  {'Single', 'Pregnant','With cubs','With yearlings'};
    legend(p, lbl,'Orientation','horizontal',...
        'Position',[0.485 0.95 0.4 0.02],'FontSize',12);
    
    set(gca,'fontsize',12) 
end
end
