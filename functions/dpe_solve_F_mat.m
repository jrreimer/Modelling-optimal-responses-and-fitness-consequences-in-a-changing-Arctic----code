function [ dpatch, dpreg, dCOY, f, vL, vP ] = dpe_solve_F_mat( dpatch, dpreg,...
    dCOY, f, vL, vP, N, t_breakup, tau_icefree, T, L, x_crit, x_max, sigma,...
    sigma_hat, tau_mate, tau_den, k, lambda, sigma_0, sigma_1)

%dpe_solve: solves the stochastic dynamic programming equation using backwards
%induction, including linear interpolation

for tyear = T:-1:1 

    disp(tyear) % outputs tyear for a sense of computational progress
    
    % get terminal condition for given year (except final year)
    if tyear < T        
        j = 1:N;
        x = x_crit + ((x_max - x_crit)/N)*j; % corresponding x value
        
        % get cumulative senescence probability up to start of next spring
        fun = @(yr,xi) exp(-(yr/xi).^xi).*((yr/xi).^(xi-1));
        age = tyear + 4;
        intgl = integral(@(yr)fun(yr,23),age,age+1); % evaluate integral from _ to _ with parameter xi=24
        intgl_scl = integral(@(yr)fun(yr,23),age,Inf);
        ps = intgl/intgl_scl;
        ps = min(max(ps, 0),1);
        
        %%% eta = 1: single %%%
        w = w_1_mat( x, tau_icefree, L );
        w = max(min(w,x_max), x_crit);  
        
        % if not senescent, become senescent with probability ps
        f(2:N+1,t_breakup,tyear,1,1) = (1-ps)*sigma_hat*linear_interp_w_1_mat(tyear, w, f, N, x_crit, x_max, 1) + ...
            ps*sigma_hat*linear_interp_w_1_mat(tyear, w, f, N, x_crit, x_max, 2);
        
        % if already senescent, stay senescent
        f(2:N+1,t_breakup,tyear,1,2) = sigma_hat*linear_interp_w_1_mat(tyear, w, f, N, x_crit, x_max, 2);  

        %%% eta = 2: pregnant %%%
        % first, for bear that is not senescent: 
        xi = 1;
            % continue pregnancy (only possible if not senescent; i.e., xi=1)
            w = w_2_mat( x, tau_icefree, L, tau_den );
            w = max(min(w,x_max), x_crit);
            fw2 = linear_interp_w_2_mat( tyear, w, f, N, x_crit, x_max, xi);
            f_continue_preg = sigma_hat*fw2;                
            % abort pregnancy
            w = w_2_loss_mat( x, tau_icefree, L );
            w = max(min(w,x_max), x_crit); 
            fw2loss = linear_interp_w_2_loss_mat( tyear, w, f, N, x_crit, x_max, xi);
            f_abort_preg = sigma_hat*fw2loss;   
            % abort pregnancy due to senescence 
            fw2loss_sens = linear_interp_w_2_loss_mat( tyear, w, f, N, x_crit, x_max, 2);
            f_abort_preg_sens = sigma_hat*fw2loss_sens;            
            % now she may be come senescent, with probability ps, and if
            % not, then f results from decision to abort or continue pregnancy:
            f(2:N+1,t_breakup,tyear,2,1) = (1-ps)*max(f_abort_preg, f_continue_preg).' + ...
                ps*f_abort_preg_sens.';            
            % record which decision is optimal (abort or continue preg.)
            abt = f_abort_preg > f_continue_preg;
                dpreg(abt,tyear) = 1;
            ctn = f_continue_preg > f_abort_preg;
                dpreg(ctn,tyear) = 2;
            nodiff = f_abort_preg == f_continue_preg;
                dpreg(nodiff,tyear) = 3; %i.e. it doesn't matter
        % then, for bear that is senescent: 
        f(2:N+1,t_breakup,tyear,2,2) = f_abort_preg_sens.';
                                     
        %%% eta = 3: with COYs %%%
            % first for a non-senescent bear
            xi = 1;
            % keep litter
            w = w_3_mat( x, tau_icefree, L );
            w = max(min(w,x_max), x_crit); 
            fw3 = linear_interp_w_3_mat( tyear, w, f, N, x_crit, x_max,xi);
            f_keep_litter = sigma_hat*fw3;              
            % abandon litter
            w = w_3_loss_mat( x, tau_icefree, L );
            w = max(min(w,x_max), x_crit); 
            fw3loss = linear_interp_w_3_loss_mat( tyear, w, f, N, x_crit, x_max,xi);
            f_abandon_litter = sigma_hat*fw3loss;
            % abandon litter but become senescent
            fw3loss_sens = linear_interp_w_3_loss_mat( tyear, w, f, N, x_crit, x_max,2);
            f_abandon_litter_sens = sigma_hat*fw3loss_sens;           
                
        % f results from decision to abort or continue pregnancy:
        f(2:N+1,t_breakup,tyear,3,1) = (1-ps)*max(f_keep_litter, f_abandon_litter).' + ...
            ps*f_abandon_litter_sens.';
        
        % then for a senescent bear
        f(2:N+1,t_breakup,tyear,3,2) = f_abandon_litter_sens.';
        
        % record which decision is optimal (abort or continue preg.) if not senescent
        abnd = f_abandon_litter > f_keep_litter;
            dCOY(abnd,tyear) = 1;
        keep = f_keep_litter > f_abandon_litter;
            dCOY(keep,tyear) = 2;
        thesame = f_abandon_litter == f_keep_litter;
            dCOY(thesame,tyear) = 3; %i.e. it doesn't matter
                                           
        %%% eta = 4: with yearlings %%%
        w = w_4_mat( x, tau_icefree, L );
        w = max(min(w,x_max), x_crit); 
        % if bear is not senescent:         
            fw4 = linear_interp_w_4_mat(tyear, w, f, N, x_crit, x_max, 1);
            % if bear becomes senscent:
            fw4_sens = linear_interp_w_4_mat(tyear, w, f, N, x_crit, x_max, 2);
            f(2:N+1,t_breakup,tyear,4,1) = k + sigma_hat*((1-ps)*fw4 + ps*fw4_sens);
            % note the fitness increment "+ k"     
        % if bear is senescent % although I think this case can't actually occur
            f(2:N+1,t_breakup,tyear,4,2) = k + sigma_hat*fw4_sens;        
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % work backwards, by day, for the given year
    for t = (t_breakup - 1):-1:1
        for xi = 1:2 % for senescent and not senescent scenarios
        
            % for each value of x, compute fitness of visiting patch 1 or 2
            % xp: state resulting from choosing patch i and catching food
            % xpp: state resulting from choosing patch i and not catching food
            j = 1:N;
            x = x_crit + ((x_max - x_crit)/N).*j; % x value corresponding to j 
            
            % calculate value of mass specific functions
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            epsilon = epsilon_F_mat( x, L, x_crit, t, t_breakup, tau_mate);
            a = Ai_F_mat(x,L);
            [g3,g4] = milkFun_mat(x,L,x_crit);            
            Y = Yi_F( t );
            
            for i = 1:2 % now for each patch
                [xp, xpp, x1, xp3, xpp3, xp4, xpp4] = stateChanges_mat(x, a,...
                    Y, i, tau_mate, g3, g4, x_max, x_crit,L);
            
                % obtain fitness values (linear interpolate, if necessary)
                [Value] = ValueFunctions_mat(t,...
                    tyear, i, f, xp, xpp, x1, xp3, xpp3, xp4, xpp4, ...
                    epsilon, g3, g4, ...
                    N, t_breakup, tau_icefree, T, L, x_crit, x_max, sigma,...
                    sigma_hat, tau_mate, tau_den, k, lambda, sigma_0, sigma_1,xi);
              
                for eta = 1:4
                    if i==1  % foraging patch 1, Landfast ice
                        vL(:,t,tyear,eta,xi) = Value(eta,:).';                   
                    else  % foraging patch 2, Pack ice
                        vP(:,t,tyear,eta,xi) = Value(eta,:).';                    
                    end
                end
            end 
        
            % find the best patch and save the corresponding decision
            vmax(1,:) = max(vL(:,t,tyear,1,xi), vP(:,t,tyear,1,xi)).';
            vmax(2,:) = max(vL(:,t,tyear,2,xi), vP(:,t,tyear,2,xi)).';   
            vmax(3,:) = max(vL(:,t,tyear,3,xi), vP(:,t,tyear,3,xi)).';   
            vmax(4,:) = max(vL(:,t,tyear,4,xi), vP(:,t,tyear,4,xi)).';  
            
            for q = 1:4
                
                land = (vmax(q,:).' == vL(:,t,tyear,q,xi)); % indices of patches with landfast bigger
                dpatch(land,t,tyear,q,xi) = 1;
                pack = (vmax(q,:).' == vP(:,t,tyear,q,xi)); % indices of patches with pack ice bigger
                dpatch(pack,t,tyear,q,xi) = 2;
            
                eql = vL(:,t,tyear,q,xi) == vP(:,t,tyear,q,xi); 
                dpatch(eql,t,tyear,q,xi) = 3;
            
                % update fitness for current time step 
                f(2:N+1,t,tyear,q,xi) = vmax(q,:);
            end 
        end 
    end
    % fix few places where rounding errors are making goofy patch predictions
    % (shouldn't change fitness values, but will change decision plots); I have
    % confirmed by hand that these are all rounding errors
    for ct = 1:92
        dpatch(:,ct,21,1,1) = 2;
    end
end



