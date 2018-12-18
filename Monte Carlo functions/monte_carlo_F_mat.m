function [MC_x, MC_eta, f_MC, dpatch_MC, dpreg_MC, dCOY_MC ] = ...
    monte_carlo_F_mat(  f, x_init, eta_init,...
    N, t_breakup, tau_icefree, T, L, x_crit, x_max, sigma, sigma_hat, ...
        tau_mate, tau_den, k, lambda, sigma_0, sigma_1)
    % Monte Carlo simulations of individual(s) assuming optimal decisions 
    % made at each time step 
    cd 'C:\Users\Jody\Google Drive\PhD\Optimal Predation Project\Matlab code\Adult females with variable senescence\Used Files\functions'    
    % initialize arrays
    sensyear = 99; % preliminary value; if outputs 99 at end, bad sign
    MC_x = zeros(1, t_breakup, T);   % values of x in MC simulation (assuming optimal decisions made at each point)
    MC_eta = zeros(1, t_breakup, T); % values of eta in MC simulation    
    f_MC = zeros(1, t_breakup, T); % fitness at each time
    
    vL_MC = zeros(1, t_breakup-1, T); % value functions in each patch
    vP_MC = zeros(1, t_breakup-1, T); 
    
    dpatch_MC = zeros(1, t_breakup-1, T); % optimal decision
    dpreg_MC = zeros(1, T-1); % optimal pregnancy overwinter decision
    dCOY_MC = zeros(1, T-1); % optimal litter of COYs overwinter decision
        
    num_offspring = 0; % keeps track of total offspring
    MC_x(1,1,1) = x_init; % initial x condition
    MC_eta(1,1,1) = eta_init;
    x = x_init;
    eta = eta_init; 
    xi = 1; % senescence variable
    breakloop = 0; % helps get out of loops
          
for tyear = 1:T
    if breakloop == 1
        break
    else
        breakloop = 0;
    end
    disp('tyear =')
    disp(tyear)
    t = 1;
      
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%% Spring Starting condition depends on over winter %%%%%%%%%%%%
    if tyear > 1 % state is w(x) from over previous winter
        % get senescence probability for over spring from tyear-1 to tyear
        fun = @(yr,xi) exp(-(yr/xi).^xi).*((yr/xi).^(xi-1));
        age = tyear + 4;
        intgl = integral(@(yr)fun(yr,23),age-1,age); % evaluate integral from _ to _ with parameter xi=24
        intgl_scl = integral(@(yr)fun(yr,23),age-1,Inf);
        ps = intgl/intgl_scl;
        ps = min(max(ps, 0),1);
        
        if eta == 1 %single
            w = w_1_mat( x, tau_icefree, L );
            w = max(min(w,x_max), x_crit);              
            % if not senescent, become senescent with probability ps
            if xi == 1
                futuref = (1-ps)*sigma_hat*linear_interp_w_1_mat(tyear-1, w, f, N, x_crit, x_max, 1) + ...
                    ps*sigma_hat*linear_interp_w_1_mat(tyear-1, w, f, N, x_crit, x_max, 2);
            elseif xi == 2 % if already senescent, stay senescent
                futuref = sigma_hat*linear_interp_w_1_mat(tyear-1, w, f, N, x_crit, x_max, 2);  
            end
            f_MC(1,t_breakup,tyear-1) = futuref(1);
            if rand < ps % become senescent
                xi = 2;
                sensyear = tyear;
            end            
            if rand > sigma_hat % die over winter
                disp('Died over winter from random causes in year:')
                %disp(tyear)
                break
            else
                x = w;
                eta = 1;
                disp('single at start of spring')
            end            
        elseif eta == 2 % pregnant
            % if aborts pregnancy due to senescence 
            w_abt = w_2_loss_mat( x, tau_icefree, L );
            w_abt = max(min(w_abt,x_max), x_crit); 
            f_abort_preg_sens = linear_interp_w_2_loss_mat( tyear-1, w_abt, f, N, x_crit, x_max, 2);
            f_abort_preg_sens = f_abort_preg_sens(1);            
            if xi == 1 % 
                % first, determine optimal overwinter strategy
                % continue pregnancy
                w = w_2_mat( x, tau_icefree, L, tau_den );
                w = max(min(w,x_max), x_crit); 
                f_continue_preg = linear_interp_w_2_mat( tyear-1, w, f, N, x_crit, x_max, xi);
                f_continue_preg = f_continue_preg(1);
                % abort pregnancy
                f_abort_preg = linear_interp_w_2_loss_mat( tyear-1, w_abt, f, N, x_crit, x_max, xi);
                f_abort_preg = f_abort_preg(1);                
                % now she may be come senescent, with probability ps, and if
                % not, then f results from decision to abort or continue pregnancy:
                f_MC(1,t_breakup,tyear-1) = sigma_hat*((1-ps)*max(f_abort_preg, f_continue_preg).' + ...
                    ps*f_abort_preg_sens.');
                % record which decision is optimal (abort or continue preg.)
                if rand < ps % become senescent
                    dpreg_MC(1,tyear-1) = 1;
                    xi = 2;
                    sensyear = tyear;
                    x = w_abt;
                    eta = 1;
                else
                    if f_abort_preg > f_continue_preg % if aborted pregnancy
                        dpreg_MC(1,tyear-1) = 1;
                        disp('aborted pregnancy')
                        if rand > sigma_hat % die over winter
                            disp('Died over winter from random causes in year')
                            %disp(tyear)
                            break
                        else
                            x = w_abt;
                            eta = 1;
                        end  
                    else 
                        % if continued pregnancy
                        dpreg_MC(1,tyear-1) = 2;
                        disp('continued pregnancy')
                        if rand > sigma_hat % die over winter
                            disp('Died over winter from random causes in year')
                            disp(tyear)
                            break
                        else
                        x = w;
                        eta = 3; 
                        end
                    end
                end
            elseif xi == 2
                f_MC(1,t_breakup,tyear-1) = sigma_hat*f_abort_preg_sens;
                dpreg_MC(1,tyear-1) = 1;
                x = w_abt;
                eta = 1;
            end               
        elseif eta == 3 % with COYS
            % keep litter
            w = w_3_mat( x, tau_icefree, L );
            w = max(min(w,x_max), x_crit); 
            f_keep_litter = linear_interp_w_3_mat( tyear-1, w, f, N, x_crit, x_max, 1);
            f_keep_litter = f_keep_litter(1);
            % abandon litter
            w_abdn = w_3_loss_mat( x, tau_icefree, L );
            w_abdn = max(min(w_abdn,x_max), x_crit); 
            f_abandon_litter = linear_interp_w_3_loss_mat( tyear-1, w_abdn, f, N, x_crit, x_max,1);
            f_abandon_litter = f_abandon_litter(1);            
            % abandon litter because senescent
            f_abandon_litter_sens = linear_interp_w_3_loss_mat( tyear-1, w_abdn, f, N, x_crit, x_max,2);
            f_abandon_litter_sens = f_abandon_litter_sens(1); 
            if xi == 1
                % f results from decision to abort or continue pregnancy:
                f_MC(1,t_breakup,tyear-1) = sigma_hat*((1-ps)*max(f_keep_litter, f_abandon_litter) +...
                    ps*f_abandon_litter_sens);                
                if rand < ps % become senescent and lose litter
                    dCOY_MC(1,tyear-1) = 1;
                    xi = 2;
                    sensyear = tyear;
                    x = w_abdn;
                    eta = 1;
                else
                    % record which decision is optimal (abort or continue preg.)      
                    if f_abandon_litter > f_keep_litter % abandoned litter
                        disp('abandoned litter of COYS')
                        dCOY_MC(1,tyear-1) = 1;
                        if rand > sigma_hat  % die over winter
                            disp('Died over winter from random causes in year')
                            disp(tyear)
                            break
                        else % survive the winter
                            x = w_abdn;
                            eta = 1;
                        end
                    else
                        dCOY_MC(1,tyear-1) = 2;
                        disp('kept litter of COYs')
                        if rand > sigma_hat % die over winter
                            disp('Died over winter from random causes in year')
                            disp(tyear)
                            break
                        else % survive the winter
                            x = w;
                            eta = 4;
                        end
                    end
                end
            elseif xi == 2
                f_MC(1,t_breakup,tyear-1) = sigma_hat*f_abandon_litter_sens;
                x = w_abdn;
                eta = 1;  
            end   
        elseif eta == 4
            num_offspring = num_offspring + k;
            w = w_4_mat( x, tau_icefree, L );
            w = max(min(w,x_max), x_crit); 
            % if not senescent:
            fw4 = linear_interp_w_4_mat(tyear, w, f, N, x_crit, x_max, 1);
            fw4 = fw4(1);
            % if bear becomes senscent:
            fw4_sens = linear_interp_w_4_mat(tyear, w, f, N, x_crit, x_max, 2);
            fw4_sens = fw4_sens(1);             
            % if bear is not senescent: 
            if xi == 1
                f(1,t_breakup,tyear-1,1) = k + sigma_hat*((1-ps)*fw4 + ps*fw4_sens);
                if rand < ps % become senescent (so no more litters)
                    dCOY_MC(1,tyear-1) = 1;
                    xi = 2;
                    sensyear = tyear;
                end
                if rand > sigma_hat % die over winter
                    disp('Died over winter from random causes in year')
                    disp(tyear)
                    break
                else % survive the winter
                x = w;
                eta = 1;
                end
            elseif xi == 2
                % if bear is senescent (although I think this can't happen)
                f(1,t_breakup,tyear-1,2) = k + sigma_hat*fw4_sens;
                x = w;
                eta = 1;
            end
        end      
        MC_x(1,t,tyear) = x;
        MC_eta(1,t,tyear) = eta;
        if x == 0 % if died from random causes, get out of the loop for this individual
            disp('died from random causes')
            break
        end        
        if x == x_crit
            disp('Individual died of starvation over winter in year')
            disp(tyear)
            break % gets us out of while loop
        end
    end % if tyear > 1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    while (t <= t_breakup-1) 
        % calculate state/time dependent parameters
        epsilon = epsilon_F_mat( x, L, x_crit, t, t_breakup, tau_mate);         
        a = Ai_F_mat( x , L );
        Y = Yi_F( t );
        [ g3, g4 ] = milkFun_mat(x, L, x_crit);           
            
        for q = 1:2
            [xp, xpp, x1, xp3, xpp3, xp4, xpp4] = stateChanges_mat(x, a,...
                    Y, q, tau_mate, g3, g4, x_max, x_crit,L);
            
            % obtain fitness values (linear interpolate, if necessary)
                [Value] = ValueFunctions_mat(t, tyear, q, ...
                    f, xp, xpp, x1, xp3, xpp3, xp4, xpp4, epsilon, g3, g4,...
                    N, t_breakup, tau_icefree, T, L, x_crit, x_max, sigma,...
                    sigma_hat, tau_mate, tau_den, k, lambda, sigma_0, sigma_1,xi);
%                 if Value(eta) == 0
%                     break
%                 end
            if q==1  % foraging patch 1, Landfast ice
                vL_MC(1,t,tyear) = Value(eta); % get state specific value
            else  % foraging patch 2, Pack ice
                vP_MC(1,t,tyear) = Value(eta);                    
            end                       
        end
        
        % now find the best patch and save the corresponding decision
        vmax = max(vL_MC(1,t,tyear), vP_MC(1,t,tyear));   
        
        % first check if they're the same
        if vL_MC(1,t,tyear) == vP_MC(1,t,tyear) 
            %disp('Value function is equal in both patches')
            if t == 1
                if tyear == 1
                    disp('Note, arbitrary starting patch as it didnt matter')
                    if rand(1) < 0.5
                        dpatch_MC(1,t,tyear) = 1;
                    else
                        dpatch_MC(1,t,tyear) = 2;
                    end
                else
                    dpatch_MC(1,t,tyear) = dpatch_MC(1, t_breakup-1,tyear-1); % if it's all the same, stay put
                    %disp('Arbitrary; same patch as last step')
                end
            else
                %disp('Arbitrary; same patch as last step')
                if dpatch_MC(1,t-1,tyear) ~= 0
                    dpatch_MC(1,t,tyear) = dpatch_MC(1,t-1,tyear); % if it's all the same, stay put
                else
                    dpatch_MC(1,t,tyear) = dpatch_MC(1,max(t-tau_mate,1),tyear); % in case just got pregnant
                end
            end
        elseif vmax == vL_MC(1,t,tyear)
            dpatch_MC(1,t,tyear) = 1;
        elseif vmax == vP_MC(1,t,tyear)
            dpatch_MC(1,t,tyear) = 2;
        end
        
        % update fitness for current time step 
        f_MC(1,t,tyear) = vmax;      
        
        %%%%%%% STOCHASTIC SIMULATION OF EVENTS %%%%%%%%%%%%%%%%%%%%%%%%%%%
        q = dpatch_MC(1,t,tyear);
        
        if q == 0 % if pregnant, this causes grief here, so fix it
            disp('t is')
            disp(t)
            if t == 1
                if rand(1) < 0.5
                    dpatch_MC(1,t,tyear) = 1;
                else
                    dpatch_MC(1,t,tyear) = 2;
                end
                q = dpatch_MC(1,t,tyear);
            else
                dpatch_MC(1,t,tyear) = dpatch_MC(1,t-1,tyear);
                disp(dpatch_MC(1,t-1,tyear));
                q = dpatch_MC(1,t-1,tyear);
            end
        end
        %disp('q is')       
        %disp(q)
        %disp(lambda(q))
              
        if rand > sigma  % individual doesn't survive
            disp('Individual died from random chance');
            disp(tyear)
            breakloop = 1;
            break            
        else % individual survives
            if MC_eta(1,t,tyear) == 1 % single
                % calculate change in x due to mating
                mate_costs = 0;
                xmate = x;
                for j = 1:tau_mate
                    mate_costs = mate_costs + Ai_F_mat( xmate , L );
                    xmate = xmate - Ai_F_mat( xmate , L ) ;                                   
                end               
                if x-mate_costs <= x_crit % if too thin, ensure don't mate
                    epsilon = 0;
                end
                if rand <= epsilon % mate
                    xm = x - mate_costs;   
                    xm = max(min(xm, x_max),x_crit);
                    for c = (t+1):(t+tau_mate-1)
                        MC_x(1,c,tyear) = MC_x(1,c-1,tyear)-Ai_F_mat(MC_x(1,c-1,tyear),L); % declining mass until last day of mating
                        MC_eta(1,c,tyear) = eta; % same state until last day of mating
                    end
                    MC_x(1,t+tau_mate,tyear) = xm;
                    MC_eta(1,t+tau_mate,tyear) = 2; % now pregnant
                    x = xm;
                    eta = 2;
                    t = t+tau_mate;
                else % don't mate
                    if rand > lambda(q)  % no seal encounter
                        xpp = x - a;      
                        MC_x(1,t+1,tyear) = max(min(xpp, x_max), x_crit);
                        MC_eta(1,t+1,tyear) = 1;
                    else % seal encounter
                        xp = x - a + Y(q);
                        MC_x(1,t+1,tyear) = max(min(xp, x_max), x_crit);
                        MC_eta(1,t+1,tyear) = 1;
                    end
                    if MC_x(1,t+1,tyear) == x_crit
                        disp('Individual died of starvation in year:')
                        disp(tyear)
                        break
                    end
                    x = MC_x(1,t+1,tyear);
                    eta = 1;
                    t=t+1;
                end
                    
            elseif MC_eta(1,t,tyear) == 2 % pregnant %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if rand > lambda(q)  % no seal caught
                    xpp = x - a;      
                    MC_x(1,t+1,tyear) = max(min(xpp, x_max), x_crit);
                    MC_eta(1,t+1,tyear) = 2;
                else % seal caught
                    xp = x - a + Y(q);
                    MC_x(1,t+1,tyear) = max(min(xp, x_max), x_crit);
                    MC_eta(1,t+1,tyear) = 2;
                end
                if MC_x(1,t+1,tyear) == x_crit
                    disp('Individual died of starvation in year:')
                    disp(tyear)
                    break
                end
                x = MC_x(1,t+1,tyear);
                eta = 2;
                t=t+1;
                
            elseif MC_eta(1,t,tyear)==3 % with COYs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if rand > sigma_0(q) % lose litter
                    if rand > lambda(q)  % no seal caught
                        xpp = x - a;      
                        MC_x(1,t+1,tyear) = max(min(xpp, x_max), x_crit);
                        MC_eta(1,t+1,tyear) = 1;
                    else % seal caught
                        xp = x - a + Y(q);
                        MC_x(1,t+1,tyear) = max(min(xp, x_max), x_crit);
                        MC_eta(1,t+1,tyear) = 1;
                    end
                    if MC_x(1,t+1,tyear) == x_crit
                        disp('Individual died of starvation in year:')
                        disp(tyear)
                        break
                    end
                    x = MC_x(1,t+1,tyear);
                    eta = 1;
                    t=t+1;
                    
                else % keep litter
                    if rand > lambda(q)  % no seal caught
                        xpp = x - a - g3;      
                        MC_x(1,t+1,tyear) = max(min(xpp, x_max), x_crit);
                        MC_eta(1,t+1,tyear) = 3;
                    else % seal caught
                        xp = x - a + Y(q) - g3;
                        MC_x(1,t+1,tyear) = max(min(xp, x_max), x_crit);
                        MC_eta(1,t+1,tyear) = 3;
                    end
                    if MC_x(1,t+1,tyear) == x_crit
                        disp('Individual died of starvation in year:')
                        disp(tyear)
                        break
                    end
                    x = MC_x(1,t+1,tyear);
                    eta = 3;
                    t=t+1;
                end
                       
            elseif MC_eta(1,t,tyear) == 4 % with yearlings %%%%%%%%%%%%%%%%%%%%%%%%%%%
                if rand > sigma_1(q) % lose litter
                    if rand > lambda(q)  % no seal caught
                        xpp = x - a;      
                        MC_x(1,t+1,tyear) = max(min(xpp, x_max), x_crit);
                        MC_eta(1,t+1,tyear) = 1;
                    else % seal caught
                        xp = x - a + Y(q);
                        MC_x(1,t+1,tyear) = max(min(xp, x_max), x_crit);
                        MC_eta(1,t+1,tyear) = 1;
                    end
                    if MC_x(1,t+1,tyear) == x_crit
                        disp('Individual died of starvation in year:')
                        disp(tyear)
                        break
                    end
                    x = MC_x(1,t+1,tyear);
                    eta = 1;
                    t=t+1;
                    
                else % keep litter
                    if rand > lambda(q)  % no seal caught
                        xpp = x - a - g4;      
                        MC_x(1,t+1,tyear) = max(min(xpp, x_max), x_crit);
                        MC_eta(1,t+1,tyear) = 4;
                    else % seal caught
                        xp = x - a + Y(q) - g4;
                        MC_x(1,t+1,tyear) = max(min(xp, x_max), x_crit);
                        MC_eta(1,t+1,tyear) = 4;
                    end
                    if MC_x(1,t+1,tyear) == x_crit
                        disp('Individual died of starvation in year:')
                        disp(tyear)
                        break
                    end
                    x = MC_x(1,t+1,tyear);
                    eta = 4;
                    t=t+1;
                end
            end % end eta class statement
        end  % end survival if statement 
        
    end % end while statement
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if tyear == T
        if eta == 4
            f_MC(1, t , tyear) = k;
            num_offspring = num_offspring + k;
        else
            f_MC(1, t, tyear) = 0;
        end
    end    
end % end years loop

disp('simulation over')
disp('number of offspring:')
disp(num_offspring)
disp('senescence occured at age: ')
disp(sensyear-1+4)
end % End function