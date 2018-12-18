function[xp, xpp, x1, xp3, xpp3, xp4, xpp4] = stateChanges_mat(x, a, Y, i, ...
    tau_mate, g3, g4, x_max, x_crit, L)

                % GENERAL STATE CHANGE
                xp = x - a + Y(i); % state change if catch a seal
                xpp = x - a;      % state change if no catch
                
                % SINGLE BEAR, ETA = 1  
                mate_costs = 0;
                xmate = x; % initial state if start mating "now"
                for j = 1:tau_mate
                    mate_costs = mate_costs + Ai_F_mat( xmate , L );
                    xmate = xmate - Ai_F_mat( xmate , L ) ;                                   
                end
                x1 = x - mate_costs;% state change over mating event
                
                % PREGNANT BEAR, ETA = 2 - N/A, covered by xp and xpp
                
                % BEAR WITH COYS, ETA = 3
                xp3 = xp - g3; % additional constraints of milk production
                xpp3 = xpp - g3;
                
                % BEAR WITH YEARLINGS, ETA = 4
                xp4 = xp - g4; 
                xpp4 = xpp - g4;
            
                % ensure x_crit <= all new states <= x_max
                xp = max(min(xp, x_max), x_crit);  
                xpp = max(min(xpp, x_max), x_crit);  
                x1 = max(min(x1, x_max), x_crit);
                xp3 = max(min(xp3, x_max), x_crit);
                xpp3 = max(min(xpp3, x_max), x_crit);
                xp4 = max(min(xp4, x_max), x_crit);
                xpp4 = max(min(xpp4, x_max), x_crit);
end
