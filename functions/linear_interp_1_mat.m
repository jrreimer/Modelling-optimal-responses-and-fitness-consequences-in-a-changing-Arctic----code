function [ f1_xp, f1_xpp, f2_x1 ] = linear_interp_1_mat( t, tyear, xp, xpp,...
    x1, f, N, x_crit, x_max, tau_mate, t_breakup,xi)

    % takes "new" x values, and gets resulting fitness values
    % linear interpolates if necessary
    
    %put real values into computer values ranging from 0 to N
    xp_c = (xp-x_crit)*N/(x_max - x_crit);
    xpp_c = (xpp-x_crit)*N/(x_max - x_crit);
    x1_c = (x1-x_crit)*N/(x_max - x_crit);
            
    %find base integers
    jp = floor(xp_c);
    jpp = floor(xpp_c);
    j1 = floor(x1_c);
        
    %define delta x_c
    change_p = (xp_c - jp).';
    change_pp = (xpp_c - jpp).';
    change_1 = (x1_c - j1).';
            
    %linear interpolated values of f
    % IF BEAR MATES...
        % recall: for mating bears, go up to t+tau rather than t+1   
        if t+tau_mate > t_breakup
            tau_step = 0; % temporary tau_mate to avoid error 
            % when calculate fitness function below; 
            % in reality, if t+tau_mate > t_breakup, epsilon = 0 so mating
            % won't happen
        else
            tau_step = tau_mate;
        end
        
    m1 = j1 == N;
        f2_x1(m1) = f(N+1,t+tau_step,tyear,2,xi); 
    m2 = j1 ~= N;
        f2_x1(m2) = (1-change_1(m2)).*f(j1(m2)+1,t+tau_step,tyear,2,xi) + change_1(m2).*f(j1(m2)+2,t+tau_step,tyear,2,xi);
        
    % DON'T MATE: CATCH SEAL
    m3 = jp == N;
        f1_xp(m3) = f(N+1,t+1,tyear,1,xi);
    m4 = jp ~= N;
        f1_xp(m4) = (1-change_p(m4)).*f(jp(m4)+1,t+1,tyear,1,xi) + (change_p(m4)).*f(jp(m4)+2,t+1,tyear,1,xi);
    
    
    % DON'T MATE, DON'T CATCH SEAL
    m5 = jpp == N;
        f1_xpp(m5) = f(N+1,t+1,tyear,1,xi);
    m6 = jpp ~= N;
        f1_xpp(m6) = (1-change_pp(m6)).*f(jpp(m6)+1,t+1,tyear,1,xi) + (change_pp(m6)).*f(jpp(m6)+2,t+1,tyear,1,xi);
    
end
