function [ epsilon ] = epsilon_F_mat( x, L, x_crit, t, t_breakup, tau_mate)
    
    epsval = 0.05;
    % if too close to end of spring, can't mate
    if t <= t_breakup-tau_mate
                epsilon = ones(1,length(x))*epsval;
    else
                epsilon = zeros(1,length(x));
    end    
    
end