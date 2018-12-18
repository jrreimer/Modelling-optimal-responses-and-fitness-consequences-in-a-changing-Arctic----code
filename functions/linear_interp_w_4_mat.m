function [ f_wx ] = linear_interp_w_4_mat( tyear, w, f, N, x_crit, x_max, xi)

    % litter is weaned and she returns to state 1, single
    f_wx = linear_interp_w_1_mat( tyear, w, f, N, x_crit, x_max, xi); 
    
    
end
