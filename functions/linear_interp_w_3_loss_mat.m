function [ f_wx ] = linear_interp_w_3_loss_mat( tyear, w, f, N, x_crit, x_max, xi)

    % if she loses the litter, state is as though she were single
    f_wx = linear_interp_w_1_mat( tyear, w, f, N, x_crit, x_max, xi); 
    
    
end
