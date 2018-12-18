function [ w ] = w_3_loss_mat( x, tau_icefree, L )
% overwinter energy dynamics of a bear who has a litter of COYs at end of 
% spring but abandons them halfway through the summer, ceasing milk
% production.
    
    
    % through first half of summer, produce milk
    % note: round(tau_icefree/2) in case tau_icefree is an odd number
    for j = 1:round(tau_icefree/2) % only for half of summer; assume litter loss is halfway through
        x = x - 0.392*mass_mat(x,L).^0.813 - 0.17*mass_mat(x,L).^0.75;
    end
    
    % through second half of summer, do not produce milk
    for k = (round(tau_icefree/2) + 1):tau_icefree
        x = x - 0.392*mass_mat(x,L).^0.813; 
    end
   
    w = x;
end