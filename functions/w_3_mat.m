function [ w ] = w_3_mat( x, tau_icefree, L )
% overwinter energy dynamics of a bear who has a litter of COYs at end of
% spring and keeps her litter through to the subsequent spring
    
    xloss = zeros(1,length(x)); 
    
    for j = 1:tau_icefree
        xloss = xloss + 0.392*mass_mat(x,L).^0.813 + 0.17*mass_mat(x,L).^0.75; % xloss + RMR + milk production
        x = x - 0.392*mass_mat(x,L).^0.813 - 0.17*mass_mat(x,L).^0.75;
    end

    w = x; 
    
end