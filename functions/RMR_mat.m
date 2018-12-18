function [ rmr ] = RMR_mat(x,L,tau_icefree)
% calculate the total energetic losses a female bear experiences over a
% summer of length tau_icefree, dependent on initial energy reserves x
% on the first day of summer. 

    rmr = zeros(1,length(x));
    for j = 1:tau_icefree
        rmr = rmr + 0.392*(mass_mat(x,L).^0.813); % running total of losses due to RMR
        x = x - 0.392*(mass_mat(x,L).^0.813); % update x, accounting for these losses
        %rmr = rmr + 0.392*vpa(mass_mat(x,L).^0.813); % running total of losses due to RMR
        %x = x - 0.392*vpa(mass_mat(x,L).^0.813); % update x, accounting for these losses
    end
    
end

