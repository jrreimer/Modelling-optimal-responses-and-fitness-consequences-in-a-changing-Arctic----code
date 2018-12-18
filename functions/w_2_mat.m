function [ w ] = w_2_mat( x, tau_icefree, L, tau_den)
% overwinter energy dynamics of a bear who is pregnant at end of spring and
% carries her litter to term, producing a litter of COYs by the following
% spring

% define denning metabolic rate
 DMR = zeros(1,length(x)); % running total of energy lost in den
 x_den = x - RMR_mat(x, L, tau_icefree); % enter den after summer with reduced x
 
 for j = 1:tau_den
    DMR = DMR + 0.02*mass_mat(x_den,L).^1.09; % update total energy lost daily
    x_den = x_den - 0.02*mass_mat(x_den,L).^1.09; % update energetic state daily
 end

% overwinter energy dynamics
w = x - RMR_mat(x, L, tau_icefree) - DMR;
% note w = x_den. 

end