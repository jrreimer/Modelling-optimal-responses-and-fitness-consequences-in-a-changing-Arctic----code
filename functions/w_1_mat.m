function [ w ] = w_1_mat( x, tau_icefree, L )
% overwinter energy dynamics of a bear who is single at end of spring
    w = x - RMR_mat(x, L, tau_icefree);
end