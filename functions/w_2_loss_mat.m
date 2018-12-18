function [ w ] = w_2_loss_mat( x, tau_icefree, L )
% overwinter energy dynamics of a bear who is pregnant at end of spring but
% loses her litter before denning; as we assume no additional energy is
% invested in the pregnancy until she enters the den, her energetic change
% is the same as if she were not pregnant, i.e.: 

w = w_1_mat( x, tau_icefree, L ); 

end