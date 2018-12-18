% daily energy requirements 

function [ A ]= Ai_F_mat( x , L )

    % x is the bear's energy reserves    
    Mass = mass_mat( x, L );
    A = 0.0002*Mass.^2.41;
end

