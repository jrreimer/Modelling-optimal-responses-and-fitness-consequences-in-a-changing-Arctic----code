function [ m ] = mass_mat( x, L )
    % estimate the bear's mass (kg) from storage energy (MJ) and length (m) 
    m = (x + 390.53*L^3)/26.14;
end
