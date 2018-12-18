function [ f2_xp, f2_xpp] = linear_interp_2_mat( t, tyear, xp, xpp,...
    f, N, x_crit, x_max, xi)
    
    % takes new x values, and gets resulting fitness values
    % linear interpolates if necessary
    
    %put real values into computer values
    xp_c = (xp-x_crit)*N/(x_max - x_crit);
    xpp_c = (xpp-x_crit)*N/(x_max - x_crit);
            
    %find base integers
    jp = floor(xp_c);
    jpp = floor(xpp_c);
        
    %define delta x_c
    change_p = ((xp_c - jp).');
    change_pp = ((xpp_c - jpp).');
            
    %linear interpolated values of f
    % CATCH SEAL
    m1 = jp == N;
        f2_xp(m1) = f(N+1,t+1,tyear,2,xi);
    m2 = jp ~= N;
        f2_xp(m2) = (1-change_p(m2)).*f(jp(m2)+1,t+1,tyear,2,xi) + change_p(m2).*f(jp(m2)+2,t+1,tyear,2,xi);
    
    
    % DON'T CATCH SEAL
    m3 = jpp == N;
        f2_xpp(m3) = f(N+1,t+1,tyear,2,xi);
    m4 = jpp ~= N;
        f2_xpp(m4) = (1-change_pp(m4)).*f(jpp(m4)+1,t+1,tyear,2,xi) + change_pp(m4).*f(jpp(m4)+2,t+1,tyear,2,xi);
   
end
