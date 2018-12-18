function [ f3_xp3, f3_xpp3, f1_xp, f1_xpp] = linear_interp_3_mat( t, tyear, ...
    xp, xpp, xp3, xpp3, f, N, x_crit, x_max,xi)
    
    % takes new x values, and gets resulting fitness values
    % linear interpolates if necessary
    
    %put real values into computer values
    xp_c = (xp-x_crit)*N/(x_max - x_crit);
    xpp_c = (xpp-x_crit)*N/(x_max - x_crit);
    xp3_c = (xp3-x_crit)*N/(x_max - x_crit);
    xpp3_c = (xpp3-x_crit)*N/(x_max - x_crit);
            
    %find base integers
    jp = floor(xp_c);
    jpp = floor(xpp_c);
    jp3 = floor(xp3_c);
    jpp3 = floor(xpp3_c);
        
    %define delta x_c
    change_p = ((xp_c - jp).');
    change_pp = ((xpp_c - jpp).');
    change_p3 = ((xp3_c - jp3).');
    change_pp3 = ((xpp3_c - jpp3).');
            
    %linear interpolated values of f
    % KEEP LITTER & CATCH SEAL
    m1 = jp3 == N;
        f3_xp3(m1) = f(N+1,t+1,tyear,3,xi);
    m2 = jp3 ~= N;
        f3_xp3(m2) = (1-change_p3(m2)).*f(jp3(m2)+1,t+1,tyear,3,xi) + change_p3(m2).*f(jp3(m2)+2,t+1,tyear,3,xi);
        
    % KEEP LITTER & DON'T CATCH SEAL
    m3 = jpp3 == N;
        f3_xpp3(m3) = f(N+1,t+1,tyear,3,xi);
    m4 = jpp ~= N;
        f3_xpp3(m4) = (1-change_pp3(m4)).*f(jpp3(m4)+1,t+1,tyear,3,xi) + change_pp3(m4).*f(jpp3(m4)+2,t+1,tyear,3,xi);
    
    
    % LOSE LITTER & CATCH SEAL
    m5 = jp == N;
        f1_xp(m5) = f(N+1,t+1,tyear,1,xi);
    m6 = jp ~= N;
        f1_xp(m6) = (1-change_p(m6)).*f(jp(m6)+1,t+1,tyear,1,xi) + change_p(m6).*f(jp(m6)+2,t+1,tyear,1,xi);
    
    
    % LOSE LITTER & DON'T CATCH SEAL
    m7 = jpp == N;
        f1_xpp(m7) = f(N+1,t+1,tyear,1,xi);
    m8 = jpp ~= N;
        f1_xpp(m8) = (1-change_pp(m8)).*f(jpp(m8)+1,t+1,tyear,1,xi) + change_pp(m8).*f(jpp(m8)+2,t+1,tyear,1,xi);
    
end
