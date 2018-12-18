function [ f4_xp4, f4_xpp4, f1_xp, f1_xpp] = linear_interp_4_mat( t, ...
                    tyear, xp, xpp, xp4, xpp4, f, N, x_crit, x_max, xi)

    % takes new x values, and gets resulting fitness values
    % linear interpolates if necessary
    Nval = N;
    % takes "new" x values, and gets resulting fitness values
    % linear interpolates if necessary
    x_crit = ones(1,N)*x_crit;
    x_max = ones(1,N)*x_max;
    N = ones(1,N)*Nval;
    
    %put real values into computer values
    xp_c = (xp-x_crit).*N./(x_max - x_crit);
    xpp_c = (xpp-x_crit).*N./(x_max - x_crit);
    xp4_c = (xp4-x_crit).*N./(x_max - x_crit);
    xpp4_c = (xpp4-x_crit).*N./(x_max - x_crit);
            
    %find base integers
    jp = floor(xp_c);
    jpp = floor(xpp_c);
    jp4 = floor(xp4_c);
    jpp4 = floor(xpp4_c);
        
    %define delta x_c
    change_p = ((xp_c - jp).');
    change_pp = ((xpp_c - jpp).');
    change_p4 = ((xp4_c - jp4).');
    change_pp4 = ((xpp4_c - jpp4).');
            
    %linear interpolated values of f
    % KEEP LITTER & CATCH SEAL
    m1 = jp4 == Nval;
        f4_xp4(m1) = f(Nval+1,t+1,tyear,4, xi);
    m2 = jp4 ~= Nval;
        f4_xp4(m2) = (1-change_p4(m2)).*f(jp4(m2)+1,t+1,tyear,4,xi) + change_p4(m2).*f(jp4(m2)+2,t+1,tyear,4,xi);
   
    % KEEP LITTER & DON'T CATCH SEAL
    m3 = jpp4 == Nval;
        f4_xpp4(m3) = f(Nval+1,t+1,tyear,4,xi);
   m4 = jpp4 ~= Nval;
        f4_xpp4(m4) = (1-change_pp4(m4)).*f(jpp4(m4)+1,t+1,tyear,4,xi) + change_pp4(m4).*f(jpp4(m4)+2,t+1,tyear,4,xi);
   
    
    % LOSE LITTER & CATCH SEAL
    m5 = jp == Nval;
        f1_xp(m5) = f(Nval+1,t+1,tyear,1,xi);
    m6 = jp ~= Nval;
        f1_xp(m6) = (1-change_p(m6)).*f(jp(m6)+1,t+1,tyear,1,xi) + change_p(m6).*f(jp(m6)+2,t+1,tyear,1,xi);
        
    % LOSE LITTER & DON'T CATCH SEAL
    m7= jpp == Nval;
        f1_xpp(m7) = f(Nval+1,t+1,tyear,1,xi);
    m8 = jpp ~= Nval;
        f1_xpp(m8) = (1-change_pp(m8)).*f(jpp(m8)+1,t+1,tyear,1,xi) + change_pp(m8).*f(jpp(m8)+2,t+1,tyear,1,xi);
 
end
