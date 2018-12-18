function[ g3, g4 ] = milkFun_mat(x, L, x_crit)
      
    % calculate target milk production
    m3 = 0.24*mass_mat(x,L).^0.75;
    m4 = 0.1*mass_mat(x,L).^0.75;
    
    % check if she is in good enough condition to continue milk production
    good3 = x > x_crit + m3;
        g3(good3) = m3(good3);
    toolow3 = x <= x_crit + m3;
        g3(toolow3) = 0;
    
    good4 = x > x_crit + m4;
    g4(good4) = m4(good4);
    toolow4 = x <= x_crit + m4;
    g4(toolow4) = 0;
end