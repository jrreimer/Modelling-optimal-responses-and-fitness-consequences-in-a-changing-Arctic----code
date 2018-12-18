function [ f_wx ] = linear_interp_w_2_mat( tyear, w, f, N, x_crit, x_max,xi)

    % takes new x values (w), and gets resulting fitness values
    % linear interpolates if necessary
    
    Nval = N;
    % takes "new" x values, and gets resulting fitness values
    % linear interpolates if necessary
    x_crit = ones(1,N)*x_crit;
    x_max = ones(1,N)*x_max;
    N = ones(1,N)*Nval;
       
    %put real values into computer values, 0:N
    xval_c = (w-x_crit).*N./(x_max - x_crit);
    
    %find base integers
    jbase = floor(xval_c); 
    
    %define delta x_c
    change = ((xval_c - jbase).');
       
    %linear interpolated values of f
    m1 = jbase == Nval;
        f_wx(m1) = f(Nval+1, 1, tyear+1, 3, xi);
    m2 = jbase ~= Nval;
        f_wx(m2) = (1-change(m2)).*f(jbase(m2)+1,1,tyear+1,3) + ...
               change(m2).*f(jbase(m2)+2,1,tyear+1,3,xi);
    
end
