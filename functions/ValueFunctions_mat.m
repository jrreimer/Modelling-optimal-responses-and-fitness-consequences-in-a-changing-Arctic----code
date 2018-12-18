function[Value] = ValueFunctions_mat(t, tyear, i, ...
    f, xp, xpp, x1, xp3, xpp3, xp4, xpp4, epsilon, g3, g4,...
    N, t_breakup, tau_icefree, T, L, x_crit, x_max, sigma,...
    sigma_hat, tau_mate, tau_den, k, lambda, sigma_0, sigma_1,xi)

                % SINGLE BEAR, ETA = 1
                [ f1_xp, f1_xpp, f2_x1 ] = linear_interp_1_mat( t, tyear, xp,...
                    xpp, x1, f, N, x_crit, x_max, tau_mate, t_breakup, xi);
                
                Value(1,:) = sigma*(... % survive - doesn't depend on patch
                    epsilon.*f2_x1 +... % mate - doesn't depend on patch
                    (1-epsilon).*(... % don't mate 
                    lambda(i)*f1_xp +... % catch seal - depends on patch
                    (1-lambda(i))*f1_xpp)); % don't catch seal - depends on patch
                
                % PREGNANT BEAR, ETA = 2
                [ f2_xp, f2_xpp] = linear_interp_2_mat( t, tyear, xp, xpp,...
                    f, N, x_crit, x_max, xi);
                
                Value(2,:) = sigma*(... % survive
                    lambda(i)*f2_xp +... % catch seal
                    (1-lambda(i))*f2_xpp); % don't catch seal
                
                % BEAR WITH COYS, ETA = 3 
                [ f3_xp3, f3_xpp3, f1_xp, f1_xpp] = linear_interp_3_mat( t, ...
                    tyear, xp, xpp, xp3, xpp3, f, N, x_crit, x_max, xi);
                
                milkprod = g3 > 0; % i.e. able to produce milk
                Value(3,milkprod) = sigma*(... % survive
                        sigma_0(i)*(... % litter survives
                            lambda(i)*f3_xp3(milkprod) +... % catch seal
                            (1-lambda(i))*f3_xpp3(milkprod)) +... % don't catch seal
                        (1-sigma_0(i))*(... % lose litter
                            lambda(i)*f1_xp(milkprod) +... % catch seal
                            (1-lambda(i))*f1_xpp(milkprod))); % don't catch seal
                nomilk = g3 <= 0; % if can't produce milk, litter starves
                    Value(3,nomilk) = sigma*(... % survive
                        lambda(i)*f1_xp(nomilk) +... % catch seal
                        (1-lambda(i))*f1_xpp(nomilk)); % don't catch seal
                                       
                    
                % BEAR WITH YEARLINGS, ETA = 4
                [ f4_xp4, f4_xpp4, f1_xp, f1_xpp] = linear_interp_4_mat( t, ...
                    tyear, xp, xpp, xp4, xpp4, f, N, x_crit, x_max, xi);

                milkprod = g4>0; % able to produce milk
                    Value(4,milkprod) = sigma*(... % survive
                        sigma_1(i)*(... % litter survives
                            lambda(i)*f4_xp4(milkprod) +... % catch seal
                            (1-lambda(i))*f4_xpp4(milkprod)) + ... % don't catch seal
                        (1-sigma_1(i))*(... % lose litter
                            lambda(i)*f1_xp(milkprod) +... % catch seal
                            (1-lambda(i))*f1_xpp(milkprod))); % don't catch seal
                nomilk = g4<=0; % if can't produce milk, assume litter starves
                    Value(4,nomilk) = sigma*(... % survive
                        lambda(i)*f1_xp(nomilk) +... % catch seal
                        (1-lambda(i))*f1_xpp(nomilk)); % don't catch seal
                
                    
end