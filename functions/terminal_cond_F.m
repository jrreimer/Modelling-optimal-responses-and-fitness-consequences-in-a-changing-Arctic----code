function [ f ] = terminal_cond_F( t_breakup, f, N, T)
% sets terminal condition of SDP model; for this model, the terminal
% condition for bears in all states at time T is 0

for i = 2:N+1
   f(i,t_breakup,T,4,1) = 0; 
   f(i,t_breakup,T,4,2) = 0; %senescent bears have same terminal fitness
end

end

