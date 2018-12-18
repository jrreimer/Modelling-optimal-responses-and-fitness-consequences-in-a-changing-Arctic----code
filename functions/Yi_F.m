% Yi values used in SDP model; depends on date in the feeding season
% returns a vector Y = [y fast ice, y pack ice]

function [ Y ]= Yi_F( t )
% t ranges from 1 to tmax-1

    %%% Pack ice %%%
    val = [165908, 268399, 392935];    
    val = val/1000; % convert KJ to MJ
    slope1 = (val(2)-val(1))/30; % Apr. 1-30
    slope2 = (val(3)-val(2))/15; % May 1-15

   if t <= 30
    yp = slope1*(t-1) + val(1);
   elseif t <= 45
    yp = slope2*(t-31) + val(2);
   else
    yp = val(3);
   end

%%% Fast ice %%%
val = [219159, 295744, 394909];
val = val/1000; % convert KJ to MJ
slope1 = (val(2)-val(1))/30; % April 1-30
slope2 = (val(3)-val(2))/15; % May 1-15

if t <= 30
    yf = slope1*(t-1) + val(1);
elseif t <= 45
    yf = slope2*(t-31) + val(2);
else
    yf = val(3);
end

Y = [yf yp];

end
