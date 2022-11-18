function [y] = computeTimeProfile(u,Ts,dsys1,dsys2,dsys3,He1,He2,He3)
    variance =0.01;

    time = 1000;
    T = 0:Ts:time-1; 
    y1 = lsim(dsys1, u(1:1000), T')+lsim(He1,sqrt(variance)*randn(time,1),T'); 
    y2 = lsim(dsys2, u(1001:2000), T')+lsim(He2,sqrt(variance)*randn(time,1),T'); 
    y3 = lsim(dsys3, u(2001:3000), T')+lsim(He3,sqrt(variance)*randn(time,1),T');
    
    y = cat(1,y1,y2,y3);     %% time-profile

end

