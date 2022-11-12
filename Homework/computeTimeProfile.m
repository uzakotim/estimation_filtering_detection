function [y] = computeTimeProfile(u,Ts)
    
    He = tf(1, 1, Ts);  
    time = 1000;
    T = 0:Ts:time-1;
    T_all = 0:Ts:3*time-1;
    tau1 = 10; % for 0 ... 1000
    tau2 = 20; % for 1000 ... 2000
    tau3 = 5;  % for 2000 ... 3000
    variance =0.01;


    tau = tau1;
    sys = tf(1,[tau^3 3*tau^2 3*tau 1]);
    dsys1 = c2d(sys,Ts,'tustin');
    y1 = lsim(dsys1, u, T');
    
    tau = tau2;
    sys = tf(1,[tau^3 3*tau^2 3*tau 1]);
    dsys2 = c2d(sys,Ts,'tustin');
    y2 = lsim(dsys2, u, T');
    y2 = y2+y1(end)*ones(1000,1);
    
    tau = tau3;
    sys = tf(1,[tau^3 3*tau^2 3*tau 1]);
    dsys3 = c2d(sys,Ts,'tustin');
    y3 = lsim(dsys3, u, T');
    y3 = y3+y2(end)*ones(1000,1);
    
    y = cat(1,y1,y2,y3);     %% time-profile
    y = y+lsim(He,sqrt(variance)*randn(3*time,1),T_all); %% noisy time-profile

end

