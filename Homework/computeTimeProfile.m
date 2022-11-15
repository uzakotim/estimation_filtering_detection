function [y] = computeTimeProfile(u,Ts)
    
    
    variance =0.01;

    time = 1000;
    T = 0:Ts:time-1;
    T_all = 0:Ts:3*time-1;
    tau1 = 10; % for 0 ... 1000
    tau2 = 20; % for 1000 ... 2000
    tau3 = 5;  % for 2000 ... 3000

    tau = tau1;
    sys = tf(1,[tau^3 3*tau^2 3*tau 1]);
    dsys1 = c2d(sys,Ts,'tustin');
    a = [(8*tau^3+12*tau^2+6*tau+1) (-24*tau^3-12*tau^2+6*tau+3) (24*tau^3-12*tau^2-6*tau+3) (-8*tau^3+12*tau^2-6*tau+1)];
    He = tf(1, a, Ts,'Variable','z^-1');  
    y1 = lsim(dsys1, u, T')+lsim(He,sqrt(variance)*randn(time,1),T');
    
    tau = tau2;
    sys = tf(1,[tau^3 3*tau^2 3*tau 1]);
    dsys2 = c2d(sys,Ts,'tustin');
    a = [(8*tau^3+12*tau^2+6*tau+1) (-24*tau^3-12*tau^2+6*tau+3) (24*tau^3-12*tau^2-6*tau+3) (-8*tau^3+12*tau^2-6*tau+1)];
    He = tf(1, a, Ts,'Variable','z^-1');  
    y2 = lsim(dsys2, u, T')+lsim(He,sqrt(variance)*randn(time,1),T');
    y2 = y2+y1(end)*ones(1000,1);
    
    tau = tau3;
    sys = tf(1,[tau^3 3*tau^2 3*tau 1]);
    dsys3 = c2d(sys,Ts,'tustin');
    a = [(8*tau^3+12*tau^2+6*tau+1) (-24*tau^3-12*tau^2+6*tau+3) (24*tau^3-12*tau^2-6*tau+3) (-8*tau^3+12*tau^2-6*tau+1)];
    He = tf(1, a, Ts,'Variable','z^-1');  
    y3 = lsim(dsys3, u, T')+lsim(He,sqrt(variance)*randn(time,1),T');
    y3 = y3+y2(end)*ones(1000,1);
    
    y = cat(1,y1,y2,y3);     %% time-profile

end

