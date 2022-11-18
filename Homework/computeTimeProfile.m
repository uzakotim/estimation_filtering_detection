function [y] = computeTimeProfile(u,Ts,tau1,tau2,tau3)
    
    
    variance =0.01;

    time = 1000;
    T = 0:Ts:time-1;

    tau1 = 10;
    tau = tau1;
    b1_ = [1 3 3 1]/(8*tau^3+12*tau^2+6*tau+1);
    a1_ = [(8*tau^3+12*tau^2+6*tau+1) (-24*tau^3-12*tau^2+6*tau+3) (24*tau^3-12*tau^2-6*tau+3) (-8*tau^3+12*tau^2-6*tau+1)]/(8*tau^3+12*tau^2+6*tau+1);
    c1_ = 1/(8*tau^3+12*tau^2+6*tau+1);
    dsys1 = tf(b1_,a1_,Ts,'Variable','z^-1');
    He = tf(c1_, a1_, Ts,'Variable','z^-1');  
    y1 = lsim(dsys1, u(1:1000), T')+lsim(He,sqrt(variance)*randn(time,1),T');
    
    tau2 = 20;
    tau = tau2;
    b2_ = [1 3 3 1]/(8*tau^3+12*tau^2+6*tau+1);
    a2_ = [(8*tau^3+12*tau^2+6*tau+1) (-24*tau^3-12*tau^2+6*tau+3) (24*tau^3-12*tau^2-6*tau+3) (-8*tau^3+12*tau^2-6*tau+1)]/(8*tau^3+12*tau^2+6*tau+1);
    c2_ = 1/(8*tau^3+12*tau^2+6*tau+1);
    dsys2 = tf(b2_,a2_,Ts,'Variable','z^-1');
    He = tf(c2_, a2_, Ts,'Variable','z^-1');  
    y2 = lsim(dsys2, u(1001:2000), T')+lsim(He,sqrt(variance)*randn(time,1),T');
    
    tau3 = 5;
    tau = tau3;
    b3_ = [1 3 3 1]/(8*tau^3+12*tau^2+6*tau+1);
    a3_ = [(8*tau^3+12*tau^2+6*tau+1) (-24*tau^3-12*tau^2+6*tau+3) (24*tau^3-12*tau^2-6*tau+3) (-8*tau^3+12*tau^2-6*tau+1)]/(8*tau^3+12*tau^2+6*tau+1);
    c3_ = 1/(8*tau^3+12*tau^2+6*tau+1);
    
    dsys3 = tf(b3_,a3_,Ts,'Variable','z^-1');
    He = tf(c3_, a3_, Ts,'Variable','z^-1');  
    y3 = lsim(dsys3, u(2001:3000), T')+lsim(He,sqrt(variance)*randn(time,1),T');
    
    y = cat(1,y1,y2,y3);     %% time-profile

end

