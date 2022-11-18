function [MSE] = calculateMSE(time_profile_est,output)
    total_error = 0;
    for i=1:length(output)
        total_error = total_error + (time_profile_est(i)-output(i))^2;
    end
    MSE = sqrt(total_error/length(output));
end

