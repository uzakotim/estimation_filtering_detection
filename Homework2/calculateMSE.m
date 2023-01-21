function [MSE] = calculateMSE(output1,output2)
    total_error = 0;
    for i=1:length(output1)
        total_error = total_error + (output1(i)-output2(i))^2;
    end
    MSE = sqrt(total_error/length(output1));
end