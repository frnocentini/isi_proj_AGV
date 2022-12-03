
function meas = getMeasurement(curr_X, sens_pos, w_sens)
    meas = [curr_X(5,1) + w_sens(1,1)^2*randn(1,1); 
            curr_X(4,1) + w_sens(2,1)^2*randn(1,1);
            curr_X(1,1) + w_sens(3,1)^2*randn(1,1);
            sqrt((sens_pos(1,1)-curr_X(1,1))^2 + curr_X(2,1)^2) + w_sens(4,1)^2*randn(1,1);]
end