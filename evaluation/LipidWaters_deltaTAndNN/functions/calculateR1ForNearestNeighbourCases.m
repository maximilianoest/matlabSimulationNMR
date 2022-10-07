function r1_NN_theta_phi = calculateR1ForNearestNeighbourCases( ...
    corrFuncFirstOrder,corrFuncSecondOrder,omega0,deltaT ...
    ,dipoleDipoleConstant)
r1_NN_theta_phi = [];
for nNCase = fieldnames(corrFuncFirstOrder)'
    corrFuncFirstOrder_theta_phi = corrFuncFirstOrder.(nNCase{1});
    corrFuncSecondOrder_theta_phi = corrFuncSecondOrder.(nNCase{1});
    
    r1_NN_theta_phi(end+1,:,:) = calculateR1ForDifferentAngles( ...
        corrFuncFirstOrder_theta_phi,corrFuncSecondOrder_theta_phi ...
        ,omega0,deltaT,dipoleDipoleConstant); %#ok<AGROW>
 
end

end
