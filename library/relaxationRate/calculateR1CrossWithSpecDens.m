function r1Cross = calculateR1CrossWithSpecDens(specDensZerothOrder ...
    ,specDensSecondOrder,dipoleDipoleConstant)

r1Cross = dipoleDipoleConstant*(3/4*abs(real(mean( ...
    specDensSecondOrder))) - 1/12*abs(real(mean(specDensZerothOrder))));

end