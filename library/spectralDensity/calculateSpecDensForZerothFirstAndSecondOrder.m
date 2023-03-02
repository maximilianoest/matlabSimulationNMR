function [specDensZero,specDensFirst,specDensSecond]  ...
    = calculateSpecDensForZerothFirstAndSecondOrder( ...
    zerothOrderCorrFunc,firstOrderCorrFunc,secondOrderCorrFunc ...
    ,omega0,deltaT,averagingRegion)

corrFuncLength = length(zerothOrderCorrFunc);
specDensZeroAllTimeSteps = 2*(deltaT*cumsum(zerothOrderCorrFunc));
specDensFirstAllTimeSteps = 2*(deltaT*cumsum(firstOrderCorrFunc ...
    .*exp(-1i*omega0*deltaT*(0:corrFuncLength-1))));
specDensSecondAllTimeSteps = 2*(deltaT*cumsum(secondOrderCorrFunc ...
    .*exp(-1i*omega0*deltaT*(0:corrFuncLength-1))));

specDensLength = length(specDensZeroAllTimeSteps);
firstPoint = round(averagingRegion(1) * specDensLength);
lastPoint = round(averagingRegion(2) * specDensLength);
if firstPoint == 0
    firstPoint = 1;    
end

specDensZero = specDensZeroAllTimeSteps(firstPoint:lastPoint);
specDensFirst = specDensFirstAllTimeSteps(firstPoint:lastPoint);
specDensSecond = specDensSecondAllTimeSteps(firstPoint:lastPoint);



end