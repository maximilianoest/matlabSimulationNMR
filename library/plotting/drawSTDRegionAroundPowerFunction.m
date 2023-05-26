function drawSTDRegionAroundPowerFunction(factor,exponent ...
    ,covarianceMatrix,fieldStrengthAxis,rgbColorArray)

stdMatrix = sqrt(covarianceMatrix);
factor_plusSTD = factor + stdMatrix(1,1);
factor_minusSTD = factor - stdMatrix(1,1);
exponent_plusSTD = exponent + stdMatrix(2,2);
exponent_minusSTD = exponent - stdMatrix(2,2);

powerFunction_plusSTD = factor_plusSTD*fieldStrengthAxis.^( ...
    -exponent_plusSTD);
powerFunction_minusSTD = factor_minusSTD*fieldStrengthAxis.^( ...
    -exponent_minusSTD);
fill([fieldStrengthAxis fliplr(fieldStrengthAxis)] ...
    ,[powerFunction_plusSTD fliplr(powerFunction_minusSTD)] ...
    ,rgbColorArray,'FaceAlpha',0.3,'EdgeColor','none'); 

end