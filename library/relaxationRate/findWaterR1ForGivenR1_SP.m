function [waterR1,waterMagnetization] ...
    = findWaterR1ForGivenR1_SP(r1_SPArray,r1_SL,r1_water ...
    ,waterFraction,solidLipidFraction,solidProteinFraction ...
    ,exchange_waterToSolid,exchange_lipidToWater,exchange_proteinToWater ...
    ,timeStepCount,deltaT)

initialWaterMagnetization = 0;
initialSLMagnetization = 0;
initialSPMagnetization = 0;

waterR1 = zeros(1,length(r1_SPArray));
timeAxis = 0:deltaT:(timeStepCount - 1)*deltaT;

equilibriumSLMagnetization = solidLipidFraction;
slMagnetization = zeros(1,timeStepCount);

equilibriumSPMagnetization = solidProteinFraction;
spMagnetization = zeros(1,timeStepCount);

equilibriumWaterMagnetization = waterFraction;
waterMagnetization = zeros(1,timeStepCount);

opts = optimset('Display','off');

for r1_SPCounter = 1:length(r1_SPArray)
    r1_SP = r1_SPArray(r1_SPCounter);
    
    slMagnetization(1) = initialSLMagnetization;
    spMagnetization(1) = initialSPMagnetization;
    waterMagnetization(1) = initialWaterMagnetization;
    
    oldStep = 1;
    % solving Bloch equations numerically
    for newStep = 2:timeStepCount
        slMagnetization(newStep) = slMagnetization(oldStep) + deltaT  ...
            * ((equilibriumSLMagnetization - slMagnetization(oldStep)) ...
            *r1_SL - exchange_lipidToWater*slMagnetization(oldStep) ...
            + exchange_waterToSolid*waterMagnetization(oldStep));
        
        spMagnetization(newStep) = spMagnetization(oldStep) + deltaT  ...
            * ((equilibriumSPMagnetization - spMagnetization(oldStep)) ...
            *r1_SP - exchange_proteinToWater*spMagnetization(oldStep) ...
            + exchange_waterToSolid*waterMagnetization(oldStep));
        
        waterMagnetization(newStep) = waterMagnetization(oldStep)  ...
            + deltaT *((equilibriumWaterMagnetization  ...
            - waterMagnetization(oldStep))*r1_water ...
            - exchange_waterToSolid*waterMagnetization(oldStep) ...
            - exchange_waterToSolid*waterMagnetization(oldStep) ...
            + exchange_lipidToWater*slMagnetization(oldStep) ...
            + exchange_proteinToWater*spMagnetization(oldStep));
        oldStep = newStep;
    end
    
    waterR1(r1_SPCounter) = lsqcurvefit(@relaxationFunction ...
        ,1,{timeAxis equilibriumWaterMagnetization} ...
        ,waterMagnetization,[],[],opts);
    
end
end

function magnetization = relaxationFunction(relaxationRate,inputArguments)
magnetization = ...
    inputArguments{2}*(1 - 1*exp(-relaxationRate*inputArguments{1}));
end