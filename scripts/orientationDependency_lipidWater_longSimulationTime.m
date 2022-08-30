%% script start
% some information:
% due to the long simulation time and the fast movement of the hydrogen
% atoms in the water molecules, the trajectories of these molecules will
% reach across the whole water part of the lipid model. Therefore it won't
% be usefull to determine the location dependency here.

logMessage(sprintf('Starting the external script: %s.m', mfilename) ...
    ,logFilePath,true);

%% simulation specific information

lags = round(configuration.fractionForLags*timeSteps);
logMessage(sprintf(['The lag is set to %d time steps, which ' ...
    'is equivalent to %.2f %%. This configuration only shortens the '...
    'correlation functions and NOT the simulation time.'],lags ...
    ,(configuration.fractionForLags)*100),logFilePath,false);


nearestNeighbours = getNumberOfNearestNeighbourNumber(configuration ...
    ,numberOfHs,logFilePath);

[fibreAnglesTheta,fibreAnglesThetaCount,fibreAnglesPhi ...
    ,fibreAnglesPhiCount] = getThetaAndPhiFromConfiguration( ...
    configuration,logFilePath);

%% Preallocation

logMessage('Preallocation of arrays.',logFilePath,true);

r1Estimation_theta_phi_atomCounter = zeros(fibreAnglesThetaCount ...
    ,fibreAnglesPhiCount,atomsToCalculate);
r1Estimation_theta_phi = zeros(fibreAnglesThetaCount,fibreAnglesPhiCount);
meanPositions = single([mean(trajectoryX,2) mean(trajectoryY,2) ...
    mean(trajectoryZ,2)]);

atomCounter = 1;
calculatedAtomIndices = zeros(1,atomsToCalculate);
calculationSteps = fibreAnglesThetaCount*fibreAnglesPhiCount;

randomSequenceOfAtoms = randperm(numberOfHs);

nearestNeighboursIDs = zeros(1,nearestNeighbours);
nearestNeighbourDistancesPow3 = zeros(nearestNeighbours,timeSteps ...
    ,'single');

rotatedX = zeros(nearestNeighbours,timeSteps,'single');
rotatedY = zeros(nearestNeighbours,timeSteps,'single');
rotatedZ = zeros(nearestNeighbours,timeSteps,'single');

% THINK AGAIN: first transform to spherical harmonics and then rotate by
% add the thetas and phis. Calculate the commutators to determine wether
% this is possible or not.
polarAngle = zeros(nearestNeighbours,timeSteps,'single');
azimuthAngle = zeros(nearestNeighbours,timeSteps,'single');

sphericalHarmonicZerothOrder = zeros(nearestNeighbours,timeSteps ...
    ,'like',single(1j));
sphericalHarmonicFirstOrder = zeros(nearestNeighbours,timeSteps ...
    ,'like',single(1j));
sphericalHarmonicSecondOrder = zeros(nearestNeighbours,timeSteps ...
    ,'like',single(1j));

% THINK AGAIN: maybe try different NN here by using different sums
correlationFunctionZerothOrder = zeros(1,lags,'like',single(1j));
correlationFunctionFirstOrder = zeros(1,lags,'like',single(1j));
correlationFunctionSecondOrder = zeros(1,lags,'like',single(1j));

sumCorrelationFunctionSaverZerothOrder = zeros(fibreAnglesThetaCount ...
    ,fibreAnglesPhiCount,lags,'like',single(1i));
sumCorrelationFunctionSaverFirstOrder = zeros(fibreAnglesThetaCount ...
    ,fibreAnglesPhiCount,lags,'like',single(1i));
sumCorrelationFunctionSaverSecondOrder = zeros(fibreAnglesThetaCount ...
    ,fibreAnglesPhiCount,lags,'like',single(1i));

atomTimer = zeros(1,atomsToCalculate);
%% Start simulation
printBreakLineToLogFile(logFilePath);
logMessage('Start simulation.',logFilePath,true);
for atomNumber = randomSequenceOfAtoms(atomCounter:atomsToCalculate)
    atomTimerStart = tic;
    logMessage(sprintf('Selected atom number %i',atomNumber),logFilePath ...
        ,false);
    calculatedAtomIndices(atomCounter) = atomNumber;
    
    % this seems to be a bit faster than writing the relative coordinates
    % to another variable and additionally saves RAM
    [trajectoryX,trajectoryY,trajectoryZ] = calculateRelativePositions( ...
        trajectoryX,trajectoryY,trajectoryZ,atomNumber);
    
    [nearestNeighboursIDs,nearestNeighbourDistancesPow3] =  ...
        findNearestNeighboursIDs(nearestNeighbours,trajectoryX ...
        ,trajectoryY,trajectoryZ);
    thetaNullCalculated = false;
    for phiNumber = 1:fibreAnglesPhiCount
        phi = fibreAnglesPhi(phiNumber);
        zAxis = [0 0 1];
        rotationMatrixPhi = get3DRotationMatrix(phi,zAxis);
        for thetaNumber = 1:fibreAnglesThetaCount 
            theta = fibreAnglesTheta(thetaNumber);
            if ~thetaNullCalculated || theta > 0
                yAxis = [0 1 0];
                rotationMatrixTheta = get3DRotationMatrix(theta,yAxis);
                totalRotationMatrix = ...
                    rotationMatrixTheta*rotationMatrixPhi;
                [rotatedX,rotatedY,rotatedZ]  ...
                    = rotateTrajectoriesWithRotationMatrix( ...
                    totalRotationMatrix,trajectoryX(nearestNeighboursIDs,:) ...
                    ,trajectoryY(nearestNeighboursIDs,:) ...
                    ,trajectoryZ(nearestNeighboursIDs,:));
                
                [polarAngle,azimuthAngle] = ...
                    transformToSphericalCoordinates(rotatedX,rotatedY ...
                    ,rotatedZ);
                
                [sphericalHarmonicZerothOrder,sphericalHarmonicFirstOrder ...
                    ,sphericalHarmonicSecondOrder] ...
                    = calculateSphericalHarmonics(polarAngle,azimuthAngle ...
                    ,nearestNeighbourDistancesPow3);
                
                correlationFunctionZerothOrder ...
                    = calculateCorrelationFunction( ...
                    sphericalHarmonicZerothOrder,lags);
                sumCorrelationFunctionSaverZerothOrder(thetaNumber,phiNumber,:) ...
                    = squeeze(sumCorrelationFunctionSaverZerothOrder( ...
                    thetaNumber,phiNumber,:))' ...
                    + correlationFunctionZerothOrder;
                
                correlationFunctionFirstOrder ...
                    = calculateCorrelationFunction( ...
                    sphericalHarmonicFirstOrder,lags);
                sumCorrelationFunctionSaverFirstOrder(thetaNumber,phiNumber,:) ...
                    = squeeze(sumCorrelationFunctionSaverFirstOrder( ...
                    thetaNumber,phiNumber,:))' + correlationFunctionFirstOrder;
                
                correlationFunctionSecondOrder ...
                    = calculateCorrelationFunction( ...
                    sphericalHarmonicSecondOrder,lags);
                sumCorrelationFunctionSaverSecondOrder(thetaNumber,phiNumber,:) ...
                    = squeeze(sumCorrelationFunctionSaverSecondOrder( ...
                    thetaNumber,phiNumber,:))' ...
                    + correlationFunctionSecondOrder;
                
                [spectralDensityFirstOrder,spectralDensitySecondOrder] = ...
                    calculateSpectralDensities(correlationFunctionFirstOrder ...
                    ,correlationFunctionSecondOrder,omega0,deltaTInS ...
                    ,lags);
                
                r1Estimation_theta_phi_atomCounter(thetaNumber,phiNumber ...
                    ,atomCounter) = ...
                    calculateR1WithSpectralDensity(spectralDensityFirstOrder ...
                    ,spectralDensitySecondOrder,dipolDipolConstant);
                
                [averageSpectralDensityFirstOrder ...
                    ,averageSpectralDensitySecondOrder] ...
                    = calculateSpectralDensities( ...
                    squeeze(sumCorrelationFunctionSaverFirstOrder(thetaNumber ...
                    ,phiNumber,:))'/atomCounter ...
                    ,squeeze(sumCorrelationFunctionSaverSecondOrder(thetaNumber ...
                    ,phiNumber,:))'/atomCounter,omega0,deltaTInS,lags);
                
                r1Estimation_theta_phi(thetaNumber,phiNumber) ...
                    = calculateR1WithSpectralDensity( ...
                    averageSpectralDensityFirstOrder ...
                    ,averageSpectralDensitySecondOrder,dipolDipolConstant);
                
                logMessage(sprintf(['=> theta: %i, phi: %i\n' ...
                    '    Estimated relaxation rate for atom: %.4f\n'...
                    '    Estimated relaxation rate for calculated atoms:' ...
                    ' %.4f\n'],rad2deg(theta),rad2deg(phi) ...
                    ,r1Estimation_theta_phi_atomCounter(thetaNumber,phiNumber ...
                    ,atomCounter),r1Estimation_theta_phi(thetaNumber ...
                    ,phiNumber)),logFilePath,false);
                thetaNullCalculated = true;
            else
                sumCorrelationFunctionSaverZerothOrder(thetaNumber ...
                    ,phiNumber,:)  ...
                    = sumCorrelationFunctionSaverZerothOrder(1,1,:);
                sumCorrelationFunctionSaverFirstOrder(thetaNumber ...
                    ,phiNumber,:) ...
                    = sumCorrelationFunctionSaverFirstOrder(1,1,:);
                sumCorrelationFunctionSaverSecondOrder(thetaNumber ...
                    ,phiNumber,:) ...
                    = sumCorrelationFunctionSaverSecondOrder(1,1,:);
                r1Estimation_theta_phi_atomCounter(thetaNumber,phiNumber ...
                    ,atomCounter) = ...
                    r1Estimation_theta_phi_atomCounter(1,1 ...
                    ,atomCounter);
                r1Estimation_theta_phi(thetaNumber,phiNumber) = ...
                    r1Estimation_theta_phi(1,1);
                logMessage(sprintf(['R1(tehta=0,phi=%.2f) is used from theta = 0' ...
                    ' and phi = 0, because R1 is the same for different phi at' ...
                    ' theta = 0'],phi),logFilePath,false);
                logMessage(sprintf(['=> theta: %i, phi: %i\n' ...
                    '    Estimated relaxation rate for atom: %.4f\n'...
                    '    Estimated relaxation rate for calculated atoms:' ...
                    ' %.4f\n'],rad2deg(theta),rad2deg(phi) ...
                    ,r1Estimation_theta_phi_atomCounter(thetaNumber,phiNumber ...
                    ,atomCounter),r1Estimation_theta_phi(thetaNumber ...
                    ,phiNumber)),logFilePath,false);
            end
        end
    end
    printDottedBreakLineToLogFile(logFilePath);
    if mod(atomCounter,configuration.savingIntervall) == 0
        lastSavingDate = datestr(now,'yyyymmdd_HHMM');
        createDataSavingObject();
        save(resultsFileSavingPath,'-struct','dataSavingObject','-v7.3');
        logMessage('Saved data',logFilePath);
    end
    logMessage(sprintf('Calculated %i atom(s)',atomCounter),logFilePath ...
        ,false);
    atomTimer(atomCounter) = toc(atomTimerStart);
    averageTimeForOneAtom = seconds(mean(atomTimer(1:atomCounter)));
    logMessage(sprintf([' ---> Average time for one atom: %s \n' ...
        '        Average time for one phi/theta: %s \n'...
        '        Approximately ready on: %s.'] ...
        ,datestr(averageTimeForOneAtom,'HH:MM:SS') ...
        ,datestr(averageTimeForOneAtom/calculationSteps,'HH:MM:SS') ...
        ,datetime('now')+averageTimeForOneAtom ...
        *(atomsToCalculate-atomCounter)),logFilePath,false);
    printEqualSignBreakLineToLogFile(logFilePath);
    atomCounter = atomCounter + 1;
end

