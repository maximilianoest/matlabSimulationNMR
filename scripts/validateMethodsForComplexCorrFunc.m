%% script start
% In the last analysis I found out that there was a mistake in the
% correlation functions. These have to imaginary to give the right results.
% Nevertheless the imaginary part of the correlation functions seems to be
% very small, so that this won't have a too high impact and the methods
% used (NN-scaling, reduction of simulation duration etc...)maybe still
% work. This should be tested with this script.

logMessage(sprintf('Starting the external script: %s.m', mfilename) ...
    ,logFilePath,true);

%% Preallocation
if ~exist('nearestNeighbourCases','var')
    logMessage('Preallocation of arrays as a first part of the script' ...
        ,logFilePath,true);
    nearestNeighbourCases = getNearestNeighbourCases(configuration ...
        ,numberOfHs,logFilePath);
    nearestNeighbourCasesCount = size(nearestNeighbourCases,2);
    
    [fibreAnglesTheta,fibreAnglesThetaCount,fibreAnglesPhi ...
        ,fibreAnglesPhiCount] = getThetaAndPhiFromConfiguration( ...
        configuration,logFilePath);
    
    atomCounter = 1;
    calculatedAtomIndices = zeros(1,atomsToCalculate);
    calculationSteps = fibreAnglesThetaCount*fibreAnglesPhiCount;
    
    load(sprintf('matFiles%srandomSequenceOfAtoms.mat',filesep));
    logMessage(sprintf('%s:%s', 'Random sequence of atoms: ' ...
        ,sprintf(' %i',randomSequenceOfAtoms(1:15))),logFilePath,false);
    
    nearestNeighboursIDs = zeros(1,max(nearestNeighbourCases));
    nearestNeighbourDistancesPow3 = zeros(max(nearestNeighbourCases) ...
        ,timeSteps,'single');
    
    rotatedX = zeros(max(nearestNeighbourCases),timeSteps,'single');
    rotatedY = zeros(max(nearestNeighbourCases),timeSteps,'single');
    rotatedZ = zeros(max(nearestNeighbourCases),timeSteps,'single');
    
    polarAngle = zeros(max(nearestNeighbourCases),timeSteps,'single');
    azimuthAngle = zeros(max(nearestNeighbourCases),timeSteps,'single');
    
    sphericalHarmonicZerothOrder = zeros(max(nearestNeighbourCases) ...
        ,timeSteps,'like',single(1j));
    sphericalHarmonicFirstOrder = zeros(max(nearestNeighbourCases) ...
        ,timeSteps,'like',single(1j));
    sphericalHarmonicSecondOrder = zeros(max(nearestNeighbourCases) ...
        ,timeSteps,'like',single(1j));
    
    for nearestNeighbours = nearestNeighbourCases
        key = sprintf('nearestNeighbours%g',nearestNeighbours);
        
        r1Estimation_theta_phi_atomCounter.(key) = zeros( ...
            fibreAnglesThetaCount,fibreAnglesPhiCount,atomsToCalculate);
        r1Estimation_theta_phi.(key) = zeros( ...
            fibreAnglesThetaCount,fibreAnglesPhiCount);
        
        correlationFunctionsZerothOrder.(key) = zeros(1,timeSteps ...
            ,'like',double(1j));
        correlationFunctionsFirstOrder.(key) = zeros(1,timeSteps ...
            ,'like',double(1j));
        correlationFunctionsSecondOrder.(key) = zeros(1,timeSteps ...
            ,'like',double(1j));
        
        sumCorrelationFunctionsSaverZerothOrder.(key) = zeros( ...
            fibreAnglesThetaCount,fibreAnglesPhiCount,timeSteps,'like' ...
            ,double(1j));
        sumCorrelationFunctionSaverFirstOrder.(key) = zeros( ...
            fibreAnglesThetaCount,fibreAnglesPhiCount,timeSteps,'like' ...
            ,double(1j));
        sumCorrelationFunctionSaverSecondOrder.(key) = zeros( ...
            fibreAnglesThetaCount,fibreAnglesPhiCount,timeSteps,'like' ...
            ,double(1j));
    end
    logMessage('Arrays successfully preallocated',logFilePath,false);
    return;
end
meanPositions = single([mean(trajectoryX,2) mean(trajectoryY,2) ...
        mean(trajectoryZ,2)]);
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
        findNearestNeighboursIDs(max(nearestNeighbourCases),trajectoryX ...
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
                
                correlationFunctionsZerothOrder ...
                    = calculateCorrelationFunctionForDifferentNNCases( ...
                    sphericalHarmonicZerothOrder,nearestNeighbourCases);
                sumCorrelationFunctionsSaverZerothOrder ...
                    = addTwoCorrelationFunctionStructs( ...
                    sumCorrelationFunctionsSaverZerothOrder ...
                    ,correlationFunctionsZerothOrder,thetaNumber ...
                    ,phiNumber);
                
                correlationFunctionsFirstOrder ...
                    = calculateCorrelationFunctionForDifferentNNCases( ...
                    sphericalHarmonicFirstOrder,nearestNeighbourCases);
                sumCorrelationFunctionSaverFirstOrder ...
                    = addTwoCorrelationFunctionStructs( ...
                    sumCorrelationFunctionSaverFirstOrder ...
                    ,correlationFunctionsFirstOrder,thetaNumber ...
                    ,phiNumber);
                
                correlationFunctionsSecondOrder ...
                    = calculateCorrelationFunctionForDifferentNNCases( ...
                    sphericalHarmonicSecondOrder,nearestNeighbourCases);
                sumCorrelationFunctionSaverSecondOrder ...
                    = addTwoCorrelationFunctionStructs( ...
                    sumCorrelationFunctionSaverSecondOrder ...
                    ,correlationFunctionsSecondOrder,thetaNumber ...
                    ,phiNumber);
                
                
                for fieldName = string(fieldnames( ...
                        correlationFunctionsFirstOrder))'
                    correlationFunctionFirstOrder = ...
                        correlationFunctionsFirstOrder.(fieldName);
                    correlationFunctionSecondOrder = ...
                        correlationFunctionsSecondOrder.(fieldName);
                    
                    [averageSpectralDensityFirstOrder ...
                        ,averageSpectralDensitySecondOrder] ...
                        = calculateSpectralDensities( ...
                        squeeze(sumCorrelationFunctionSaverFirstOrder...
                        .(fieldName)(thetaNumber,phiNumber,:))' ...
                        /atomCounter,squeeze( ...
                        sumCorrelationFunctionSaverSecondOrder ...
                        .(fieldName)(thetaNumber,phiNumber,:))' ...
                        /atomCounter,omega0,deltaTInS,timeSteps);
                    
                    r1Estimation_theta_phi.(fieldName)( ...
                        thetaNumber,phiNumber) ...
                        = calculateR1WithSpectralDensity( ...
                        averageSpectralDensityFirstOrder ...
                        ,averageSpectralDensitySecondOrder ...
                        ,dipolDipolConstant);
                    
                    logMessage(sprintf(['=> %s, theta: %.2f, phi: %.2f\n' ...
                        '    Estimated relaxation rate for calculated' ...
                        ' atoms: %.4f\n'],fieldName,rad2deg(theta) ...
                        ,rad2deg(phi),r1Estimation_theta_phi.( ...
                        fieldName)(thetaNumber,phiNumber)),logFilePath ...
                        ,false);
                end
                thetaNullCalculated = true;
            else
                for fieldName = string(fieldnames( ...
                        correlationFunctionsFirstOrder))'
                    sumCorrelationFunctionsSaverZerothOrder.(fieldName)( ...
                        thetaNumber,phiNumber,:)  ...
                        = sumCorrelationFunctionsSaverZerothOrder ...
                        .(fieldName)(1,1,:);
                    sumCorrelationFunctionSaverFirstOrder.(fieldName)( ...
                        thetaNumber,phiNumber,:) ...
                        = sumCorrelationFunctionSaverFirstOrder ...
                        .(fieldName)(1,1,:);
                    sumCorrelationFunctionSaverSecondOrder.(fieldName)( ...
                        thetaNumber,phiNumber,:) ...
                        = sumCorrelationFunctionSaverSecondOrder ...
                        .(fieldName)(1,1,:);
                    r1Estimation_theta_phi.(fieldName)( ...
                        thetaNumber,phiNumber) = r1Estimation_theta_phi ...
                        .(fieldName)(1,1);
                    logMessage(sprintf(['=> theta: %s, %.2f, phi: %.2f\n' ...
                        '    Estimated relaxation rate for calculated atoms:' ...
                        ' %.4f\n'],fieldName,rad2deg(theta),rad2deg(phi) ...
                        ,r1Estimation_theta_phi.(fieldName)(thetaNumber ...
                        ,phiNumber)),logFilePath,false);
                    logMessage(sprintf(['R1(theta=0,phi=%.2f) used from' ...
                        'theta = 0 and phi = 0'],rad2deg(phi)) ...
                        ,logFilePath,false);
                end
            end
            printDottedBreakLineToLogFile(logFilePath);
        end
    end
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



