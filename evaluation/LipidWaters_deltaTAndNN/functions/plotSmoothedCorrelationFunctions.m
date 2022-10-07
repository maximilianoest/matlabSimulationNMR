function plottedFigures = plotSmoothedCorrelationFunctions(...
    relevantResults,smoothingStep,figPosAndSize,theta,phi)
plottedFigures = [];
for whichLipid = fieldnames(relevantResults)'
    lipidData = relevantResults.(whichLipid{1});
    fig = initializeFigure('posAndSize',figPosAndSize);
    for datasetNr = 1:length(lipidData)
        dataset = lipidData(datasetNr);
        corrFuncFirstOrder = squeeze(dataset.corrFuncFirstOrder( ...
            theta,phi,:));
        corrFuncSecondOrder = squeeze(dataset.corrFuncSecondOrder( ...
            theta,phi,:));
        corrFuncPart = 1;
        smoothedCorrFuncFirstOrder = [];
        smoothedCorrFuncSecondOrder = [];
        initializeSubplot(fig,length(lipidData),1,datasetNr);
        if datasetNr == 1
            title(sprintf('Lipid: %s, simDur: %.2d, deltaT: %.2d' ...
                ,whichLipid{1},dataset.simulationDurationInS ...
                ,dataset.deltaTInS));
        else
            title(sprintf('SimDur: %.2d, deltaT: %.2d' ...
                ,dataset.simulationDurationInS,dataset.deltaTInS));
        end
        while corrFuncPart <= length(corrFuncFirstOrder) - smoothingStep
            smoothedCorrFuncFirstOrder(end+1) = ...
                mean(corrFuncFirstOrder(corrFuncPart : corrFuncPart ...
                + smoothingStep));  %#ok<AGROW>
            smoothedCorrFuncSecondOrder(end+1) = ...
                mean(corrFuncSecondOrder(corrFuncPart : corrFuncPart ...
                + smoothingStep));  %#ok<AGROW>
            corrFuncPart = corrFuncPart + smoothingStep;
        end
        smoothedCorrFuncFirstOrder(end+1) = ...
            mean(corrFuncFirstOrder(corrFuncPart : end)); %#ok<AGROW>
        smoothedCorrFuncSecondOrder(end+1) = ...
            mean(corrFuncSecondOrder(corrFuncPart : end));  %#ok<AGROW>
        timeAxis = 0:dataset.deltaTInS*smoothingStep:(length( ...
            smoothedCorrFuncFirstOrder)-1)*dataset.deltaTInS*smoothingStep;
        plot(timeAxis,smoothedCorrFuncFirstOrder);
        plot(timeAxis,smoothedCorrFuncSecondOrder);
        legend('First Order', 'Second Order');
    end
    
    plottedFigures(end+1).whichLipid = whichLipid{1}; %#ok<AGROW>
    plottedFigures(end).figure = gcf;
    plottedFigures(end).matlabSimulationDate = ...
        lipidData(1).matlabSimulationDate;
    
end


end