function plottedFigures = ...
    plotCorrFuncDistributionsAtZero(relevantResults ...
    ,variationsOfCuttedCorrFunc,theta,phi ...
    ,histogramPlotRegion,binWidths,figPosAndSize)
variationsCount = length(variationsOfCuttedCorrFunc);
plottedFigures = [];
for whichLipid = fieldnames(relevantResults)'
    lipidData = relevantResults.(whichLipid{1});
    numberOfDataSetsInLipidData = length(lipidData);
    fig = initializeFigure('posAndSize',figPosAndSize,'legend',false);
    ylabel('Value of Corr. Func.');
    subplotCounter = 1;
    for datasetNr = 1:numberOfDataSetsInLipidData
        dataset = lipidData(datasetNr);
        for variationNr = 1:variationsCount
            correlationFunctionDuration = ...
                variationsOfCuttedCorrFunc(variationNr);
            timeSteps = round(correlationFunctionDuration ...
                /dataset.deltaTInS);
            
            initializeSubplot(fig,numberOfDataSetsInLipidData ...
                ,length(variationsOfCuttedCorrFunc),subplotCounter);
            subplotCounter = subplotCounter + 1;
            lgd = legend();
            if variationNr == 1
                if datasetNr == 1
                    title(sprintf(['Lipid: %s, $\\theta$: %.2f,' ...
                        ' $\\phi$: %.2f simDur: %.2d, deltaT: %.2d'] ...
                        ,whichLipid{1},dataset.fibreAnglesTheta(theta) ...
                        ,dataset.fibreAnglesPhi(phi) ...
                        ,dataset.simulationDurationInS,dataset.deltaTInS));
                    lgd.Visible = 'On';
                else
                    title(sprintf('simDur: %.2d, deltaT: %.2d' ...
                        ,dataset.simulationDurationInS,dataset.deltaTInS));
                    lgd.Visible = 'Off';
                end
            else
                lgd.Visible = 'Off';
            end
            xlabel(sprintf('Corr Func: %.2d' ...
                ,variationsOfCuttedCorrFunc(variationNr)));
            if correlationFunctionDuration > dataset.simulationDurationInS
                plot(1:10,1:10,'r')
                continue
            end
            lastPartCorrFuncFirstOrder = dataset.corrFuncFirstOrder( ...
                theta,phi,round(histogramPlotRegion(1)*timeSteps):round( ...
                histogramPlotRegion(2)*timeSteps));
            lastPartCorrFuncSecondOrder = dataset.corrFuncSecondOrder( ...
                theta,phi,round(histogramPlotRegion(1)*timeSteps):round( ...
                histogramPlotRegion(2)*timeSteps));
            histogram(lastPartCorrFuncFirstOrder,'BinWidth',binWidths(2));
            histogram(lastPartCorrFuncSecondOrder,'BinWidth',binWidths(2));
            
            set(lgd,'String',[{'First Order'} {'Second Order'}])
        end
    end
    plottedFigures(end+1).whichLipid = whichLipid{1}; %#ok<AGROW>
    plottedFigures(end).figure = gcf;
    plottedFigures(end).matlabSimulationDate = ...
        lipidData(1).matlabSimulationDate;
    
end


end