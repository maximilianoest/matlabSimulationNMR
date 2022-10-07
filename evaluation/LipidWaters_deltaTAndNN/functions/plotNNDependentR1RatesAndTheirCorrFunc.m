function plottedFigures = plotNNDependentR1RatesAndTheirCorrFunc( ...
    relevantResults,datasetsToPlot,nearestNeighbourCase,figPosAndSize)
plottedFigures = [];
region = [ 0.2 0.25];
for whichLipid = fieldnames(relevantResults)'
    lipidData = relevantResults.(whichLipid{1});
    fig = initializeFigure('posAndSize',figPosAndSize);
    initializeSubplot(fig,2,2,1:2);
    title(sprintf('Lipid: %s',whichLipid{1}));
    xlabel('Nearest neighours');
    ylabel('R1 [Hz]');
    legendEntries = {};
    for datasetNr = datasetsToPlot
        dataset = lipidData(datasetNr);
        for thetaNr = 1:length(dataset.fibreAnglesTheta)
            for phiNr = 1:length(dataset.fibreAnglesPhi)
                plot(dataset.nearestNeighbourCases ...
                    ,dataset.r1_NN_theta_phi(:,thetaNr,phiNr));
                legendEntries{end+1} = sprintf(['simDur: %.2d, deltaT:' ...
                    ' %.2d, $\\theta$: %.2f, $\\phi$: %.2f'] ...
                    ,dataset.simulationDurationInS ...
                    ,dataset.deltaTInS,dataset.fibreAnglesTheta(thetaNr) ...
                    ,dataset.fibreAnglesPhi(phiNr)); %#ok<AGROW>
            end
        end
    end
    lgd = legend(legendEntries);
    lgd.Location = 'northwest';
    corrFuncNames = ["shortenedCorrFuncFirstOrder" ...
        ,"shortenedCorrFuncSecondOrder"];
    for name = corrFuncNames
        if strcmp(name,"shortenedCorrFuncFirstOrder")
            subplotPosition = 3;
            titleName = sprintf(['Corr. func. first order; Derivation' ...
                ' from mean (Region: %.2f - %.2f)'],region(1), region(2));
            yLabelName = "Frequency";
            binWidth = 1;
            legendVisibility = 'on';
        else
            subplotPosition = 4;
            titleName = sprintf(['Corr. func. second order; Derivation' ...
                ' from mean (Region: %.2f - %.2f)'],region(1), region(2));
            yLabelName = "";
            binWidth = 3;
            legendVisibility = 'off';
        end
        initializeSubplot(fig,2,2,subplotPosition);
        legendEntries = {};
        title(titleName);
        xlabel("Value Corr. Func");
        ylabel(yLabelName);
        for datasetNr = datasetsToPlot
            dataset = lipidData(datasetNr);
            for thetaNr = 1%:length(dataset.fibreAnglesTheta)
                for phiNr = 1:length(dataset.fibreAnglesPhi)
                    corrFuncPart = squeeze(dataset.(name).( ...
                        nearestNeighbourCase)(thetaNr,phiNr ...
                        ,round(0.2*end):round(0.25*end)))';
                    corrFuncPart = corrFuncPart - mean(corrFuncPart);
                    histogram(corrFuncPart,'BinWidth',binWidth);
                    legendEntries{end+1} = sprintf(['simDur: %.2d, deltaT:' ...
                        ' %.2d, $\\theta$: %.2f, $\\phi$: %.2f'] ...
                        ,dataset.simulationDurationInS ...
                        ,dataset.deltaTInS,dataset.fibreAnglesTheta(thetaNr) ...
                        ,dataset.fibreAnglesPhi(phiNr)); %#ok<AGROW>
                end
            end
        end
        lgd = legend(legendEntries);
        lgd.Visible = legendVisibility;
    end
    plottedFigures(end+1).whichLipid = whichLipid{1}; %#ok<AGROW>
    plottedFigures(end).figure = gcf; 
    plottedFigures(end).matlabSimulationDate = lipidData(1).matlabSimulationDate;
end
end
