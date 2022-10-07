function relevantResultsForAllNNCases = ...
    getRelevantResultsForAllNNCases(allMatFilesInResultsDirectory ...
    ,resultsPath,omega0)
relevantResultsForAllNNCases = struct();
for fileCounter = 1:length(allMatFilesInResultsDirectory)
    fileName = allMatFilesInResultsDirectory(fileCounter,1).name;
    if strcmp(fileName,"")
        error("Matfile not found.");
    end
    filePath = sprintf('%s%s%s',resultsPath,filesep,fileName);
    simulationResults = load(filePath);
    relevantResultsForOneLipid = ...
        addMatlabSimulationParametersToRelevantResults(simulationResults);
    
    for nNCase = fieldnames( ...
            simulationResults.sumCorrelationFunctionSaverFirstOrder)'
        [corrFuncFirstOrder.(nNCase{1}),corrFuncSecondOrder.(nNCase{1})] =  ...
            getFirstAndSecondOrderCorrFunctFromMatlabSimResults( ...
            simulationResults,nNCase{1});
    end
    relevantResultsForOneLipid.corrFuncFirstOrder = corrFuncFirstOrder;
    relevantResultsForOneLipid.corrFuncSecondOrder = corrFuncSecondOrder;
    relevantResultsForOneLipid.omega0 = omega0;
    
    relevantResultsForOneLipid.nNCases = ...
        simulationResults.nearestNeighbourCases;
    whichLipid = simulationResults.whichLipid;
    
    if ~isfield(relevantResultsForAllNNCases,whichLipid)
        relevantResultsForAllNNCases.(whichLipid) = ...
            relevantResultsForOneLipid;
    else
        relevantResultsForAllNNCases.(whichLipid)(end+1) = ...
            relevantResultsForOneLipid;
    end
end

end
