function relevantResults = calculateR1ForDifferentLipids(relevantResults)

for whichLipid = fieldnames(relevantResults)'
    lipidData = relevantResults.(whichLipid{1});
    for datasetNr = 1:length(lipidData)
        dataset = lipidData(datasetNr);
        relevantResults.(whichLipid{1})(datasetNr).r1_NN_theta_phi ...
            = calculateR1ForNearestNeighbourCases( ...
            dataset.shortenedCorrFuncFirstOrder ...
            ,dataset.shortenedCorrFuncSecondOrder ...
            ,dataset.omega0,dataset.deltaTInS ...
            ,dataset.dipoleDipoleConstant);
    end
end
end