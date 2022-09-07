function nearestNeighbourCases = getNearestNeighbourCases(configuration ...
    ,numberOfHs,logFilePath)
% nearestNeighbourCases = getNearestNeighbourCases(configuration ...
%     ,numberOfHs,logFilePath)
% Gets nearest neighbour cases from configuration file and sorts it from
% large to small values. Additionally, the function writes to the log file.

nearestNeighbourCases = sort(getValuesFromStringEnumeration( ...
    configuration.nearestNeighbourCases,';','numeric'),'descend');

nearestNeighboursString = '';
for nearestNeighbours = nearestNeighbourCases
    if nearestNeighbours >= numberOfHs
        logMessage(['The number of nearest neighbours is higher than '
            'the number of possible atoms. PLEASE CHECK YOUR CONFIG ' ...
            'FILE!'],path2LogFile);
        error(['The number of nearest neighbours is higher than the ' ...
            'number of possible atoms. Please check your config file!']);
    end
    nearestNeighboursString = sprintf('%s %i',nearestNeighboursString ...
        ,nearestNeighbours);
end

logMessage(sprintf(['Analysing [ %s ] nearest neighbour cases of ' ...
    'overall %i hydrogen atoms'],nearestNeighboursString,numberOfHs) ...
    ,logFilePath,false);
end
