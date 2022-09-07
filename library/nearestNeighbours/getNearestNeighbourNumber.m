function nearestNeighbours = getNearestNeighbourNumber( ...
    configuration,numberOfHs,logFilePath)

nearestNeighbours = configuration.nearestNeighbours;
if nearestNeighbours >= numberOfHs
    logMessage(['The number of nearest neighbours is higher than '
        'the number of possible atoms. PLEASE CHECK YOUR CONFIG ' ...
        'FILE!'],path2LogFile);
    error(['The number of nearest neighbours is higher than the ' ...
        'number of possible atoms. Please check your config file!']);
end

logMessage(sprintf(['Analysing %.f nearest neighbours of ' ...
    'overall %.f hydrogen atoms'],nearestNeighbours,numberOfHs) ...
    ,logFilePath,false);

end
