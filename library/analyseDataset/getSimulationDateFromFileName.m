function simulationDate = getSimulationDateFromFileName(fileName)
fileName = strsplit(fileName,'_');
simulationDate = fileName{1};
end
