function simTime = getSimulationTimeFromFileName(fileName)
fileName = strsplit(fileName,'_');
simTime = fileName{10};
simTime = simTime(8:end-2);
simTime = num2str(simTime);
end
