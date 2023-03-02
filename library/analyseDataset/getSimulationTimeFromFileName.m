function simTime = getSimulationTimeFromFileName(fileName)
fileName = strsplit(fileName,'_');
simTime = fileName{10};
simTime = simTime(8:end-2);
if simTime(1) == '0'
    simTime = sprintf('%s.%s',simTime(1),simTime(2:end));
end
simTime = str2double(simTime);
end
