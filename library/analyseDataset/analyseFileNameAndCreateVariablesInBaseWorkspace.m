function analyseFileNameAndCreateVariablesInBaseWorkspace(configuration ...
    ,logFilePath)

constants = readConstantsFile('constants.txt');

fileName = configuration.fileName;
saveVariableToBaseWorkspace(fileName)

simulationDate = getSimulationDateFromFileName(fileName);
saveVariableToBaseWorkspace(simulationDate);

whichLipid = getLipidNameFromFileName(fileName);
saveVariableToBaseWorkspace(whichLipid);

waterModel = getWaterModelFromFileName(fileName);
saveVariableToBaseWorkspace(waterModel);

formOfLayer = getFormOfLayerFromFileName(fileName);
saveVariableToBaseWorkspace(formOfLayer);

waterMoleculesCount = getWaterMoleculesCountFromFileName(fileName);
saveVariableToBaseWorkspace(waterMoleculesCount);

constituent = getConstituentFromFileName(fileName);
saveVariableToBaseWorkspace(constituent);

composingMode = getComposingModeFromFileName(fileName);
saveVariableToBaseWorkspace(composingMode);

timeBetweenTimeStepsInPs = getSamplingFrequencyFromFileName(fileName);
saveVariableToBaseWorkspace(timeBetweenTimeStepsInPs);
timeBetweenTimeStepsInS = timeBetweenTimeStepsInPs * constants.picoSecond;
saveVariableToBaseWorkspace(timeBetweenTimeStepsInS);

simulationTimeInNs = getSimulationTimeFromFileName(fileName);
saveVariableToBaseWorkspace(simulationTimeInNs);
simulationTimeInS = str2double(simulationTimeInNs)*constants.nanoSecond;
saveVariableToBaseWorkspace(simulationTimeInS);

logMessage(sprintf(['Data was simulated with the following ' ...
    'information: \n' ...
    '    GROMACS Simulation Date: %s \n' ...
    '    Lipid that is simulated: %s \n' ...
    '    Used Water Model: %s \n' ...
    '    Layer Form: %s \n' ...
    '    Number of Water Molecules: %s \n' ...
    '    Constituent of Simulation: %s \n' ...
    '    Composing mode in postprocessing: %s \n' ...
    '    Time between two time steps: %.3f ps \n' ...
    '    Simulation Time: %s ns\n'],simulationDate,whichLipid,waterModel ...
    ,formOfLayer,waterMoleculesCount,constituent,composingMode ...
    ,timeBetweenTimeStepsInPs,simulationTimeInNs),logFilePath,false);

end
