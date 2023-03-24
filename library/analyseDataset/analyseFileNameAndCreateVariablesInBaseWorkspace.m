function analyseFileNameAndCreateVariablesInBaseWorkspace(configuration ...
    ,logFilePath)

constants = readConstantsFile('constants.txt');

fileName = configuration.fileName;
saveVariableToBaseWorkspace(fileName)

gromacsSimulationDate = getSimulationDateFromFileName(fileName);
saveVariableToBaseWorkspace(gromacsSimulationDate);

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

deltaTInPs = getSamplingFrequencyFromFileName(fileName);
saveVariableToBaseWorkspace(deltaTInPs);
deltaTInS = deltaTInPs * constants.picoSecond;
saveVariableToBaseWorkspace(deltaTInS);

simulationDurationInNs = getSimulationTimeFromFileName(fileName);
saveVariableToBaseWorkspace(simulationDurationInNs);
simulationDurationInS = simulationDurationInNs*constants.nanoSecond;
saveVariableToBaseWorkspace(simulationDurationInS);

logMessage(sprintf(['Data was simulated with the following ' ...
    'configuration: \n' ...
    '    GROMACS Simulation Date: %s \n' ...
    '    Lipid that is simulated: %s \n' ...
    '    Used Water Model: %s \n' ...
    '    Layer Form: %s \n' ...
    '    Number of Water Molecules: %s \n' ...
    '    Constituent of Simulation: %s \n' ...
    '    Composing mode in postprocessing: %s \n' ...
    '    Time between two time steps: %.3f ps \n' ...
    '    Simulation duration: %.3f ns\n'],gromacsSimulationDate ...
    ,whichLipid,waterModel,formOfLayer,waterMoleculesCount ...
    ,constituent,composingMode,deltaTInPs ...
    ,simulationDurationInNs),logFilePath,false);

end
