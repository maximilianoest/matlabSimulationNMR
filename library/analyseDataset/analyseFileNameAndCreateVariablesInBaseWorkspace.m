function analyseFileNameAndCreateVariablesInBaseWorkspace(configuration ...
    ,logFilePath)

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

samplingFrequency = getSamplingFrequencyFromFileName(fileName);
saveVariableToBaseWorkspace(samplingFrequency);

simulationTime = getSimulationTimeFromFileName(fileName);
saveVariableToBaseWorkspace(simulationTime);

logMessage(sprintf(['Data was simulated with the following ' ...
    'information: \n' ...
    '    GROMACS Simulation Date: %s \n' ...
    '    Lipid that is simulated: %s \n' ...
    '    Used Water Model: %s \n' ...
    '    Layer Form: %s \n' ...
    '    Number of Water Molecules: %s \n' ...
    '    Constituent of Simulation: %s \n' ...
    '    Composing mode in postprocessing: %s \n' ...
    '    Basic Sampling Frequency(changes during sim): %.3f ps \n' ...
    '    Simulation Time: %s ns\n'],simulationDate,whichLipid,waterModel ...
    ,formOfLayer,waterMoleculesCount,constituent,composingMode ...
    ,samplingFrequency,simulationTime),logFilePath,false);

end
