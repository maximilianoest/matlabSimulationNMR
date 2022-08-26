function composingMode = getComposingModeFromFileName(fileName)
fileName = strsplit(fileName,'_');
composingMode = fileName{8};
end
