function formOfLayer = getFormOfLayerFromFileName(fileName)
fileName = strsplit(fileName,'_');
formOfLayer = fileName{4};
end
