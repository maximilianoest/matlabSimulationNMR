function waterModel = getWaterModelFromFileName(fileName)

fileName = strsplit(fileName,'_');
waterModel = fileName{3};

end
