function waterMoleculesCount = getWaterMoleculesCountFromFileName( ...
    fileName)

fileName = strsplit(fileName,'_');
waterMoleculesCount = fileName{5};

end
