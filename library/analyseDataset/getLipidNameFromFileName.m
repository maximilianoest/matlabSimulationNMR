function whichLipid = getLipidNameFromFileName(fileName)
fileName = strsplit(fileName,'_');
whichLipid = fileName{2};
end
