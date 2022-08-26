function fileNameArray = getFileNameArray(configuration)
fileNames = getValuesFromStringEnumeration(configuration.fileNames,';' ...
    ,'string');
fileNameArray = string();
for fileNameNr = 1:length(fileNames)
    fileName = fileNames(fileNameNr);
    fileName = strsplit(fileName,'_');
    fileName = fileName(2)+fileName(6);
    if ~contains(fileNameArray,fileName)
        fileNameArray = fileNameArray + fileName;
    end
end
fileNameArray = convertStringsToChars(fileNameArray);
end