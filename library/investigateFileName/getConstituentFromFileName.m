function constituent = getConstituentFromFileName(fileName)
fileName = strsplit(fileName,'_');
constituent = fileName{6};
end
