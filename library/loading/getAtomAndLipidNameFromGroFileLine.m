function [lipidID,lipidName,atomName] ...
    = getAtomAndLipidNameFromGroFileLine(groFileLine)

splittedGroFileLine = strsplit(groFileLine," ");
lipidNameColumn = splittedGroFileLine{2};
separatedLipidColumn = textscan(lipidNameColumn,'%d%s');
lipidID = separatedLipidColumn{1};
lipidName = string(separatedLipidColumn{2}{1});

atomNameColumn = splittedGroFileLine{3};
atomName = string(atomNameColumn(1));

end