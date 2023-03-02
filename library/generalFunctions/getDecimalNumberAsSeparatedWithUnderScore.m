function numberWithUnderScore = ...
    getDecimalNumberAsSeparatedWithUnderScore(number)

numberString = num2str(number);
numberWithUnderScore = replace(numberString,".","_");


end