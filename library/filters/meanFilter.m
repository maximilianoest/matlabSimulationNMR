function meanFilteredArray = meanFilter(arrayToFilter,windowSize)

if ~mod(windowSize,2) 
   error("Window size must be chosen as an odd number.");
end

halfWindowSize = (windowSize-1)/2;
arrayLength = length(arrayToFilter);
meanFilteredArray = zeros(1,arrayLength);

for elementNr = 1:length(arrayToFilter)
    lowerLimit = elementNr - halfWindowSize;
    upperLimit = elementNr + halfWindowSize;
    
    if lowerLimit < 1
        lowerLimit = 1;
    end
    if upperLimit > arrayLength
        upperLimit = arrayLength;
    end
    
    summedValue = sum(arrayToFilter(lowerLimit:upperLimit));
    meanFilteredArray(elementNr) = summedValue ...
        /(upperLimit - lowerLimit +1);
 
end

end