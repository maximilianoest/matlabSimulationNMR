function values = getValuesFromStringEnumeration(Enumeration ...
    ,separator,dataType)

stringEnumeration = num2str(Enumeration);
elementsInEnumeration = split(stringEnumeration,separator);

if isempty(elementsInEnumeration)
   warning('No Field Names where given'); 
end

elementsCount = length(elementsInEnumeration);

for elementNr = 1:elementsCount
    element = elementsInEnumeration{elementNr};
    if element ~= ""
        switch (dataType)
            case "string"
                values(elementNr) = string(element); %#ok<AGROW>
            case "numeric"
                values(elementNr) = str2double(element); %#ok<AGROW>
            otherwise
                warning('No such data type known.')
        end
    end
end

end
