

%% Test 1: return value is struct
configuration = readConfigurationFile("testConfiguration.txt");
assert(isstruct(configuration))

%% Test 2: return right data types
configuration = readConfigurationFile("testConfiguration.txt");
assert(configuration.firstValue == 1, ...
    'Number not read the right way')
assert(strcmp(configuration.secondValue,'asdf'), ...
    'String not read the right way')
assert(strcmp(configuration.thirdValue,''), ...
    'Empty value not read the right way')
assert(strcmp(configuration.fourthValue,'1;2;3;4;5'), ...
    'String enumeration not read the right way')