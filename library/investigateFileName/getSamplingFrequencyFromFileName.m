function samplingFrequency = getSamplingFrequencyFromFileName(fileName)
fileName = strsplit(fileName,'_');
samplingFrequency = fileName{9};
samplingFrequency = samplingFrequency(3:end-2);
if samplingFrequency(1) == '0'
   samplingFrequency = sprintf('%s.%s',samplingFrequency(1) ...
       ,samplingFrequency(2:end));
end
samplingFrequency = str2double(samplingFrequency);
end
