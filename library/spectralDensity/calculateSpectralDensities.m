function [spectralDensity1,spectralDensity2] = ...
    calculateSpectralDensities(correlationFunction1,correlationFunction2 ...
    ,omega0,deltaT,lags)
% This function calculates the spectral densities of the given
% correlation functions

spectralDensity1AllTimeSteps = 2*(deltaT*cumsum(correlationFunction1 ... 
        .*exp(-1i*omega0*deltaT*(0:lags-1))));
spectralDensity2AllTimeSteps = 2*(deltaT*cumsum(correlationFunction2 ...
        .*exp(-2i*omega0*deltaT*(0:lags-1))));

spectralDensity1 = spectralDensity1AllTimeSteps(end);
spectralDensity2 = spectralDensity2AllTimeSteps(end);

end

