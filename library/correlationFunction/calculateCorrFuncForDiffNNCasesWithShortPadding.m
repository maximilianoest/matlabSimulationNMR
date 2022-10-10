function [sumCorrelationFunction] = ...
    calculateCorrFuncForDiffNNCasesWithShortPadding(sphericalHarmonic ...
    ,nearestNeighbourCases)

% EXPLANATIONS OF CALCULATION:
% The symmetric flag for ifft can't be set here because the imaginary part
% of the correlation function is necessary! Therefore, the return value of
% this function is a complex correlation function.
% The reason why in frequency domain the squared real and imaginary part
% are summed is, because the autocorrelation require the "convolution" (not
% directly the convolution, but more the correlation) of a signal and its 
% complex conjugate. This is a multiplication of the signal
% and its complex conjugate in frequency domain and therefore the real and
% imaginary part is squared and then added up.
%
% fft creates a vector with higher amplitudes => to get a valid result from
% fft you have to devide by the length of the vector
% ifft reduces to the original amplitudes i.e. fft(x) = higher amplitudes ,
% ifft(fft(x)) = original amplitudes => with ifft(fft(x) *fft(y)) you need
% to devide by the length but only once
%
% it is NOT possible to just cut the correlation function with lags to get
% a location dependent correlation function. The reason for that is that
% if there is a lot of movement in space, the spherical harmonics contain
% also information about fluctutations different from the location
% considered. Therefore, the simulation time has to be shortended to get a
% loaction dependent correlation function.
%
% NOTE: There is no offset suppression implemented. This have to be made in
% the postprocessing part.

[~,timeSteps] = size(sphericalHarmonic);
zeroPaddingLength = 2^(nextpow2(timeSteps));

% This is faster than double(fft(...
% double for higher precision in correlation functions
fftSphericalHarmonic = fft(double(sphericalHarmonic),zeroPaddingLength,2);

% This is faster than real^2 + img^2
% due to double, this should not change any precision
correlationFunction = ifft(abs(fftSphericalHarmonic).^2 ,[],2);

clear fftSphericalHarmonic

for nearestNeighbours = nearestNeighbourCases
    fieldName = sprintf('nearestNeighbours%g',nearestNeighbours);
    sumCorrelationFunction.(fieldName) = sum( ...
        correlationFunction(1:nearestNeighbours,1:timeSteps),1)/timeSteps;
end

end

