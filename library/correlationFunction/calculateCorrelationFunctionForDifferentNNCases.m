function [sumCorrelationFunction] = ...
    calculateCorrelationFunctionForDifferentNNCases(sphericalHarmonic ...
    ,numLags,nearestNeighbourCases)

% EXPLANATIONS OF CALCULATION:
% IFFt is slower than FFT, because FFT has a build in function which
% detects whether the input data is real or complex and therefore chooses
% an optimized way to calculate the inverse FT. Because the input data of
% IFFT is real, this data can be seen as conjugate symmetric. Therefore,
% the 'symmetric' flag at IFFT is set.
% The reason why in frequency domain the squared real and imaginary part
% are summed is, because the autocorrelation require the "convolution" (not
% directly the convolution, but more the correlation) of a signal and its 
% complex conjugate. This is a multiplication of the signal
% and its complex conjugate in frequency domain and therefore the real and
% imaginary part summed up.
%
% fft creates a vector with higher amplitudes => to get a valid result from
% fft you have to devide by the length of the vector
% ifft reduces to the original amplitudes i.e. fft(x) = higher amplitudes ,
% ifft(fft(x)) = original amplitudes => with ifft(fft(x) *fft(y)) you need
% to devide by the length but only once
%
% it is not possible to just cut the correlation function with lags to get
% a location dependent correlation function. The reason for that is that
% if there is a lot of movement in space, the spherical harmonics contain
% also information about fluctutations different from the location
% considered. Therefore, the simulation time has to be shortended.
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
correlationFunction = ifft(abs(fftSphericalHarmonic).^2 ...
    ,[],2,'symmetric');

clear fftSphericalHarmonic

for nearestNeighbours = nearestNeighbourCases
    fieldName = sprintf('nearestNeighbours%g',nearestNeighbours);
    sumCorrelationFunction.(fieldName) = sum( ...
        correlationFunction(1:nearestNeighbours,1:numLags),1)/timeSteps;
    % not posssible because arrays do not have the same length...
end

 end



