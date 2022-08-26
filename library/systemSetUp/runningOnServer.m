function runOnServer = runningOnServer()
plattform = computer;
switch plattform
    case {'PCWIN64'}
        runOnServer = 0;
    case {'GLNXA64'}
        runOnServer = 1;
    otherwise
        error('configuration:checkIfRunOnServer', ...
            'plattform: "%s" is not part of plattform list', ...
            plattform);
end

end
