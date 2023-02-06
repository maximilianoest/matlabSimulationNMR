
dockerIsReady = false;
formatSpec = '%s';
while ~dockerIsReady
    fileID = fopen('readyWaterFile.txt','r');
    information = fscanf(fileID,formatSpec);
    if information == "dockerIsReady"
        dockerIsReady = true;
    end 
    fclose(fileID);
    if ~dockerIsReady
        disp("Docker isn't ready, yet");
        pause(120);
    end
end
