This is the README file to the existing Simulation.

Following information is important:
- the given mainMagnetigField value in the configuration file is only
    important to calculate a rough estimate of R1 and other constants. For
    better results, please use the calculated and summed correlation 
    functions to determine valid results.
- all set up parameters can be adjusted in the configMain.txt
- other scripts can be added to the scripts directory to simulate other 
    parameters. This script name should then be given in the configMain.txt
- nearestNeigbours is used if one case for NN should be used and 
    nearestNeighbourCases is used if multiple NN cases should be considered.
    But be careful: the script file should consider the right parameter
- in some files a RESULTS directory is used. This directory will not be 
    shipped with a github clone but can be found on the server and must be
    downloaded separately (too high data amount)


