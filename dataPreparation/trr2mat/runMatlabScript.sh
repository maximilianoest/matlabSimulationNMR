#!/bin/bash

echo "Sleeping till GROMACS is ready"
sleep 7h
echo "Awakening to run Matlab script"
echo "Running matlab script"
matlab -nosplash -nodisplay -nodesktop -r 'try; createMatFileFromTrrFile; catch e; disp(e.message); disp(e.stack); save(sprintf("%s_Error",datestr(now,"yyyymmdd")),"e");  end; quit;' > Monolayer_wholeDatasetCreation.log