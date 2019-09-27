function data_echo = testrun_yarra_GRASP()

%Pass measurementID as argument... 
measurementID = '113';

% Load the paths to required modules 
IncludePathToModules();
workPath = '~/project/temp/basic_grasp';
tempPath = workPath;
outputPath = '~/project/temp/basic_grasp'; 
modeFile = 'modes/GRASP_basic.mode';

directoryContents = dir(workPath);
for i = 1:length(directoryContents)
    if contains(directoryContents(i).name, measurementID)
        fname = directoryContents(i).name;
    end
end

dataEcho = yarra_GRASP_xdse(workPath, fname, outputPath, tempPath, modeFile);
end

