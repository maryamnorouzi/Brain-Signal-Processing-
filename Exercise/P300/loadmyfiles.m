clc;
clear;

% Define the directory where your .dat files are located
datFilesDir = 'C:\Users\Pankaj\Documents\MARYAM\BCI2000\P300\pxb010\pxb010\';

% Get a list of all .dat files in the directory
datFiles = dir(fullfile(datFilesDir, '*.dat'));

% Loop through each .dat file
for k = 1:length(datFiles)
    % Full path to the current .dat file
    datFileName = fullfile(datFilesDir, datFiles(k).name);
    
    % Use the BCI2000 load function to read the data
    [data, states, parameters] = load_bcidat(datFileName);
    
    % Construct the name of the output .mat file
    [pathstr, name, ext] = fileparts(datFileName);
    matFileName = fullfile(pathstr, [name '.mat']);
    
    % Save the data, states, and parameters to a .mat file
    save(matFileName, 'data', 'states', 'parameters');
end
