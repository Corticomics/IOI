% Intrinsic Signal Imaging Analysis Code by AKF and JJC (c) @ 2024

% Steps:
% Converts Binary Files from OiS200 Microscope to unit16
% Extracts Imaging metadata From info.txt
% Deletes Image Header Data
% Extracts Each Frame of the Imaging Session
% Computes change in absorption for stimOn relative to stimOff 

% Next Steps: Convert to function with inputs
% subjectID
% nStimLocs
% nRepsPerLoc
% stimDurS
% baselineDurS
%% Get User Input for location of data files
selectedDirectory = uigetdir('Select the directory containing trial folders');
% If user clicked 'Cancel'
if selectedDirectory == 0
    error('User canceled folder selection');
end
%% Convert to Function Inputs Later
nStimLocs   = 5;
nRepsPerLoc = 15;
% Can also Use Dir to Figure out number of trials = number of folders. 
% S = dir(selectedDirectory);
% nTrials = nnz(~ismember({S.name},{'.','..'})&[S.isdir]);
%% Parse info.txt for relevant info

% All datasets will have at least 1 folder, use info.txt from 1st trial.
cd([selectedDirectory,'/',num2str(1)]);
% read text file as table
t = readtable('info.txt','ReadVariableNames', false);
% Extract Relevant Inputs for Analysis
frameRateHz = table2array(t(1,2));
xPix        = table2array(t(2,2));
yPix        = table2array(t(3,2));
binning     = table2array(t(4,2));
expTimeMS   = table2array(t(5,2));
trialTimeS  = table2array(t(25,2));

% Compute Basics for Analysis
nTrials          = nStimLocs*nRepsPerLoc;
framesPerTrial   = frameRateHz*trialTimeS;
framesPerSession = framesPerTrial*nTrials;
% We used equal durations for stim vs. baseline
nStimFrames       = framesPerTrial/2;
nBlFrames         = framesPerTrial/2;

%% Extract Imaging Frames

% Navigate Back to Session Directory Level
cd(selectedDirectory);
% Initialize a matrix to save the size of each trial
trialFrames = zeros(nTrials,1);
% Initialize the 4D matrix
session = zeros(xPix, yPix, framesPerTrial, nTrials);

% Loop through all trials (all folders) 
for trial = 1:nTrials
    % Construct the folder path for the current trial
    folder_path = [selectedDirectory, filesep, num2str(trial)];
    % Construct the full file path for the current trial
    file_path = [folder_path, filesep, 'img_00000.bin'];
    % Open the binary file for reading
    fId = fopen(file_path, 'rb');
    % Check if the file was successfully opened
    if fId == -1
       error(['Could not open the file for trial ', num2str(trial)]);
    end

    % Read the file
    image_data = fread(fId, 'uint16');

    % Exclude the first 22 values from image_data array (header of first frame)
    image_data = image_data(23:end);
       
    % Loop through "image_data" to delete the header of each subsequent frame
    for i = 1:framesPerTrial-1
        startIndx = i * xPix * yPix + 1;
        image_data(startIndx : startIndx + 11) = []; % Exclude header from image_data
    end

    % Initialize a 3d zero matrix to store all frames from this trial 
    current_trial = zeros(xPix, yPix, framesPerTrial);

    % Loop through all frames and store them in the 3d matrix (current_trial)
    for i = 1:framesPerTrial
        startIndx = (i - 1) * xPix * yPix + 1;
        frame_subset = image_data(startIndx : startIndx + xPix*yPix-1);
        frame = reshape(frame_subset, [xPix, yPix]); 
        current_trial(:,:,i) = frame; % Add the frame to the 3d matrix
    end

    % Save the size of the trial 
    trialFrames(trial,:) = size(current_trial,3); 

    % Concatenate frames along the fourth dimension
    session(:, :, :, trial) = current_trial;
end  

% clean up workspace to free RAM
clear image_data frame_subset frame current_trial;

%% Compute Df/F for stimOn vs. StimOff

% For extracting stim v bl frames
blFrames     = 1:nBlFrames;                  % 5 Second BL  
stimFrames   = nStimFrames+1:framesPerTrial; % 5 Second Stimulus 

% Initialize a 3D matrix to store differences for each trial
deltaF = zeros(xPix,yPix,nTrials);

% Compute Change in Signal for On vs Off frames for each trial
for trialNum = 1:nTrials
    stimOn  = mean(session(:,:,stimFrames, trialNum),3);  % mean of stim on frames
    stimOff = mean(session(:,:,blFrames, trialNum),3);    % mean of stim off frames
    deltaF(:,:,trialNum) = (stimOn - stimOff) ./ stimOff; % deltaF/F
end

%% Average The Data Across Trials From Each Stimulus Location

% mapData: xPix-by-yPix-by-nStimLocs
mapData = zeros(xPix, yPix, nStimLocs);

for locNum = 1:nStimLocs
    % indices of deltaF trials corresponding to this stimulus location
    idx = locNum + (0:nRepsPerLoc-1) * nStimLocs;

    % average across repeats (3rd dim of deltaF is "trial")
    mapData(:,:,locNum) = mean(deltaF(:,:,idx), 3);
end

%% Show Results Along with Average Image of Window

% Average Image Across All Frames
avgImage = mean(mean(session, 3), 4); 

% Show Results in 2x3 Grid
figure('Position',[10 10 1500 750]);
hold on;
for locNum = 1:nStimLocs
    subplot(2, 3, locNum);
    imagesc(deltaF(:,:,locNum));
    colormap;
    colorbar('viridis');
    clim([-0.02 0.01])
    title(['Stimulus Location ', num2str(locNum)]);
    xlabel('X-axis'); % Will need to update these to A-P and M-L
    ylabel('Y-axis');
end

% Add the reference image (avgImage)
subplot(2, 3, 6);
imagesc(avgImage);
colormap; 
colorbar('viridis');
title('Reference Image');
xlabel('X-axis'); % Will need to update these to A-P and M-L
ylabel('Y-axis');
hold off; 

% saveas(gcf, 'ISOI_results.png'); % Will need to update this to include subjectID
