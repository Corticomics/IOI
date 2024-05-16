function [mapData, session] = retinotopy(subjectID, nStimLocs, nRepsPerLoc, stimDurS, baselineDurS, respWindowS, blWindowS)
    % Retinotopic Mapping Function
    % Inputs:
    % - subjectID: Identifier for the subject
    % - nStimLocs: Number of stimulus locations
    % - nRepsPerLoc: Number of repetitions per location
    % - stimDurS: Duration of stimulus in seconds
    % - baselineDurS: Duration of baseline in seconds
    % - stimSeconds: Number of seconds of the stimulus to be analyze (at the begining  of stimulus)
    % - blSeconds: Number of seconds of the baseline to be analyze (at the end of baseline)
    

    % Get User Input for location of data files
    selectedDirectory = uigetdir('Select the directory containing trial folders');
    % If user clicked 'Cancel'
    if selectedDirectory == 0
        error('User canceled folder selection');
    end


    % Parse info.txt for relevant info
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
    nTrials          = nStimLocs * nRepsPerLoc;
    framesPerTrial   = frameRateHz * trialTimeS;
    framesPerSession = framesPerTrial * nTrials;
    nStimFrames      = frameRateHz * stimDurS;
    nBlFrames        = frameRateHz * baselineDurS;


    % Preallocate Memory for session data
    session = zeros(yPix, xPix, framesPerTrial, nTrials);

  
    % Extract Imaging Frames
    cd(selectedDirectory);

    % Loop through all trials (all folders) 
    for trial = 1:nTrials
        folder_path = [selectedDirectory, filesep, num2str(trial)];
        file_path = [folder_path, filesep, 'img_00000.bin'];
        fId = fopen(file_path, 'rb');
        if fId == -1
           error(['Could not open the file for trial ', num2str(trial)]);
        end
       
        % Read the file, skipping the first 22 values (header of the first frame)
        image_data = fread(fId, 'uint16');
        fclose(fId);
        image_data = image_data(23:end);
       

       % Remove the first 12 elements (header) from each frame 
        for i = 1:framesPerTrial-1
            startIndx = i * xPix * yPix + 1;
            image_data(startIndx : startIndx + 11) = [];
        end
 
        % Initialize a 3D matrix to store all frames from this trial
        current_trial = zeros(yPix, xPix, framesPerTrial);

        % Loop through all frames and store them in the 3D matrix
        for i = 1:framesPerTrial
            startIndx = (i - 1) * xPix * yPix + 1;
            frame_subset = image_data(startIndx : startIndx + xPix*yPix-1);
            frame = reshape(frame_subset, [xPix, yPix]); 
            frame = rot90(frame); % Rotate 90 degrees counterclockwise
            frame = flipud(frame); % Flip upside down
            current_trial(:,:,i) = frame; % Add the frame to the 3d matrix
        end

        session(:, :, :, trial) = current_trial;
    end  

    clear image_data current_trial

    % Compute Df/F for stimOn vs. StimOff
    blFrames = (nBlFrames+1)-(blWindowS*frameRateHz) : nBlFrames; 
    stimFrames = nBlFrames+1 : nBlFrames+(respWindowS*frameRateHz); 

    deltaF = zeros(yPix, xPix, nTrials);

    for trialNum = 1:nTrials
        stimOn  = mean(session(:,:,stimFrames, trialNum),3);  % mean of stim on frames
        stimOff = mean(session(:,:,blFrames, trialNum),3);    % mean of stim off frames
        deltaF(:,:,trialNum) = (stimOn - stimOff) ./ stimOff; % deltaF/F
    end

    % Average the data across trials for each stimulus location
    mapData = zeros(yPix, xPix, nStimLocs);

    for locNum = 1:nStimLocs
        % initialize the master array for this location
        stimArray = zeros(yPix, xPix, nRepsPerLoc);
        for n = 1:nRepsPerLoc 
            % indices for deltaF's that correspond to the current locNum
            idx = (locNum + (n - 1) * nStimLocs);   
            % Combine data for this location
            stimArray(:,:,n) = deltaF(:,:,idx);
        end 
        % Assign Mean of the trials to Output
        mapData(:,:,locNum) = mean(stimArray, 3);
    end
end
