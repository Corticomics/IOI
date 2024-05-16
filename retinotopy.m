function mapData = retinotopy(subjectID, nStimLocs, nRepsPerLoc, stimDurS, baselineDurS, stimSeconds, blSeconds)
    % Intrinsic Signal Imaging Analysis Function
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
    % Use provided durations for stim and baseline
    nStimFrames      = frameRateHz * stimDurS;
    nBlFrames        = frameRateHz * baselineDurS;

    % Extract Imaging Frames
    cd(selectedDirectory);
    trialFrames = zeros(nTrials,1);
    session = zeros(yPix, xPix, framesPerTrial, nTrials);

    % Loop through all trials (all folders) 
    for trial = 1:nTrials
        folder_path = [selectedDirectory, filesep, num2str(trial)];
        file_path = [folder_path, filesep, 'img_00000.bin'];
        fId = fopen(file_path, 'rb');
        if fId == -1
           error(['Could not open the file for trial ', num2str(trial)]);
        end

        image_data = fread(fId, 'uint16');
        image_data = image_data(23:end);
       
        for i = 1:framesPerTrial-1
            startIndx = i * xPix * yPix + 1;
            image_data(startIndx : startIndx + 11) = [];
        end

        current_trial = zeros(yPix, xPix, framesPerTrial);

        for i = 1:framesPerTrial
            startIndx = (i - 1) * xPix * yPix + 1;
            frame_subset = image_data(startIndx : startIndx + xPix*yPix-1);
            frame = reshape(frame_subset, [xPix, yPix]); 
            frame = rot90(frame); % Rotate 90 degrees counterclockwise
            frame = flipud(frame); % Flip upside down
            current_trial(:,:,i) = frame; % Add the frame to the 3d matrix
        end

        trialFrames(trial,:) = size(current_trial,3); 
        session(:, :, :, trial) = current_trial;
    end  

    clear current_trial clear image_data i trialFrames

    % Compute Df/F for stimOn vs. StimOff
    blFrames = (nBlFrames+1)-(blSeconds*frameRateHz) : nBlFrames; 
    stimFrames = nBlFrames+1 : nBlFrames+(stimSeconds*frameRateHz); 

    deltaF = zeros(yPix, xPix, nTrials);

    for trialNum = 1:nTrials
        stimOn  = mean(session(:,:,stimFrames, trialNum),3);  % mean of stim on frames
        stimOff = mean(session(:,:,blFrames, trialNum),3);    % mean of stim off frames
        deltaF(:,:,trialNum) = (stimOn - stimOff) ./ stimOff; % deltaF/F
    end

    % Average The Data Across Trials From Each Stimulus Location
    mapData = zeros(yPix, xPix, nStimLocs);

    for locNum = 1:nStimLocs
        stimArray = zeros(xPix, yPix, nStimLocs);
        n1 =  locNum;
        n2 =  (locNum + 1 * nStimLocs);
        n3 =  (locNum + 2 * nStimLocs);
        n4 =  (locNum + 3 * nStimLocs);
        n5 =  (locNum + 4 * nStimLocs);      
        n6 =  (locNum + 5 * nStimLocs); 
        n7 =  (locNum + 6 * nStimLocs); 
        n8 =  (locNum + 7 * nStimLocs); 
        n9 =  (locNum + 8 * nStimLocs); 
        n10 = (locNum + 9 * nStimLocs);
        n11 = (locNum + 10 * nStimLocs);
        n12 = (locNum + 11 * nStimLocs);
        n13 = (locNum + 12 * nStimLocs);
        n14 = (locNum + 13 * nStimLocs);
        n15 = (locNum + 14 * nStimLocs);
        
        stimArray = cat(3, deltaF(:,:,n1), deltaF(:,:,n2), deltaF(:,:,n3), deltaF(:,:,n4), deltaF(:,:,n5), ...
            deltaF(:,:,n6), deltaF(:,:,n7), deltaF(:,:,n8), deltaF(:,:,n9), deltaF(:,:,n10), ...
            deltaF(:,:,n11), deltaF(:,:,n12), deltaF(:,:,n13), deltaF(:,:,n14), deltaF(:,:,n15));
        
        mapData(:,:,locNum) = mean(stimArray,3);
    end

    % Show Results Along with Average Image of Window
    avgImage = mean(mean(session, 3), 4); 
    avgImage = imgaussfilt(avgImage);

    figure('Position',[10 10 750 1500]);
    hold on;
    for locNum = 1:nStimLocs
        subplot(2, 3, locNum);
        imagesc(imgaussfilt(mapData(:,:,locNum),0.5));
        colormap('gray');
        colorbar;
        title(['Stimulus Location ', num2str(locNum)]);
        xlabel('X-axis'); 
        ylabel('Y-axis');
    end

    subplot(2, 3, 6);
    imagesc(avgImage);
    colormap('gray'); 
    colorbar;
    title('Reference Image');
    xlabel('X-axis');
    ylabel('Y-axis');
    hold off; 

    saveas(gcf, [subjectID, '.png']);
end
