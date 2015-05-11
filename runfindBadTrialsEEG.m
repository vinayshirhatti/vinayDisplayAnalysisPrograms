% run finBadTrialsEEG
% Vinay Shirhatti, 13 Oct 2014
%==========================================================================

clear all; clc;

% Get the details for protocols
[listProtocolName,subjectNames,expDates,protocolNames,stimTypes] = listProtocols;

% Check the OS and set paths accordingly
if ispc
    folderSourceString = 'K:\';
else
    folderSourceString = '/media/store/';
end

% define grid
gridType = 'EEG';


% Choose the protocols (the indices correspond to those in listProtocols.m)
runIndex = 68:70;

% occipitalElec = [9,10,45,46,59,60,63,64];
% parietalElec = [7,8,15,16,19,37,38,51,52];

occipitalElec = 20:24;
% parietalElec = [7,8,15,16,19,37,38,51,52];

% Main loop to determine the bad trials
for i = 1:length(runIndex)
    index = runIndex(i);
    disp(['Find bad trials for index: ' num2str(index)]);
    subjectName = subjectNames{index};
    expDate = expDates{index};
    protocolName = protocolNames{index};
    
    % check these electrodes to decide the overall bad trials
%     checkTheseElectrodes = [occipitalElec,parietalElec];
    checkTheseElectrodes = occipitalElec;
    
    % Set limits/criteria
    threshold = 6;
    maxLimit = 100; minLimit = -100;
    saveData = 1;
    showTrials = 0;
    
    [~,~,~] = findBadTrialsEEG(subjectName,expDate,...
        protocolName,folderSourceString,gridType,...
        checkTheseElectrodes,threshold,maxLimit,minLimit,saveData,showTrials);
end

