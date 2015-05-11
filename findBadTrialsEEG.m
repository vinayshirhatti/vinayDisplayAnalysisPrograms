% find bad trials for each electrode for EEG
% modified from findBadTrialsWithLFP and runFindBadTrials
%
% Vinay Shirhatti, 12 Oct 2014
%
% 24 March 2015, modified from findBadTrialsEEG in
% GammaStimDiscontinuityProject
%==========================================================================


function [allBadTrials, badTrials, nameElec] = findBadTrialsEEG(subjectName,expDate,protocolName,folderSourceString,gridType,...
    checkTheseElectrodes,threshold,maxLimit,minLimit,saveData,showTrials)

if ~exist('checkTheseElectrodes','var')     checkTheseElectrodes = 1:64;   end
if ~exist('folderSourceString','var')       folderSourceString = 'K:\';                end
if ~exist('threshold','var')                threshold = 6;                             end
if ~exist('minLimit','var')                 minLimit = -300;                          end
if ~exist('maxLimit','var')                 maxLimit = 300;                          end

folderName = fullfile(folderSourceString,'data',subjectName,gridType,expDate,protocolName);
folderSegment = fullfile(folderName,'segmentedData');
folderLFP = fullfile(folderSegment,'LFP');

load(fullfile(folderLFP,'lfpInfo'));

numElectrodes = length(analogChannelsStored);

numCheckElectrodes = length(checkTheseElectrodes);

allBadTrials = cell(1,numElectrodes);
nameElec = cell(1,numElectrodes);

for i=1:numElectrodes
    clear analogData
    electrodeNum=analogChannelsStored(i); % changed from checkTheseElectrodes
    % to calculate bad trials for each electrode irrespective of the
    % electrodes to be checked
    disp(['Processing electrode: ' num2str(electrodeNum)]);
    nameElec{i} = ['elec' num2str(electrodeNum)];
    disp(nameElec{i});
    
    % Set the limits higher for frontal electrodes 1,2 & 61(for eyeblinks)
    if (electrodeNum == 1 || electrodeNum == 2 || electrodeNum == 61)
        maxLimitElec = 300; minLimitElec = -300;
    else
        maxLimitElec = maxLimit; minLimitElec = minLimit;
    end
    
    analogData = loadAnalogData(fullfile(folderSegment,'LFP',['elec' num2str(electrodeNum) '.mat']));
    
    if showTrials
        hAllTrials = figure(11);
        subplot(8,8,i); plot(timeVals,analogData);title(['elec' num2str(i)]);
        axis('tight');
    end
    
    numTrials = size(analogData,1); %#ok<*NODEF>
    meanData = mean(analogData,2)';
    stdData  = std(analogData,[],2)';
    maxData  = max(analogData,[],2)';
    minData  = min(analogData,[],2)';
    
    clear tmpBadTrials1 tmpBadTrials2 tmpBadTrials3 tmpBadTrials4 
    tmpBadTrials1 = unique([find(maxData > meanData + threshold * stdData) find(minData < meanData - threshold * stdData)]);
    % Vinay - Ideally set maxLimit and minLimit for the below criteria to
    % be quite high/low so that only the extremely noisy trials are
    % rejected by these
    tmpBadTrials2 = unique(find(maxData > maxLimitElec));
    tmpBadTrials3 = unique(find(minData < minLimitElec));
    
    % Vinay - Set another criterion based on the deviation of each trial
    % from the mean trial signal. This way even if the actual signal shows
    % high deflections, if those deflections are consistent then the trials
    % are not rejected purely on the max/min criteria above
    meanTrialData = mean(analogData,1); % mean trial trace
    stdTrialData = std(analogData,[],1); % std across trials
%     maxTrialData = max(analogData,[],1);
%     minTrialData = min(analogData,[],1);
    tDplus = (meanTrialData + (threshold)*stdTrialData); % upper boundary/criterion
    tDminus = (meanTrialData - (threshold)*stdTrialData); % lower boundary/criterion
    
    % Check for trials that cross these boundaries and mark them as bad
    tBoolTrials = zeros(1,numTrials);
    for tr = 1:numTrials
        
        trialDeviationHigh = tDplus - analogData(tr,:); % deviation from the upper boundary
        trialDeviationLow = analogData(tr,:) - tDminus; % deviation from the lower boundary
        
        tBool = zeros(size(meanTrialData,1),size(meanTrialData,2));
        tBool(trialDeviationHigh<0) = 1; % set if the upper boundary is crossed anywhere
        tBool1 = sum(tBool); % number of times the upper boundary was crossed
        
        tBool = zeros(size(meanTrialData,1),size(meanTrialData,2));
        tBool(trialDeviationLow<0) = 1; % set if the lower boundary is crossed anywhere
        tBool2 = sum(tBool); % number of times the lower boundary was crossed

        tBoolTrials(tr) = tBool1 + tBool2; % number of times the boundaries were crossed
    end
    
    tmpBadTrials4 = find(tBoolTrials>0);
    
    allBadTrials{i} = unique([tmpBadTrials1 tmpBadTrials2 tmpBadTrials3 tmpBadTrials4]);
    
    goodTrials = 1:numTrials;
    goodTrials = setdiff(goodTrials,allBadTrials{i});
    
    clear meanData maxData minData stdData meanTrialData stdTrialData maxLimitElec minLimitElec
    clear trialDeviationHigh trialDeviationLow tBool tBool1 tBool2 tBoolTrials tDminus tDplus
    
    if showTrials && ~isempty(goodTrials)
        hGoodTrials = figure(12);
        subplot(8,8,i); plot(timeVals,analogData(goodTrials,:));title(['elec' num2str(i)]);
        axis('tight');
    end
end

% badTrials=allBadTrials{checkTheseElectrodes(1)};
badTrials=allBadTrials{1};
for i=1:numCheckElectrodes
    badTrials=intersect(badTrials,allBadTrials{i}); % in the previous case we took the union
end

disp(['total Trials: ' num2str(numTrials) ', bad trials: ' num2str(badTrials)]);

% saveData = 1;
for i=1:numElectrodes
    if length(allBadTrials{i}) ~= length(badTrials)
        disp(['Bad trials for electrode ' num2str(analogChannelsStored(i)) ': ' num2str(length(allBadTrials{i}))]);
    else
        disp(['Bad trials for electrode ' num2str(analogChannelsStored(i)) ': all (' num2str(length(badTrials)) ')']);
    end
end

if saveData
    disp(['Saving ' num2str(length(badTrials)) ' bad trials']);
    save(fullfile(folderSegment,'badTrials.mat'),'badTrials','allBadTrials','nameElec','checkTheseElectrodes','threshold','maxLimit','minLimit');
else
    disp('Bad trials will not be saved..'); %#ok<UNRCH>
end

end

function analogData = loadAnalogData(analogDataPath)
    load(analogDataPath);
end
