% Time-frequency analysis for different reference techniques
% Single, Bipolar, Average, CSD
% Vinay Shirhatti 16 Sep 2014
%
% modified from refTFAnalysisGeneric in GammaStimDicontinuityProject
% 27 March 2015
%==========================================================================


function refTFAnalysisCombined(subjectName,expDate,protocolName,folderSourceString,gridType,loadProtocolNumber)

if ~exist('folderSourceString','var')  folderSourceString='/media/store/';        end
if ~exist('gridType','var')            gridType='EEG';                  end
if ~exist('loadProtocolNumber','var')  loadProtocolNumber = 11;         end % Vinay - anything greater than 10 will read the protocolNumber from the saved data

%========================Define folder paths===============================
    folderName = fullfile(folderSourceString,'data',subjectName,gridType,expDate,protocolName);

    % Get folders
    folderExtract = fullfile(folderName,'extractedData');
    folderSegment = fullfile(folderName,'segmentedData');
    folderLFP = fullfile(folderSegment,'LFP');
    folderSpikes = fullfile(folderSegment,'Spikes');
    
    folderBipolar = fullfile(folderLFP,'bipolar');
    folderAverage = folderLFP;
    folderCSD = fullfile(folderLFP,'csd');
    
%-----

% load LFP Information
[analogChannelsStored,timeVals,~,analogInputNums] = loadlfpInfo(folderLFP);
[neuralChannelsStored,SourceUnitIDs] = loadspikeInfo(folderSpikes);

% Get Combinations
[~,aValsUnique,eValsUnique,sValsUnique,fValsUnique,oValsUnique,cValsUnique,tValsUnique,rValsUnique,pValsUnique] = loadParameterCombinations(folderExtract);

% Get properties of the Stimulus
% stimResults = loadStimResults(folderExtract);

% Vinay - read the variable loadProtocolNumber as protocolNumber if it is
% less than 11. Otherwise it will be loaded from stimResults. This is so
% that one can directly pass the protocolNumber in case it is not recorded
% in stimResults
if strncmp(protocolName,'CRS',3)
    if (loadProtocolNumber < 11)
        protocolNumber = loadProtocolNumber;
    else
        protocolNumber = getProtocolNumber(folderExtract);
    end
else
    protocolNumber = 11;
end

% Vinay - get the particular gabors that were displayed in the protocol.
gaborsDisplayed = getGaborsDisplayed(folderExtract);
disp(['Gabors displayed in this protocol:' num2str(gaborsDisplayed)]);

% Vinay - if all the gabors were not displayed then assign some default
% value to the null gabor's parameters

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Display main options
% fonts
fontSizeSmall = 10; fontSizeMedium = 12; fontSizeLarge = 16; fontSizeTiny = 8;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make Panels
panelHeight = 0.3; panelStartHeight = 0.68;
staticPanelWidth = 0.22; staticStartPos = 0.025; % [Vinay] - changed width 0.25 to 0.22
dynamicPanelWidth = 0.23; dynamicStartPos = 0.245; % [Vinay] - changed width 0.25 to 0.23
timingPanelWidth = 0.20; timingStartPos = 0.475; % [Vinay] - changed width 0.25 to 0.20
plotOptionsPanelWidth = 0.15; plotOptionsStartPos = 0.675; % [Vinay] - changed width 0.2 to 0.15
selectParamPanelWidth = 0.15; selectParamStartPos = 0.825; % [Vinay] - adding a panel to select parameters
backgroundColor = 'w';

% Vinay - define a flag for notching the line noise
notchData = 0;
% [Vinay] - define a flag to take bipolar (i.e. diff) signals: V1-V2
useBipolar = 0; % (time domain) electrode1 - electrode2
useDiff = 0; % (freq domain) spectrum1 -spectrum2
% [Vinay]- define a flag to plot SEM for the ERP plots
plotSEM = 1;
% [Vinay] - define a flag to save MP data if it is set
saveMPFlag = 1;
% [Vinay] - define a flag to save HHT data if it is set
saveHHTFlag = 1;
% [Vinay] - define a flag to use badTrials for each electrode separately
useAllBadTrials = 1;
% [Vinay] - define a flag to use intersection of good trials across
% electrodes or else to use individual electrodewise good trials
intersectTrials = 0;

% [Vinay] - set this flag if bipolar data is stored
existsBipolarData = 0;
% [Vinay] - set this flag if log of every trials is to be taken before
% averaging (for TF plots)
takeLogTrial = 0;
% [Vinay] - set this if the bipolar pairs are made using complimentary
% contralateral electrodes and not the nearest neighbour
contralateralNeighbour = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Dynamic panel %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dynamicHeight = 0.06; dynamicGap=0.015; dynamicTextWidth = 0.6;
hDynamicPanel = uipanel('Title','Parameters','fontSize', fontSizeLarge, ...
    'Unit','Normalized','Position',[dynamicStartPos panelStartHeight dynamicPanelWidth panelHeight]);

% Analog channel
[analogChannelStringList,analogChannelStringArray] = getAnalogStringFromValues(analogChannelsStored,analogInputNums);
uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'Position',[0 1-(dynamicHeight+dynamicGap) dynamicTextWidth dynamicHeight],...
    'Style','text','String','Analog Channel (1:2)','FontSize',fontSizeTiny);
hAnalogChannel = uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [dynamicTextWidth 1-(dynamicHeight+dynamicGap) (0.5*(1-dynamicTextWidth)) dynamicHeight], ...
    'Style','popup','String',analogChannelStringList,'FontSize',fontSizeTiny);
% [Vinay] - adding analog channel 2 for chosing bipolar signals V1-V2 
hAnalogChannel2 = uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [(dynamicTextWidth+(0.5*(1-dynamicTextWidth))) 1-(dynamicHeight+dynamicGap) (0.5*(1-dynamicTextWidth)) dynamicHeight], ...
    'Style','popup','String',(['none|' analogChannelStringList]),'FontSize',fontSizeTiny,'Callback',{@bipolar_Callback});

% Neural channel
uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'Position',[0 1-2*(dynamicHeight+dynamicGap) dynamicTextWidth dynamicHeight],...
    'Style','text','String','Neural Channel','FontSize',fontSizeTiny);
    
if ~isempty(neuralChannelsStored)
    neuralChannelString = getNeuralStringFromValues(neuralChannelsStored,SourceUnitIDs);
    hNeuralChannel = uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
        'BackgroundColor', backgroundColor, 'Position', ...
        [dynamicTextWidth 1-2*(dynamicHeight+dynamicGap) 1-dynamicTextWidth dynamicHeight], ...
        'Style','popup','String',neuralChannelString,'FontSize',fontSizeTiny);
else
    hNeuralChannel = uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
        'Position', [dynamicTextWidth 1-2*(dynamicHeight+dynamicGap) 1-dynamicTextWidth dynamicHeight], ...
        'Style','text','String','Not found','FontSize',fontSizeTiny);
end


if strncmp(protocolName,'GRF',3)
    % Sigma
    sigmaString = getStringFromValuesGRF(sValsUnique,1);
    uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
        'Position',[0 1-3*(dynamicHeight+dynamicGap) dynamicTextWidth dynamicHeight], ...
        'Style','text','String','Sigma (Deg)','FontSize',fontSizeTiny);
    hSigma = uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
        'BackgroundColor', backgroundColor, 'Position', ...
        [dynamicTextWidth 1-3*(dynamicHeight+dynamicGap) 1-dynamicTextWidth dynamicHeight], ...
        'Style','popup','String',sigmaString,'FontSize',fontSizeTiny);

    % Spatial Frequency
    spatialFreqString = getStringFromValuesGRF(fValsUnique,1);
    uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
        'Position',[0 1-4*(dynamicHeight+dynamicGap) dynamicTextWidth dynamicHeight], ...
        'Style','text','String','Spatial Freq (CPD)','FontSize',fontSizeTiny);
    hSpatialFreq = uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
        'BackgroundColor', backgroundColor, 'Position', ...
        [dynamicTextWidth 1-4*(dynamicHeight+dynamicGap) 1-dynamicTextWidth dynamicHeight], ...
        'Style','popup','String',spatialFreqString,'FontSize',fontSizeTiny);

    % Orientation
    orientationString = getStringFromValuesGRF(oValsUnique,1);
    uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
        'Position',[0 1-5*(dynamicHeight+dynamicGap) dynamicTextWidth dynamicHeight], ...
        'Style','text','String','Orientation (Deg)','FontSize',fontSizeTiny);
    hOrientation = uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
        'BackgroundColor', backgroundColor, 'Position', ...
        [dynamicTextWidth 1-5*(dynamicHeight+dynamicGap) 1-dynamicTextWidth dynamicHeight], ...
        'Style','popup','String',orientationString,'FontSize',fontSizeTiny);

    % Contrast
    contrastString = getStringFromValuesGRF(cValsUnique,1);
    uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
        'Position',[0 1-6*(dynamicHeight+dynamicGap) dynamicTextWidth dynamicHeight], ...
        'Style','text','String','Contrast (%)','FontSize',fontSizeTiny);
    hContrast = uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
        'BackgroundColor', backgroundColor, 'Position', ...
        [dynamicTextWidth 1-6*(dynamicHeight+dynamicGap) 1-dynamicTextWidth dynamicHeight], ...
        'Style','popup','String',contrastString,'FontSize',fontSizeTiny);

    % Temporal Frequency
    temporalFreqString = getStringFromValuesGRF(tValsUnique,1);
    uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
        'Position',[0 1-7*(dynamicHeight+dynamicGap) dynamicTextWidth dynamicHeight], ...
        'Style','text','String','Temporal Freq (Hz)','FontSize',fontSizeTiny);
    hTemporalFreq = uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
        'BackgroundColor', backgroundColor, 'Position', ...
        [dynamicTextWidth 1-7*(dynamicHeight+dynamicGap) 1-dynamicTextWidth dynamicHeight], ...
        'Style','popup','String',temporalFreqString,'FontSize',fontSizeTiny);

    % Azimuth
    azimuthString = getStringFromValuesGRF(aValsUnique,1);
    uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
        'Position',[0 1-8*(dynamicHeight+dynamicGap) dynamicTextWidth dynamicHeight],...
        'Style','text','String','Azimuth (Deg)','FontSize',fontSizeTiny);
    hAzimuth = uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
        'BackgroundColor', backgroundColor, 'Position', ...
        [dynamicTextWidth 1-8*(dynamicHeight+dynamicGap) 1-dynamicTextWidth dynamicHeight], ...
        'Style','popup','String',azimuthString,'FontSize',fontSizeTiny);

    % Elevation
    elevationString = getStringFromValuesGRF(eValsUnique,1);
    uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
        'Position',[0 1-9*(dynamicHeight+dynamicGap) dynamicTextWidth dynamicHeight], ...
        'Style','text','String','Elevation (Deg)','FontSize',fontSizeTiny);
    hElevation = uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
        'BackgroundColor', backgroundColor, 'Position',...
        [dynamicTextWidth 1-9*(dynamicHeight+dynamicGap) 1-dynamicTextWidth dynamicHeight], ...
        'Style','popup','String',elevationString,'FontSize',fontSizeTiny);
    
    %======================================================================
    % Select the two parameters whose values are to be varied across
    % plots/dimensions
    %======================================================================
    parametersString = 'NA|azimuth|elevation|sigma|spatialFreq|orientation|contrast|temporalFreq';

    uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
        'Position',[0 1-10*(dynamicHeight+dynamicGap) dynamicTextWidth dynamicHeight], ...
        'Style','text','String','vary1','FontSize',fontSizeMedium);
    hParam7 = uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
        'BackgroundColor', backgroundColor, 'Position', ...
        [dynamicTextWidth 1-10*(dynamicHeight+dynamicGap) 1-dynamicTextWidth dynamicHeight], ...
        'Style','popup','String',parametersString,'FontSize',fontSizeTiny,'Callback',{@resetVaryParams_Callback});

    uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
        'Position',[0 1-12*(dynamicHeight+dynamicGap) dynamicTextWidth dynamicHeight], ...
        'Style','text','String','vary2','FontSize',fontSizeMedium);
    hParam8 = uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
        'BackgroundColor', backgroundColor, 'Position', ...
        [dynamicTextWidth 1-12*(dynamicHeight+dynamicGap) 1-dynamicTextWidth dynamicHeight], ...
        'Style','popup','String',parametersString,'FontSize',fontSizeTiny,'Callback',{@resetVaryParams_Callback});


elseif strncmp(protocolName,'CRS',3)
    %[Vinay] - create 3 columns corresponding to the parameters of the 3 gabors
    %-Centre, Ring, Surround. Textwidth = 0.4, C,R,S popout - each 0.2 width

    labelWidth = 0.4; crsWidth = 0.2; % [Vinay] - labelWidth + 3(crsWidth) = 1

    uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
        'Position',[labelWidth 1-3*(dynamicHeight+dynamicGap) (3*crsWidth) dynamicHeight],...
        'Style','text','String','Surround        Ring        Centre','FontSize',fontSizeTiny);


    for gaborNum = 1:3

        % Azimuth
        azimuthString = getStringFromValues(aValsUnique,1, gaborNum);
        uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
        'Position',[0 1-4*(dynamicHeight+dynamicGap) labelWidth dynamicHeight],...
        'Style','text','String','Azimuth (Deg)','FontSize',fontSizeTiny);
        hAzimuth(gaborNum) = uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
        'BackgroundColor', backgroundColor, 'Position', ...
        [(labelWidth+(gaborNum-1)*crsWidth) 1-4*(dynamicHeight+dynamicGap) crsWidth dynamicHeight], ...
        'Style','popup','String',azimuthString,'FontSize',fontSizeTiny);

        % Elevation
        elevationString = getStringFromValues(eValsUnique,1, gaborNum);
        uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
        'Position',[0 1-5*(dynamicHeight+dynamicGap) labelWidth dynamicHeight], ...
        'Style','text','String','Elevation (Deg)','FontSize',fontSizeTiny);
        hElevation(gaborNum) = uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
        'BackgroundColor', backgroundColor, 'Position',...
        [(labelWidth+(gaborNum-1)*crsWidth) 1-5*(dynamicHeight+dynamicGap) crsWidth dynamicHeight], ...
        'Style','popup','String',elevationString,'FontSize',fontSizeTiny);

        % Sigma
        sigmaString = getStringFromValues(sValsUnique,1, gaborNum);
        uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
        'Position',[0 1-6*(dynamicHeight+dynamicGap) labelWidth dynamicHeight], ...
        'Style','text','String','Sigma (Deg)','FontSize',fontSizeTiny);
        hSigma(gaborNum) = uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
        'BackgroundColor', backgroundColor, 'Position', ...
        [(labelWidth+(gaborNum-1)*crsWidth) 1-6*(dynamicHeight+dynamicGap) crsWidth dynamicHeight], ...
        'Style','popup','String',sigmaString,'FontSize',fontSizeTiny);

        % Spatial Frequency
        spatialFreqString = getStringFromValues(fValsUnique,1, gaborNum);
        uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
        'Position',[0 1-7*(dynamicHeight+dynamicGap) labelWidth dynamicHeight], ...
        'Style','text','String','Spatial Freq (CPD)','FontSize',fontSizeTiny);
        hSpatialFreq(gaborNum) = uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
        'BackgroundColor', backgroundColor, 'Position', ...
        [(labelWidth+(gaborNum-1)*crsWidth) 1-7*(dynamicHeight+dynamicGap) crsWidth dynamicHeight], ...
        'Style','popup','String',spatialFreqString,'FontSize',fontSizeTiny);

        % Orientation
        orientationString = getStringFromValues(oValsUnique,1, gaborNum);
        uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
        'Position',[0 1-8*(dynamicHeight+dynamicGap) labelWidth dynamicHeight], ...
        'Style','text','String','Orientation (Deg)','FontSize',fontSizeTiny);
        hOrientation(gaborNum) = uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
        'BackgroundColor', backgroundColor, 'Position', ...
        [(labelWidth+(gaborNum-1)*crsWidth) 1-8*(dynamicHeight+dynamicGap) crsWidth dynamicHeight], ...
        'Style','popup','String',orientationString,'FontSize',fontSizeTiny);

        % Contrast
        contrastString = getStringFromValues(cValsUnique,1, gaborNum);
        uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
        'Position',[0 1-9*(dynamicHeight+dynamicGap) labelWidth dynamicHeight], ...
        'Style','text','String','Contrast (%)','FontSize',fontSizeTiny);
        hContrast(gaborNum) = uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
        'BackgroundColor', backgroundColor, 'Position', ...
        [(labelWidth+(gaborNum-1)*crsWidth) 1-9*(dynamicHeight+dynamicGap) crsWidth dynamicHeight], ...
        'Style','popup','String',contrastString,'FontSize',fontSizeTiny);

        % Temporal Frequency
        temporalFreqString = getStringFromValues(tValsUnique,1, gaborNum);
        uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
        'Position',[0 1-10*(dynamicHeight+dynamicGap) labelWidth dynamicHeight], ...
        'Style','text','String','Temporal Freq (Hz)','FontSize',fontSizeTiny);
        hTemporalFreq(gaborNum) = uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
        'BackgroundColor', backgroundColor, 'Position', ...
        [(labelWidth+(gaborNum-1)*crsWidth) 1-10*(dynamicHeight+dynamicGap) crsWidth dynamicHeight], ...
        'Style','popup','String',temporalFreqString,'FontSize',fontSizeTiny);

        % [Vinay] - adding radius and spatial phase cells
        % radius
        radiusString = getStringFromValues(rValsUnique,1, gaborNum);
        uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
        'Position',[0 1-11*(dynamicHeight+dynamicGap) labelWidth dynamicHeight], ...
        'Style','text','String','Radius (deg)','FontSize',fontSizeTiny);
        hRadius(gaborNum) = uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
        'BackgroundColor', backgroundColor, 'Position', ...
        [(labelWidth+(gaborNum-1)*crsWidth) 1-11*(dynamicHeight+dynamicGap) crsWidth dynamicHeight], ...
        'Style','popup','String',radiusString,'FontSize',fontSizeTiny);

        % spatial phase
        spatialPhaseString = getStringFromValues(pValsUnique,1, gaborNum);
        uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
        'Position',[0 1-12*(dynamicHeight+dynamicGap) labelWidth dynamicHeight], ...
        'Style','text','String','Spatial Phase (Deg)','FontSize',fontSizeTiny);
        hSpatialPhase(gaborNum) = uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
        'BackgroundColor', backgroundColor, 'Position', ...
        [(labelWidth+(gaborNum-1)*crsWidth) 1-12*(dynamicHeight+dynamicGap) crsWidth dynamicHeight], ...
        'Style','popup','String',spatialPhaseString,'FontSize',fontSizeTiny);

    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Timing panel %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
timingHeight = 0.1; timingTextWidth = 0.5; timingBoxWidth = 0.20; % [Vinay] - changed width from 0.25 to 0.20 
hTimingPanel = uipanel('Title','Timing','fontSize', fontSizeLarge, ...
    'Unit','Normalized','Position',[timingStartPos panelStartHeight timingPanelWidth panelHeight]);

% Vinay - changed these values for EEG
signalRange = [-0.6 1.0];
fftRange = [0 120];
baseline = [-0.5 0];
stimPeriod = [0.25 0.75];

% Signal Range
uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'Position',[0 1-timingHeight timingTextWidth timingHeight], ...
    'Style','text','String','Parameter','FontSize',fontSizeMedium);

uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'Position',[timingTextWidth 1-timingHeight timingBoxWidth timingHeight], ...
    'Style','text','String','Min','FontSize',fontSizeMedium);

uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'Position',[timingTextWidth+timingBoxWidth 1-timingHeight timingBoxWidth timingHeight], ...
    'Style','text','String','Max','FontSize',fontSizeMedium);

% Stim Range
uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'Position',[0 1-3*timingHeight timingTextWidth timingHeight], ...
    'Style','text','String','Stim Range (s)','FontSize',fontSizeSmall);
hStimMin = uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[timingTextWidth 1-3*timingHeight timingBoxWidth timingHeight], ...
    'Style','edit','String',num2str(signalRange(1)),'FontSize',fontSizeSmall);
hStimMax = uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[timingTextWidth+timingBoxWidth 1-3*timingHeight timingBoxWidth timingHeight], ...
    'Style','edit','String',num2str(signalRange(2)),'FontSize',fontSizeSmall);

% FFT Range
uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'Position',[0 1-5*timingHeight timingTextWidth timingHeight], ...
    'Style','text','String','FFT Range (Hz)','FontSize',fontSizeSmall);
hFFTMin = uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[timingTextWidth 1-5*timingHeight timingBoxWidth timingHeight], ...
    'Style','edit','String',num2str(fftRange(1)),'FontSize',fontSizeSmall);
hFFTMax = uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[timingTextWidth+timingBoxWidth 1-5*timingHeight timingBoxWidth timingHeight], ...
    'Style','edit','String',num2str(fftRange(2)),'FontSize',fontSizeSmall);

% Baseline
uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'Position',[0 1-6*timingHeight timingTextWidth timingHeight], ...
    'Style','text','String','Basline (s)','FontSize',fontSizeSmall);
hBaselineMin = uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[timingTextWidth 1-6*timingHeight timingBoxWidth timingHeight], ...
    'Style','edit','String',num2str(baseline(1)),'FontSize',fontSizeSmall);
hBaselineMax = uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[timingTextWidth+timingBoxWidth 1-6*timingHeight timingBoxWidth timingHeight], ...
    'Style','edit','String',num2str(baseline(2)),'FontSize',fontSizeSmall);

% Stim Period
uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'Position',[0 1-7*timingHeight timingTextWidth timingHeight], ...
    'Style','text','String','Stim period (s)','FontSize',fontSizeSmall);
hStimPeriodMin = uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[timingTextWidth 1-7*timingHeight timingBoxWidth timingHeight], ...
    'Style','edit','String',num2str(stimPeriod(1)),'FontSize',fontSizeSmall);
hStimPeriodMax = uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[timingTextWidth+timingBoxWidth 1-7*timingHeight timingBoxWidth timingHeight], ...
    'Style','edit','String',num2str(stimPeriod(2)),'FontSize',fontSizeSmall);

% Y Range
uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'Position',[0 1-8*timingHeight timingTextWidth timingHeight], ...
    'Style','text','String','Y Range','FontSize',fontSizeSmall);
hYMin = uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[timingTextWidth 1-8*timingHeight timingBoxWidth timingHeight], ...
    'Style','edit','String','0','FontSize',fontSizeSmall);
hYMax = uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[timingTextWidth+timingBoxWidth 1-8*timingHeight timingBoxWidth timingHeight], ...
    'Style','edit','String','1','FontSize',fontSizeSmall);

% STA length
fBandLow = 30; fBandHigh = 80;
uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'Position',[0 1-9*timingHeight timingTextWidth timingHeight], ...
    'Style','text','String','freqBand range','FontSize',fontSizeSmall);
hfBandLow = uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[timingTextWidth 1-9*timingHeight timingBoxWidth timingHeight], ...
    'Style','edit','String',num2str(fBandLow),'FontSize',fontSizeSmall);
hfBandHigh = uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[timingTextWidth+timingBoxWidth 1-9*timingHeight timingBoxWidth timingHeight], ...
    'Style','edit','String',num2str(fBandHigh),'FontSize',fontSizeSmall);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plot Options %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plotOptionsHeight = 0.1;
hPlotOptionsPanel = uipanel('Title','Plotting Options','fontSize', fontSizeLarge, ...
    'Unit','Normalized','Position',[plotOptionsStartPos panelStartHeight plotOptionsPanelWidth panelHeight]);

% Chose color
[colorString, colorNames] = getColorString;
uicontrol('Parent',hPlotOptionsPanel,'Unit','Normalized', ...
    'Position',[0 1-plotOptionsHeight 0.6 plotOptionsHeight], ...
    'Style','text','String','Color','FontSize',fontSizeTiny);

hChooseColor = uicontrol('Parent',hPlotOptionsPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[0.6 1-plotOptionsHeight 0.4 plotOptionsHeight], ...
    'Style','popup','String',colorString,'FontSize',fontSizeTiny);

% Vinay - adding an option to adjust lineWidth of plots
linewidths = [1,1.5,2,2.5,3];
linewidthString = '1|1.5|2|2.5|3|';
uicontrol('Parent',hPlotOptionsPanel,'Unit','Normalized', ...
    'Position',[0 1-2*plotOptionsHeight 0.6 plotOptionsHeight], ...
    'Style','text','String','lineWidth','FontSize',fontSizeTiny);

hChooseWidth = uicontrol('Parent',hPlotOptionsPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[0.6 1-2*plotOptionsHeight 0.4 plotOptionsHeight], ...
    'Style','popup','String',linewidthString,'FontSize',fontSizeTiny);


% Analysis Type
analysisTypeString = 'ERP|FFT|delta FFT|Band Power';
uicontrol('Parent',hPlotOptionsPanel,'Unit','Normalized', ...
    'Position',[0 1-3*plotOptionsHeight 0.6 plotOptionsHeight], ...
    'Style','text','String','Analysis Type','FontSize',fontSizeSmall);
hAnalysisType = uicontrol('Parent',hPlotOptionsPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [0.6 1-3*plotOptionsHeight 0.4 plotOptionsHeight], ...
    'Style','popup','String',analysisTypeString,'FontSize',fontSizeSmall);


% Vinay - adding a toggle button to use notched data
hNotchData = uicontrol('Parent',hPlotOptionsPanel,'Unit','Normalized', ...
    'Position',[0 5*plotOptionsHeight 1 plotOptionsHeight], ...
    'Style','togglebutton','String','raw/lineNoiseRemoved','FontSize',fontSizeMedium, ...
    'Callback',{@notch_Callback});

uicontrol('Parent',hPlotOptionsPanel,'Unit','Normalized', ...
    'Position',[0 4*plotOptionsHeight 1 plotOptionsHeight], ...
    'Style','pushbutton','String','cla','FontSize',fontSizeMedium, ...
    'Callback',{@cla_Callback});

hHoldOn = uicontrol('Parent',hPlotOptionsPanel,'Unit','Normalized', ...
    'Position',[0 3*plotOptionsHeight 1 plotOptionsHeight], ...
    'Style','togglebutton','String','hold on','FontSize',fontSizeMedium, ...
    'Callback',{@holdOn_Callback});

uicontrol('Parent',hPlotOptionsPanel,'Unit','Normalized', ...
    'Position',[0 2*plotOptionsHeight 1 plotOptionsHeight], ...
    'Style','pushbutton','String','rescale Y','FontSize',fontSizeMedium, ...
    'Callback',{@rescaleY_Callback});

uicontrol('Parent',hPlotOptionsPanel,'Unit','Normalized', ...
    'Position',[0 plotOptionsHeight 1 plotOptionsHeight], ...
    'Style','pushbutton','String','rescale X','FontSize',fontSizeMedium, ...
    'Callback',{@rescaleData_Callback});

uicontrol('Parent',hPlotOptionsPanel,'Unit','Normalized', ...
    'Position',[0 0 1 plotOptionsHeight], ...
    'Style','pushbutton','String','plot','FontSize',fontSizeMedium, ...
    'Callback',{@plotData_Callback});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [Vinay] - adding a panel to select the parameters for the plots
% only for CRS

if strncmp(protocolName,'CRS',3)

    selectParamHeight = 0.1; selectParamWidth = 0.5; selectParamBoxWidth = 0.20; % [Vinay] - changed width from 0.25 to 0.20
    textWidth = 0.4; selGWidth = 0.2; selPWidth = 0.4;
    hSelectParamPanel = uipanel('Title','Select Parameters','fontSize', fontSizeLarge, ...
        'Unit','Normalized','Position',[selectParamStartPos panelStartHeight selectParamPanelWidth panelHeight]);


    parametersString = 'NA|azimuth|elevation|sigma|spatialFreq|orientation|contrast|temporalFreq|radius|spatialPhase';
    gaborString = 'NA|C|R|S';

    %parameter1
    uicontrol('Parent',hSelectParamPanel,'Unit','Normalized', ...
        'Position',[0 1-selectParamHeight textWidth selectParamHeight], ...
        'Style','text','String','param1','FontSize',fontSizeSmall);
    hGabor1 = uicontrol('Parent',hSelectParamPanel,'Unit','Normalized', ...
        'BackgroundColor', backgroundColor, 'Position', ...
        [textWidth 1-selectParamHeight selGWidth selectParamHeight], ...
        'Style','popup','String',gaborString,'FontSize',fontSizeTiny);
    hParam1 = uicontrol('Parent',hSelectParamPanel,'Unit','Normalized', ...
        'BackgroundColor', backgroundColor, 'Position', ...
        [(textWidth+selGWidth) 1-selectParamHeight selPWidth selectParamHeight], ...
        'Style','popup','String',parametersString,'FontSize',fontSizeTiny);


    %parameter2
    uicontrol('Parent',hSelectParamPanel,'Unit','Normalized', ...
        'Position',[0 1-2*selectParamHeight textWidth selectParamHeight], ...
        'Style','text','String','param2','FontSize',fontSizeSmall);
    hGabor2 = uicontrol('Parent',hSelectParamPanel,'Unit','Normalized', ...
        'BackgroundColor', backgroundColor, 'Position', ...
        [textWidth 1-2*selectParamHeight selGWidth selectParamHeight], ...
        'Style','popup','String',gaborString,'FontSize',fontSizeTiny);
    hParam2 = uicontrol('Parent',hSelectParamPanel,'Unit','Normalized', ...
        'BackgroundColor', backgroundColor, 'Position', ...
        [(textWidth+selGWidth) 1-2*selectParamHeight selPWidth selectParamHeight], ...
        'Style','popup','String',parametersString,'FontSize',fontSizeTiny);


    %parameter3
    uicontrol('Parent',hSelectParamPanel,'Unit','Normalized', ...
        'Position',[0 1-3*selectParamHeight textWidth selectParamHeight], ...
        'Style','text','String','param3','FontSize',fontSizeSmall);
    hGabor3 = uicontrol('Parent',hSelectParamPanel,'Unit','Normalized', ...
        'BackgroundColor', backgroundColor, 'Position', ...
        [textWidth 1-3*selectParamHeight selGWidth selectParamHeight], ...
        'Style','popup','String',gaborString,'FontSize',fontSizeTiny);
    hParam3 = uicontrol('Parent',hSelectParamPanel,'Unit','Normalized', ...
        'BackgroundColor', backgroundColor, 'Position', ...
        [(textWidth+selGWidth) 1-3*selectParamHeight selPWidth selectParamHeight], ...
        'Style','popup','String',parametersString,'FontSize',fontSizeTiny);


    %parameter4
    uicontrol('Parent',hSelectParamPanel,'Unit','Normalized', ...
        'Position',[0 1-4*selectParamHeight textWidth selectParamHeight], ...
        'Style','text','String','param4','FontSize',fontSizeSmall);
    hGabor4 = uicontrol('Parent',hSelectParamPanel,'Unit','Normalized', ...
        'BackgroundColor', backgroundColor, 'Position', ...
        [textWidth 1-4*selectParamHeight selGWidth selectParamHeight], ...
        'Style','popup','String',gaborString,'FontSize',fontSizeTiny);
    hParam4 = uicontrol('Parent',hSelectParamPanel,'Unit','Normalized', ...
        'BackgroundColor', backgroundColor, 'Position', ...
        [(textWidth+selGWidth) 1-4*selectParamHeight selPWidth selectParamHeight], ...
        'Style','popup','String',parametersString,'FontSize',fontSizeTiny);


    %parameter5
    uicontrol('Parent',hSelectParamPanel,'Unit','Normalized', ...
        'Position',[0 1-5*selectParamHeight textWidth selectParamHeight], ...
        'Style','text','String','param5','FontSize',fontSizeSmall);
    hGabor5 = uicontrol('Parent',hSelectParamPanel,'Unit','Normalized', ...
        'BackgroundColor', backgroundColor, 'Position', ...
        [textWidth 1-5*selectParamHeight selGWidth selectParamHeight], ...
        'Style','popup','String',gaborString,'FontSize',fontSizeTiny);
    hParam5 = uicontrol('Parent',hSelectParamPanel,'Unit','Normalized', ...
        'BackgroundColor', backgroundColor, 'Position', ...
        [(textWidth+selGWidth) 1-5*selectParamHeight selPWidth selectParamHeight], ...
        'Style','popup','String',parametersString,'FontSize',fontSizeTiny);


    %parameter6
    uicontrol('Parent',hSelectParamPanel,'Unit','Normalized', ...
        'Position',[0 1-6*selectParamHeight textWidth selectParamHeight], ...
        'Style','text','String','param6','FontSize',fontSizeSmall);
    hGabor6 = uicontrol('Parent',hSelectParamPanel,'Unit','Normalized', ...
        'BackgroundColor', backgroundColor, 'Position', ...
        [textWidth 1-6*selectParamHeight selGWidth selectParamHeight], ...
        'Style','popup','String',gaborString,'FontSize',fontSizeTiny);
    hParam6 = uicontrol('Parent',hSelectParamPanel,'Unit','Normalized', ...
        'BackgroundColor', backgroundColor, 'Position', ...
        [(textWidth+selGWidth) 1-6*selectParamHeight selPWidth selectParamHeight], ...
        'Style','popup','String',parametersString,'FontSize',fontSizeTiny);



    uicontrol('Parent',hSelectParamPanel,'Unit','Normalized', ...
        'Position',[0 1-7*selectParamHeight textWidth selectParamHeight], ...
        'Style','text','String','vary1','FontSize',fontSizeMedium);
    hGabor7 = uicontrol('Parent',hSelectParamPanel,'Unit','Normalized', ...
        'BackgroundColor', backgroundColor, 'Position', ...
        [textWidth 1-7*selectParamHeight selGWidth selectParamHeight], ...
        'Style','popup','String',gaborString,'FontSize',fontSizeTiny,'Callback',{@resetVaryParams_Callback});
    hParam7 = uicontrol('Parent',hSelectParamPanel,'Unit','Normalized', ...
        'BackgroundColor', backgroundColor, 'Position', ...
        [(textWidth+selGWidth) 1-7*selectParamHeight selPWidth selectParamHeight], ...
        'Style','popup','String',parametersString,'FontSize',fontSizeTiny,'Callback',{@resetVaryParams_Callback});


    uicontrol('Parent',hSelectParamPanel,'Unit','Normalized', ...
        'Position',[0 1-8*selectParamHeight textWidth selectParamHeight], ...
        'Style','text','String','vary2','FontSize',fontSizeMedium);
    hGabor8 = uicontrol('Parent',hSelectParamPanel,'Unit','Normalized', ...
        'BackgroundColor', backgroundColor, 'Position', ...
        [textWidth 1-8*selectParamHeight selGWidth selectParamHeight], ...
        'Style','popup','String',gaborString,'FontSize',fontSizeTiny,'Callback',{@resetVaryParams_Callback});
    hParam8 = uicontrol('Parent',hSelectParamPanel,'Unit','Normalized', ...
        'BackgroundColor', backgroundColor, 'Position', ...
        [(textWidth+selGWidth) 1-8*selectParamHeight selPWidth selectParamHeight], ...
        'Style','popup','String',parametersString,'FontSize',fontSizeTiny,'Callback',{@resetVaryParams_Callback});



    % check box - 'show default'. This will load default parameters as per the
    % specific protocol being used
    useDefParams = 0; % by default don't use the default parameters
    hUseDefParams = uicontrol('Parent',hSelectParamPanel,'Unit','Normalized', ...
        'BackgroundColor', backgroundColor, 'Position', ...
        [0.1 1-9*selectParamHeight 0.8 selectParamHeight], ...
        'Style','checkbox','String','Show default parameters','FontSize',fontSizeTiny, 'Callback',{@useDefParams_Callback});

    % display the specific protocol being used
    protocolNameString = {'NoneProtocol';'RingProtocol';'ContrastRingProtocol';'DualContrastProtocol'; ...
        'DualOrientationProtocol';'DualPhaseProtocol';'OrientationRingProtocol';'PhaseRingProtocol'; ...
        'Drifting Ring Protocol';'CrossOrientationProtocol';'AnnulusRingProtocol'};
    % protocolNumber = getProtocolNumber(folderExtract);

    uicontrol('Parent',hSelectParamPanel,'Unit','Normalized', ...
        'Position',[0.1 1-10*selectParamHeight 0.8 selectParamHeight], ...
        'Style','text','String',protocolNameString(protocolNumber+1),'FontSize',fontSizeTiny,'FontWeight','bold');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get plots and message handles

% Get electrode array information
gridLayout = 1; % 64 electrodes
electrodeGridPos = [staticStartPos panelStartHeight staticPanelWidth panelHeight];
hElectrodes = showElectrodeLocations(electrodeGridPos,analogChannelsStored(get(hAnalogChannel,'val')), ...
    colorNames(get(hChooseColor,'val')),[],1,0,gridType,subjectName,gridLayout);


% Set up the plot box and its dimensions
startXPos = staticStartPos; endXPos = 0.95; startYPos = 0.05; endYPos = 0.60;
centerGap = 0.05;
plotsWidthX = (endXPos-startXPos-centerGap);
plotsWidthY = (endYPos-startYPos-centerGap);
gap = 0.02; gapSmall = 0.002;
varyGridPos=[startXPos startYPos endXPos endYPos];


uicontrol('Unit','Normalized','Position',[0 0.975 1 0.025],...
    'Style','text','String',[subjectName expDate protocolName],'FontSize',fontSizeLarge);
        

if strncmp(protocolName,'CRS',3)

    % [Vinay] - load default parameters
    set(hGabor1,'val',4); % 4 => S
    set(hGabor2,'val',4);
    set(hGabor3,'val',3); % 3 => R
    set(hGabor4,'val',3);
    set(hGabor5,'val',2); % 2 => C
    set(hGabor6,'val',2);
    set(hGabor7,'val',1); % 1 => NA
    set(hGabor8,'val',1);

    set(hParam1,'val',7); % 7 => cont (surr)
    set(hParam2,'val',6); % 6 => ori
    set(hParam3,'val',7); % 7 => cont (ring)
    set(hParam4,'val',9); % 9 => rad 
    set(hParam5,'val',7); % 7 => cont (centre)
    set(hParam6,'val',6); % 6 => ori
    set(hParam7,'val',1); % 1 => NA
    set(hParam8,'val',1);

    adjFactor = 5;

    gabor7 = adjFactor - get(hGabor7,'val');
    gabor8 = adjFactor - get(hGabor8,'val');

    param7 = get(hParam7, 'val');
    param8 = get(hParam8, 'val');


    hVaryParamPlot = getPlotHandles(length(getValsUnique(gabor7, param7)),length(getValsUnique(gabor8, param8)),varyGridPos,gapSmall);

elseif strncmp(protocolName,'GRF',3)
    
    param7 = get(hParam7, 'val');
    param8 = get(hParam8, 'val');
    
    gabor7 = []; gabor8 = [];

    hVaryParamPlot = getPlotHandles(length(getValsUniqueGRF(param7)),length(getValsUniqueGRF(param8)),varyGridPos,gapSmall);

end

%**************************************************************************
% [Vinay] - TIME-FREQUENCY ANALYSIS
% Plot the time-frequency distributions on a separate figure
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

hTFfig = figure(2);

% Define variables/flags required in TF analysis with their default values
plotStyle = 1; % pcolor, imagesc, raw
spectrumType = 1; % raw, difference
cmin = -1; cmax = 3; % caxis limits
tfMethod = 1; % MTM, MP

% Default MTM params
mtmParams.Fs = 2000;
mtmParams.tapers=[2 3]; % [1 1] case is simply stft with dpss window
mtmParams.trialave=0;
mtmParams.err=0;
mtmParams.pad=-1;

movingWin = [0.4 0.01];

% Default MP parameters
numAtomsMP = 100;

% Default HHT parameters
Nstd = 0; % for emd case as default
NE = 1; % for emd case as default
gaussFtr = 7;



%__________________________________________________________________________
% The tf plots panel for the 6 chosen parameters
% ---------------tfPlotsPanel----------------------------------------------
%--------------------------------------------------------------------------

xGap = 0.01; yGap = 0.01;
tfPlotsHeight = 0.14; tfPlotsGap =0.02;
tfPanelWidth = 0.98; tfPanelHeight = 0.81;
hTFPlotsPanel = uipanel('Title','TF Plots','fontSize', fontSizeLarge, ...
    'Unit','Normalized','Position',[xGap yGap tfPanelWidth tfPanelHeight]);

startXPos = 2*xGap; startYPos = 2*yGap; % start from left bottom corner


varyGridTFPos = [startXPos startYPos tfPanelWidth-2*xGap tfPanelHeight-5*yGap];

if strncmp(protocolName,'CRS',3)
    hVaryParamTFPlot = getPlotHandles(length(getValsUnique(gabor7, param7)),length(getValsUnique(gabor8, param8)),varyGridTFPos,gapSmall);

elseif strncmp(protocolName,'GRF',3)
    hVaryParamTFPlot = getPlotHandles(length(getValsUniqueGRF(param7)),length(getValsUniqueGRF(param8)),varyGridTFPos,gapSmall);
end

%__________________________________________________________________________
% Parameters selection panels for TF methods
% ---------------tfParamsPanel----------------------------------------------
%--------------------------------------------------------------------------

xGap = 0.002;
tfParamsHeight = 1-(tfPanelHeight + yGap); tfParamsGap = tfPlotsGap;
tfParamsWidth = 0.4 - xGap; tfParamsPanelHeight = tfParamsHeight+tfPlotsGap;
tfParamsPanelxStart = (1-tfParamsWidth);
tfParamsPanelyStart = (1-tfParamsPanelHeight);

% hTFParamsPanel = uipanel('Title','TF Parameters','fontSize', fontSizeLarge, ...
%     'Unit','Normalized','Position',[tfParamsPanelxStart tfParamsPanelyStart tfParamsWidth tfParamsPanelHeight]);

textWidth = 0.24;
textHeight = 0.18;

% ====================TF Params Panel=====================================
hTFParamsPanel = uipanel('Title','MTM, MP, HHT Parameters','fontSize', fontSizeMedium, ...
    'Unit','Normalized','Position',[(tfParamsPanelxStart) (tfParamsPanelyStart) tfParamsWidth/2 tfParamsPanelHeight]);

% Parameters to divide the panel into rows and columns
xPanel1 = 0; xPanel2 = 0.5;

% Fs, sampling frequency
uicontrol('Parent',hTFParamsPanel,'Unit','Normalized', ...
    'Position',[xPanel1 (1-textHeight) textWidth textHeight], ...
    'Style','text','String','Fs','FontSize',fontSizeSmall);

hMTMFs = uicontrol('Parent',hTFParamsPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [(xPanel1+textWidth) (1-textHeight) textWidth textHeight], ...
    'Style','edit','String',mtmParams.Fs,'FontSize',fontSizeTiny,'Callback',{@resetMTMParams_Callback}); % initialize with default Fs = 2000

%------------------MTM parameters-----------------------
% Tapers TW
uicontrol('Parent',hTFParamsPanel,'Unit','Normalized', ...
    'Position',[xPanel1 (1-2*textHeight-tfParamsGap) textWidth textHeight], ...
    'Style','text','String','TW','FontSize',fontSizeSmall);

hMTMTapersTW = uicontrol('Parent',hTFParamsPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [(xPanel1+textWidth) (1-2*textHeight-tfParamsGap) textWidth textHeight], ...
    'Style','edit','String',mtmParams.tapers(1),'FontSize',fontSizeTiny,'Callback',{@resetMTMParams_Callback}); % initialize TW = 1

% Tapers 'k'
uicontrol('Parent',hTFParamsPanel,'Unit','Normalized', ...
    'Position',[xPanel1 (1-3*textHeight-2*tfParamsGap) textWidth textHeight], ...
    'Style','text','String','k','FontSize',fontSizeSmall);

hMTMTapersK = uicontrol('Parent',hTFParamsPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [(xPanel1+textWidth) (1-3*textHeight-2*tfParamsGap) textWidth textHeight], ...
    'Style','edit','String',mtmParams.tapers(2),'FontSize',fontSizeTiny,'Callback',{@resetMTMParams_Callback}); % initialize k = 1

% Window length
uicontrol('Parent',hTFParamsPanel,'Unit','Normalized', ...
    'Position',[xPanel1 (1-4*textHeight-3*tfParamsGap) textWidth textHeight], ...
    'Style','text','String','wLen','FontSize',fontSizeSmall);

hMTMwLen = uicontrol('Parent',hTFParamsPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [(xPanel1+textWidth) (1-4*textHeight-3*tfParamsGap) textWidth textHeight], ...
    'Style','edit','String',movingWin(1),'FontSize',fontSizeTiny,'Callback',{@resetMTMParams_Callback}); % initialize wLen = 0.1 s

% Window translation step
uicontrol('Parent',hTFParamsPanel,'Unit','Normalized', ...
    'Position',[xPanel1 (1-5*textHeight-4*tfParamsGap) textWidth textHeight], ...
    'Style','text','String','wStep','FontSize',fontSizeSmall);
hMTMwStep = uicontrol('Parent',hTFParamsPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [(xPanel1+textWidth) (1-5*textHeight-4*tfParamsGap) textWidth textHeight], ...
    'Style','edit','String',movingWin(2),'FontSize',fontSizeTiny,'Callback',{@resetMTMParams_Callback}); % initialize wStep = 0.01 s

%----MP parameters-----
% number of atoms for reconstruction
uicontrol('Parent',hTFParamsPanel,'Unit','Normalized', ...
    'Position',[xPanel2 (1-textHeight) textWidth textHeight], ...
    'Style','text','String','numAtoms MP','FontSize',fontSizeSmall);
hMPnumAtoms = uicontrol('Parent',hTFParamsPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [(xPanel2+textWidth) (1-textHeight) textWidth textHeight], ...
    'Style','edit','String',numAtomsMP,'FontSize',fontSizeTiny,'Callback',{@resetMPParams_Callback}); % initialize numAtoms = 100


%----HHT parameters-----
uicontrol('Parent',hTFParamsPanel,'Unit','Normalized', ...
    'Position',[xPanel2 (1-3*textHeight) 2*textWidth textHeight], ...
    'Style','text','String','HHT parameters','FontSize',fontSizeSmall);
% Noise ratio
uicontrol('Parent',hTFParamsPanel,'Unit','Normalized', ...
    'Position',[xPanel2 (1-4*textHeight) textWidth textHeight], ...
    'Style','text','String','Noise ratio','FontSize',fontSizeTiny);
hHHTnstd = uicontrol('Parent',hTFParamsPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [(xPanel2+textWidth) (1-4*textHeight) textWidth textHeight], ...
    'Style','edit','String',Nstd,'FontSize',fontSizeTiny,'Callback',{@resetHHTParams_Callback}); % initialize Nstd = 0


% Ensemble number
uicontrol('Parent',hTFParamsPanel,'Unit','Normalized', ...
    'Position',[xPanel2 (1-5*textHeight) textWidth textHeight], ...
    'Style','text','String','Ensemble number','FontSize',fontSizeTiny);
hHHTne = uicontrol('Parent',hTFParamsPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [(xPanel2+textWidth) (1-5*textHeight) textWidth textHeight], ...
    'Style','edit','String',NE,'FontSize',fontSizeTiny,'Callback',{@resetHHTParams_Callback}); % initialize Nstd = 0


% Gaussian filtering
uicontrol('Parent',hTFParamsPanel,'Unit','Normalized', ...
    'Position',[xPanel2 (1-6*textHeight) textWidth textHeight], ...
    'Style','text','String','Gaussian filter','FontSize',fontSizeTiny);
hHHTgaussFtr = uicontrol('Parent',hTFParamsPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [(xPanel2+textWidth) (1-6*textHeight) textWidth textHeight], ...
    'Style','edit','String',gaussFtr,'FontSize',fontSizeTiny,'Callback',{@resetHHTParams_Callback}); % initialize Nstd = 0


%========= Plotting settings ==============================
hTFPlotSettingsPanel = uipanel('Title','TF Plot settings','fontSize', fontSizeMedium, ...
    'Unit','Normalized','Position',[(tfParamsPanelxStart+(tfParamsWidth/2)+xGap) (tfParamsPanelyStart) ((tfParamsWidth/2)-xGap) tfParamsPanelHeight]);

% pcolor or imagesc or raw
uicontrol('Parent',hTFPlotSettingsPanel,'Unit','Normalized', ...
    'Position',[xPanel1 (1-textHeight) 2*textWidth textHeight], ...
    'Style','text','String','Plot style','FontSize',fontSizeSmall);
hTFPlotStyle = uicontrol('Parent',hTFPlotSettingsPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [(xPanel1+2*textWidth) (1-textHeight) 2*textWidth textHeight], ...
    'Style','pop','String','pcolor|imagesc|line','FontSize',fontSizeTiny,'Callback',{@resetTFSettings_Callback}); % initialize style = pcolor

% spectrum type: raw or difference (change from baseline)
uicontrol('Parent',hTFPlotSettingsPanel,'Unit','Normalized', ...
    'Position',[xPanel1 (1-2*textHeight-tfParamsGap) 2*textWidth textHeight], ...
    'Style','text','String','Spectrum type','FontSize',fontSizeSmall);
hTFSpectrumType = uicontrol('Parent',hTFPlotSettingsPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [(xPanel1+2*textWidth) (1-2*textHeight-tfParamsGap) 2*textWidth textHeight], ...
    'Style','pop','String','raw|difference','FontSize',fontSizeTiny,'Callback',{@resetTFSettings_Callback}); % initialize style = pcolor


% color axis settings
uicontrol('Parent',hTFPlotSettingsPanel,'Unit','Normalized', ...
    'Position',[xPanel1 (1-3*textHeight-2*tfParamsGap) textWidth textHeight], ...
    'Style','text','String','cmin','FontSize',fontSizeSmall);
hTFcmin = uicontrol('Parent',hTFPlotSettingsPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [(xPanel1+textWidth) (1-3*textHeight-2*tfParamsGap) textWidth textHeight], ...
    'Style','edit','String',cmin,'FontSize',fontSizeTiny,'Callback',{@rescaleCaxis_Callback}); % initialize cmin = -1

uicontrol('Parent',hTFPlotSettingsPanel,'Unit','Normalized', ...
    'Position',[xPanel2 (1-3*textHeight-2*tfParamsGap) textWidth textHeight], ...
    'Style','text','String','cmax','FontSize',fontSizeSmall);
hTFcmax = uicontrol('Parent',hTFPlotSettingsPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [(xPanel2+textWidth) (1-3*textHeight-2*tfParamsGap) textWidth textHeight], ...
    'Style','edit','String',cmax,'FontSize',fontSizeTiny,'Callback',{@rescaleCaxis_Callback}); % initialize cmax = 3


% Method used
uicontrol('Parent',hTFPlotSettingsPanel,'Unit','Normalized', ...
    'Position',[xPanel1 (1-5*textHeight-4*tfParamsGap) 2*textWidth 1.5*textHeight], ...
    'Style','text','String','Method','FontSize',fontSizeMedium);
hTFPlotMethod = uicontrol('Parent',hTFPlotSettingsPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [(xPanel1+2*textWidth) (1-5*textHeight-4*tfParamsGap) 2*textWidth 1.5*textHeight], ...
    'Style','pop','String','MTM|MP|HHT','FontSize',fontSizeTiny,'Callback',{@resetTFSettings_Callback}); % initialize method = MTM

% ------------ Reference Selection Panel ----------------------

tfRefPanelxStart = 0.05; tfRefWidth = 0.4;

% Parameters to divide the panel into rows and columns
xPanel1 = 0; xPanel2 = 0.5;

hRefPanel = uipanel('Title','Reference Selection','fontSize', fontSizeMedium, ...
    'Unit','Normalized','Position',[(tfRefPanelxStart) (tfParamsPanelyStart) tfRefWidth/2 tfParamsPanelHeight]);

% reference type: single,bipolar,average,csd
refType = 1; % default value - single reference
uicontrol('Parent',hRefPanel,'Unit','Normalized', ...
    'Position',[xPanel1 (1-textHeight) 2*textWidth textHeight], ...
    'Style','text','String','REF scheme:','FontSize',fontSizeSmall);
hRefType = uicontrol('Parent',hRefPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [(xPanel1+2*textWidth) (1-textHeight) textWidth textHeight], ...
    'Style','pop','String','single|bipolar|average|csd','FontSize',fontSizeTiny,'Callback',{@resetRefType_Callback});

% Analog channel
uicontrol('Parent',hRefPanel,'Unit','Normalized', ...
    'Position',[xPanel1 (1-2*textHeight-tfParamsGap) 2*textWidth textHeight],...
    'Style','text','String','Analog Channel','FontSize',fontSizeTiny);
hRefAnalogChannel1 = uicontrol('Parent',hRefPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [(xPanel1+2*textWidth) (1-2*textHeight-tfParamsGap) textWidth textHeight], ...
    'Style','popup','String',analogChannelStringList,'FontSize',fontSizeTiny);
analogChannelStringList2 = ['none|' analogChannelStringList];
hRefAnalogChannel2 = uicontrol('Parent',hRefPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [(xPanel1+3*textWidth) (1-2*textHeight-tfParamsGap) textWidth textHeight], ...
    'Style','popup','String',analogChannelStringList2,'FontSize',fontSizeTiny','Callback',{@diff_Callback});

% Combine electrodes
combineString = 'none|all|occipital|parietal|temporal|frontal|central';
uicontrol('Parent',hRefPanel,'Unit','Normalized', ...
    'Position',[xPanel1 (1-3*textHeight-tfParamsGap) 2*textWidth textHeight],...
    'Style','text','String','Combine','FontSize',fontSizeTiny);
hRefCombineChannel = uicontrol('Parent',hRefPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [(xPanel1+2*textWidth) (1-3*textHeight-tfParamsGap) textWidth textHeight], ...
    'Style','popup','String',combineString,'FontSize',fontSizeTiny);

% Difference Channel
uicontrol('Parent',hRefPanel,'Unit','Normalized', ...
    'Position',[xPanel1 (1-4*textHeight-tfParamsGap) 2*textWidth textHeight],...
    'Style','text','String','Differ','FontSize',fontSizeTiny);
hRefDiffChannel = uicontrol('Parent',hRefPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [(xPanel1+2*textWidth) (1-4*textHeight-tfParamsGap) textWidth textHeight], ...
    'Style','popup','String',combineString,'FontSize',fontSizeTiny);

% --------------- REF Plot Callback function -----------------------

uicontrol('Parent',hRefPanel,'Unit','Normalized', ...
    'Position',[0 0 1 textHeight], ...
    'Style','pushbutton','String','plot(REF)','FontSize',fontSizeMedium, ...
    'Callback',{@plotRefData_Callback});

%==========================================================================


% --------- TF Plot callback button ---------------------------
uicontrol('Parent',hTFPlotSettingsPanel,'Unit','Normalized', ...
    'Position',[0 0 1 textHeight], ...
    'Style','pushbutton','String','plot','FontSize',fontSizeMedium, ...
    'Callback',{@plotTFData_Callback});


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% functions
    function plotData_Callback(~,~)
        
        if strncmp(protocolName,'GRF',3)
            a=get(hAzimuth,'val');
            e=get(hElevation,'val');
            s=get(hSigma,'val');
            f=get(hSpatialFreq,'val');
            o=get(hOrientation,'val');
            c=get(hContrast,'val');
            t=get(hTemporalFreq,'val');
            r=[];
            p=[];
        elseif strncmp(protocolName,'CRS',3)
            a=cell2mat(get(hAzimuth,'val'));
            e=cell2mat(get(hElevation,'val'));
            s=cell2mat(get(hSigma,'val'));
            f=cell2mat(get(hSpatialFreq,'val'));
            o=cell2mat(get(hOrientation,'val'));
            c=cell2mat(get(hContrast,'val'));
            t=cell2mat(get(hTemporalFreq,'val'));
            r=cell2mat(get(hRadius,'val'));
            p=cell2mat(get(hSpatialPhase,'val'));
        end
        
        analysisType = get(hAnalysisType,'val');
        plotColor = colorNames(get(hChooseColor,'val'));
        BLMin = str2double(get(hBaselineMin,'String'));
        BLMax = str2double(get(hBaselineMax,'String'));
        STMin = str2double(get(hStimPeriodMin,'String'));
        STMax = str2double(get(hStimPeriodMax,'String'));
        holdOnState = get(hHoldOn,'val');
        notchData = get(hNotchData, 'val');
        
        plotLineWidth = linewidths(get(hChooseWidth,'val')); % Vinay - for plot line width
        
        fBandLow = str2double(get(hfBandLow,'String'));
        fBandHigh = str2double(get(hfBandHigh,'String'));
        
        % Load TF Parameters
        % MTM
        mtmParams.Fs = str2double(get(hMTMFs,'String'));
        mtmParams.tapers(1) = str2double(get(hMTMTapersTW,'String'));
        mtmParams.tapers(2) = str2double(get(hMTMTapersK,'String'));
        movingWin(1) = str2double(get(hMTMwLen,'String'));
        movingWin(2) = str2double(get(hMTMwStep,'String'));
        mtmParams.trialave = 1;
        mtmParams.err=0;
        mtmParams.pad=-1; % no padding
        
        if strncmp(protocolName,'CRS',3)
            hVaryParamPlot = getPlotHandles(length(getValsUnique(gabor7, param7)),length(getValsUnique(gabor8, param8)),varyGridPos,gapSmall);
        elseif strncmp(protocolName,'GRF',3)
            hVaryParamPlot = getPlotHandles(length(getValsUniqueGRF(param7)),length(getValsUniqueGRF(param8)),varyGridPos,gapSmall);
        end
        
        % Vinay - get analogChannel info freshly here from the first figure
        % selections. Otherwise the string array contains values from the
        % reference section if that has been previously selected
        [analogChannelsStored,~,~,analogInputNums] = loadlfpInfo(folderLFP);
        [~,analogChannelStringArray] = getAnalogStringFromValues(analogChannelsStored,analogInputNums);



        % Read analog channel/electrode 1
        analogChannelPos = get(hAnalogChannel,'val');
        analogChannelString = analogChannelStringArray{analogChannelPos};

        % [Vinay] - read the analog channel 2 string for bipolar case
        analogChannelPos2 = get(hAnalogChannel2,'val');
        analogChannelString2 = ('none');
        if analogChannelPos2 ~= 1
            analogChannelString2 = analogChannelStringArray{analogChannelPos2-1};
        end


        plotLFPDataVaryParameters1Channel(hVaryParamPlot,analogChannelString,analogChannelString2,a,e,s,f,o,c,t,r,p,gabor7, param7, gabor8, param8, folderLFP,...
            analysisType,timeVals,plotColor,BLMin,BLMax,STMin,STMax,folderName, protocolNumber, notchData,useBipolar,plotSEM,holdOnState,plotLineWidth,...
            protocolName,useAllBadTrials,mtmParams,movingWin,fBandLow,fBandHigh);

        if analogChannelPos<=length(analogChannelsStored)
            channelNumber = analogChannelsStored(analogChannelPos);
        else
            channelNumber = 0;
        end


        if analysisType==1 || analysisType==4 % ERP or Band Power (vs time)
            xMin = str2double(get(hStimMin,'String'));
            xMax = str2double(get(hStimMax,'String'));
        else
            xMin = str2double(get(hFFTMin,'String'));
            xMax = str2double(get(hFFTMax,'String'));
        end
        
        rescaleData(hVaryParamPlot,xMin,xMax,getYLims(hVaryParamPlot));
        
        showElectrodeLocations(electrodeGridPos,channelNumber,plotColor,hElectrodes,holdOnState,0,gridType,subjectName,gridLayout);
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function rescaleY_Callback(~,~)

        analysisType = get(hAnalysisType,'val');
        plotStyle = get(hTFPlotStyle,'val');
        
        if analysisType==1 || analysisType==4 % ERP or Band Power (vs time)
            xMin = str2double(get(hStimMin,'String'));
            xMax = str2double(get(hStimMax,'String'));
        else
            xMin = str2double(get(hFFTMin,'String'));
            xMax = str2double(get(hFFTMax,'String'));
        end

        yLims = [str2double(get(hYMin,'String')) str2double(get(hYMax,'String'))];

        
        rescaleData(hVaryParamPlot,xMin,xMax,yLims);
        
        % [Vinay] - rescale the TF plots
        figure(2);
        
        if plotStyle ~= 3 % pcolor, imagesc
            xMin = str2double(get(hStimMin,'String'));
            xMax = str2double(get(hStimMax,'String'));
            yLims = [str2double(get(hFFTMin,'String')) str2double(get(hFFTMax,'String'))];
        else % line
            xMin = str2double(get(hFFTMin,'String'));
            xMax = str2double(get(hFFTMax,'String'));
            yLims = getYLims(hVaryParamTFPlot);
        end
        
        
        rescaleData(hVaryParamTFPlot,xMin,xMax,yLims);
        
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function rescaleData_Callback(~,~)

        analysisType = get(hAnalysisType,'val');

        if analysisType==1 || analysisType==4 % ERP or spikes
            xMin = str2double(get(hStimMin,'String'));
            xMax = str2double(get(hStimMax,'String'));
        else    
            xMin = str2double(get(hFFTMin,'String'));
            xMax = str2double(get(hFFTMax,'String'));
        end
        
        rescaleData(hVaryParamPlot,xMin,xMax,getYLims(hVaryParamPlot));
        
        % [Vinay] - rescale the TF plots
        figure(2);
        
        if plotStyle ~= 3 % pcolor, imagesc
            xMin = str2double(get(hStimMin,'String'));
            xMax = str2double(get(hStimMax,'String'));
            yLims = [str2double(get(hFFTMin,'String')) str2double(get(hFFTMax,'String'))];
        else % line
            xMin = str2double(get(hFFTMin,'String'));
            xMax = str2double(get(hFFTMax,'String'));
            yLims = getYLims(hVaryParamTFPlot);
        end
        
        
        rescaleData(hVaryParamTFPlot,xMin,xMax,yLims);
        
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function holdOn_Callback(source,~)
        holdOnState = get(source,'Value');
          
        holdOnGivenPlotHandle(hVaryParamPlot,holdOnState);
        
        % TF plots
        holdOnGivenPlotHandle(hVaryParamTFPlot,holdOnState);
        
        if holdOnState
            set(hElectrodes,'Nextplot','add');
        else
            set(hElectrodes,'Nextplot','replace');
        end

        function holdOnGivenPlotHandle(plotHandles,holdOnState)
            
            [numRows,numCols] = size(plotHandles);
            if holdOnState
                for i=1:numRows
                    for j=1:numCols
                        set(plotHandles(i,j),'Nextplot','add');

                    end
                end
            else
                for i=1:numRows
                    for j=1:numCols
                        set(plotHandles(i,j),'Nextplot','replace');
                    end
                end
            end
        end 
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function cla_Callback(~,~)
        
        claGivenPlotHandle(hVaryParamPlot);
        
        % TF plots        
        claGivenPlotHandle(hVaryParamTFPlot);
        
        
        function claGivenPlotHandle(plotHandles)
            [numRows,numCols] = size(plotHandles);
            for i=1:numRows
                for j=1:numCols
                    cla(plotHandles(i,j));
                end
            end
        end
    end
    
    function notch_Callback(hObject,~,~)
       notchData = get(hObject,'val');
    end

    function useDefParams_Callback(hObject,~,~)
       useDefParams = get(hObject,'val');
       if useDefParams == 1
           switch protocolNumber
               case 1 % Ring protocol
                   set(hGabor1,'val',4);% 4 => S
                   set(hGabor2,'val',4);
                   set(hGabor3,'val',4);
                   set(hGabor4,'val',4);
                   set(hGabor5,'val',3); % 3 => R
                   set(hGabor6,'val',1); % 1 => NA
                   set(hGabor7,'val',1); % 1 => NA
                   
                   set(hParam1,'val',5); % 5 => SF
                   set(hParam2,'val',6); % 6 => ori
                   set(hParam3,'val',7); % 7 => cont
                   set(hParam4,'val',8); % 8 => TF 
                   set(hParam5,'val',9); % 9 => rad 
                   set(hParam6,'val',1); % 1 => NA
                   set(hParam7,'val',1); % 1 => NA
               case 2 % Contrast Ring Protocol
                   set(hGabor1,'val',4);% 4 => S
                   set(hGabor2,'val',4);
                   set(hGabor3,'val',4);
                   set(hGabor4,'val',4);
                   set(hGabor5,'val',3); % 3 => R
                   set(hGabor6,'val',3);
                   set(hGabor7,'val',1); % 1 => NA
                   
                   set(hParam1,'val',5); % 5 => SF
                   set(hParam2,'val',6); % 6 => ori
                   set(hParam3,'val',7); % 7 => cont
                   set(hParam4,'val',8); % 8 => TF 
                   set(hParam5,'val',7); % 7 => cont (Ring) 
                   set(hParam6,'val',9); % 9 => rad
                   set(hParam7,'val',1); % 1 => NA
               case 3 % Dual Contrast Protocol
                   set(hGabor1,'val',2);% 2 => C
                   set(hGabor2,'val',2);
                   set(hGabor3,'val',2);
                   set(hGabor4,'val',2);
                   set(hGabor5,'val',3); % 3 => R
                   set(hGabor6,'val',3);
                   set(hGabor7,'val',1); % 1 => NA
                   
                   set(hParam1,'val',5); % 5 => SF
                   set(hParam2,'val',6); % 6 => ori
                   set(hParam3,'val',7); % 7 => cont
                   set(hParam4,'val',8); % 8 => TF 
                   set(hParam5,'val',7); % 7 => cont (Ring) 
                   set(hParam6,'val',9); % 9 => rad
                   set(hParam7,'val',1); % 1 => NA
               case 4 % Dual Orientation Protocol
                   set(hGabor1,'val',2);% 2 => C
                   set(hGabor2,'val',2);
                   set(hGabor3,'val',2);
                   set(hGabor4,'val',2);
                   set(hGabor5,'val',3); % 3 => R
                   set(hGabor6,'val',3);
                   set(hGabor7,'val',1); % 1 => NA
                   
                   set(hParam1,'val',5); % 5 => SF
                   set(hParam2,'val',6); % 6 => ori
                   set(hParam3,'val',7); % 7 => cont
                   set(hParam4,'val',8); % 8 => TF 
                   set(hParam5,'val',6); % 6 => ori (Ring) 
                   set(hParam6,'val',9); % 9 => rad
                   set(hParam7,'val',1); % 1 => NA
               case 5 % Dual Phase Protocol
                   set(hGabor1,'val',2);% 2 => C
                   set(hGabor2,'val',2);
                   set(hGabor3,'val',2);
                   set(hGabor4,'val',2);
                   set(hGabor5,'val',3); % 3 => R
                   set(hGabor6,'val',3);
                   set(hGabor7,'val',1); % 1 => NA
                   
                   set(hParam1,'val',5); % 5 => SF
                   set(hParam2,'val',6); % 6 => ori
                   set(hParam3,'val',7); % 7 => cont
                   set(hParam4,'val',8); % 8 => TF 
                   set(hParam5,'val',10); % 10 => phase (Ring) 
                   set(hParam6,'val',9); % 9 => rad
                   set(hParam7,'val',1); % 1 => NA
               case 6 %Orientation Ring Protocol
                   set(hGabor1,'val',4);% 4 => S
                   set(hGabor2,'val',4);
                   set(hGabor3,'val',4);
                   set(hGabor4,'val',4);
                   set(hGabor5,'val',3); % 3 => R
                   set(hGabor6,'val',3);
                   set(hGabor7,'val',1); % 1 => NA
                   
                   set(hParam1,'val',5); % 5 => SF
                   set(hParam2,'val',6); % 6 => ori
                   set(hParam3,'val',7); % 7 => cont
                   set(hParam4,'val',8); % 8 => TF 
                   set(hParam5,'val',6); % 6 => ori (Ring) 
                   set(hParam6,'val',9); % 9 => rad
                   set(hParam7,'val',1); % 1 => NA
               case 7 %Phase Ring Protocol
                   set(hGabor1,'val',4);% 4 => S
                   set(hGabor2,'val',4);
                   set(hGabor3,'val',4);
                   set(hGabor4,'val',4);
                   set(hGabor5,'val',3); % 3 => R
                   set(hGabor6,'val',3);
                   set(hGabor7,'val',1); % 1 => NA
                   
                   set(hParam1,'val',5); % 5 => SF
                   set(hParam2,'val',6); % 6 => ori
                   set(hParam3,'val',7); % 7 => cont
                   set(hParam4,'val',8); % 8 => TF 
                   set(hParam5,'val',7); % 10 => phase (Ring) 
                   set(hParam6,'val',9); % 9 => rad
                   set(hParam7,'val',1); % 1 => NA
               case 8 %Drifting Ring Protocol
                   set(hGabor1,'val',4);% 4 => S
                   set(hGabor2,'val',4);
                   set(hGabor3,'val',4);
                   set(hGabor4,'val',4);
                   set(hGabor5,'val',3); % 3 => R
                   set(hGabor6,'val',3);
                   set(hGabor7,'val',1); % 1 => NA
                   
                   set(hParam1,'val',5); % 5 => SF
                   set(hParam2,'val',6); % 6 => ori
                   set(hParam3,'val',7); % 7 => cont
                   set(hParam4,'val',8); % 8 => TF 
                   set(hParam5,'val',8); % 8 => TF (Ring) 
                   set(hParam6,'val',9); % 9 => rad
                   set(hParam7,'val',1); % 1 => NA
               case 9 %Cross Orientation Protocol
                   set(hGabor1,'val',4); % 4 => S
                   set(hGabor2,'val',4);
                   set(hGabor3,'val',4);
                   set(hGabor4,'val',3); % 3 => R
                   set(hGabor5,'val',3);
                   set(hGabor6,'val',3);
                   set(hGabor7,'val',1); % 1 => NA
                   
                   set(hParam1,'val',5); % 5 => SF
                   set(hParam2,'val',6); % 6 => ori
                   set(hParam3,'val',7); % 7 => cont
                   set(hParam4,'val',8); % 8 => TF 
                   set(hParam5,'val',9); % 9 => rad 
                   set(hParam6,'val',9); % 9 => rad (Ring)
                   set(hParam7,'val',1); % 1 => NA
               case 10 %Annulus Fixed Protocol
                   set(hGabor1,'val',4);% 2 => S
                   set(hGabor2,'val',4);
                   set(hGabor3,'val',4);
                   set(hGabor4,'val',4);
                   set(hGabor5,'val',4);
                   set(hGabor6,'val',3); % 3 => R
                   set(hGabor7,'val',1); % 1 => NA
                   
                   set(hParam1,'val',5); % 5 => SF
                   set(hParam2,'val',6); % 6 => ori
                   set(hParam3,'val',7); % 7 => cont
                   set(hParam4,'val',8); % 8 => TF 
                   set(hParam5,'val',9); % 9 => rad 
                   set(hParam6,'val',9); % 9 => rad (Ring)
                   set(hParam7,'val',1); % 1 => NA
               otherwise
                   set(hGabor1,'val',4); % 4 => S
                   set(hGabor2,'val',4);
                   set(hGabor3,'val',3); % 3 => R
                   set(hGabor4,'val',3);
                   set(hGabor5,'val',2); % 2 => C
                   set(hGabor6,'val',2);
                   set(hGabor7,'val',1); % 1 => NA
                   
                   set(hParam1,'val',7); % 7 => cont (surr)
                   set(hParam2,'val',6); % 6 => ori
                   set(hParam3,'val',7); % 7 => cont (ring)
                   set(hParam4,'val',9); % 9 => rad 
                   set(hParam5,'val',7); % 7 => cont (centre)
                   set(hParam6,'val',6); % 6 => ori
                   set(hParam7,'val',1); % 1 => NA
           end
       end
       adjFactor = 5;
        % If sel is C, then value = 2 (NA = 1, C = 2, R = 3, S = 4) => gabor1 = 5 - 2 = 3. If R then 2, if S then 1
       gabor7 = adjFactor - get(hGabor7,'val');
       gabor8 = adjFactor - get(hGabor8,'val');
       
       param7 = get(hParam7, 'val');
       param8 = get(hParam8, 'val');
       
       
       if strncmp(protocolName,'CRS',3)
            hVaryParamPlot = getPlotHandles(length(getValsUnique(gabor7, param7)),length(getValsUnique(gabor8, param8)),varyGridPos,gapSmall);
       elseif strncmp(protocolName,'GRF',3)
            hVaryParamPlot = getPlotHandles(length(getValsUniqueGRF(param7)),length(getValsUniqueGRF(param8)),varyGridPos,gapSmall);
       end
       
       % [Vinay] - reset the TF grid as well
       figure(2);
       
       if strncmp(protocolName,'CRS',3)
            hVaryParamTFPlot = getPlotHandles(length(getValsUnique(gabor7, param7)),length(getValsUnique(gabor8, param8)),varyGridTFPos,gapSmall);
       elseif strncmp(protocolName,'GRF',3)
            hVaryParamTFPlot = getPlotHandles(length(getValsUniqueGRF(param7)),length(getValsUniqueGRF(param8)),varyGridTFPos,gapSmall);
       end
        
    end

    function resetVaryParams_Callback(~,~)
        adjFactor = 5;
        % If sel is C, then value = 2 (NA = 1, C = 2, R = 3, S = 4) => gabor1 = 5 - 2 = 3. If R then 2, if S then 1

        if strncmp(protocolName,'CRS',3)
            gabor7 = adjFactor - get(hGabor7,'val');
            gabor8 = adjFactor - get(hGabor8,'val');
        end

        param7 = get(hParam7, 'val');
        param8 = get(hParam8, 'val');
        
        
        if strncmp(protocolName,'CRS',3)
            hVaryParamPlot = getPlotHandles(length(getValsUnique(gabor7, param7)),length(getValsUnique(gabor8, param8)),varyGridPos,gapSmall);
        elseif strncmp(protocolName,'GRF',3)
            hVaryParamPlot = getPlotHandles(length(getValsUniqueGRF(param7)),length(getValsUniqueGRF(param8)),varyGridPos,gapSmall);
        end
        
        % [Vinay] - reset the TF grid as well
        figure(2);
        
        if strncmp(protocolName,'CRS',3)
            hVaryParamTFPlot = getPlotHandles(length(getValsUnique(gabor7, param7)),length(getValsUnique(gabor8, param8)),varyGridTFPos,gapSmall);
        elseif strncmp(protocolName,'GRF',3)
            hVaryParamTFPlot = getPlotHandles(length(getValsUniqueGRF(param7)),length(getValsUniqueGRF(param8)),varyGridTFPos,gapSmall);
        end
        
    end

    function valsUnique = getValsUnique(gaborNumber, paramNumber)
        if gaborNumber > 3
            valsUnique = [];
        else
            switch (paramNumber-1)
                case 1
                    valsUnique = aValsUnique{gaborNumber};
                case 2
                    valsUnique = eValsUnique{gaborNumber};
                case 3
                    valsUnique = sValsUnique{gaborNumber};
                case 4
                    valsUnique = fValsUnique{gaborNumber};
                case 5
                    valsUnique = oValsUnique{gaborNumber};
                case 6
                    valsUnique = cValsUnique{gaborNumber};
                case 7
                    valsUnique = tValsUnique{gaborNumber};
                case 8
                    valsUnique = rValsUnique{gaborNumber};
                case 9
                    valsUnique = pValsUnique{gaborNumber};
                otherwise
                    valsUnique = [];
            end
        end
    end

%--------------------------------------------------------------------------

function valsUnique = getValsUniqueGRF(paramNumber)
        
        switch (paramNumber-1)
            case 1
                valsUnique = aValsUnique;
            case 2
                valsUnique = eValsUnique;
            case 3
                valsUnique = sValsUnique;
            case 4
                valsUnique = fValsUnique;
            case 5
                valsUnique = oValsUnique;
            case 6
                valsUnique = cValsUnique;
            case 7
                valsUnique = tValsUnique;
            otherwise
                valsUnique = [];
        end
end

%--------------------------------------------------------------------------

    function bipolar_Callback(hObject,~,~)
       if (get(hObject,'val') ~= 1)
           useBipolar = 1;
       else
           useBipolar = 0;
       end
    end
        

%---------------TF analysis Functions----------------------------
    function resetMTMParams_Callback(~,~)
        mtmParams.Fs = str2double(get(hMTMFs,'String'));
        mtmParams.tapers(1) = str2double(get(hMTMTapersTW,'String'));
        mtmParams.tapers(2) = str2double(get(hMTMTapersK,'String'));
        
        movingWin(1) = str2double(get(hMTMwLen,'String'));
        movingWin(2) = str2double(get(hMTMwStep,'String'));
    end

    function resetMPParams_Callback(~,~)
        numAtomsMP = str2double(get(hMPnumAtoms,'String'));
    end

    function resetTFSettings_Callback(~,~)
        plotStyle = get(hTFPlotStyle,'val');
        spectrumType = get(hTFSpectrumType,'val');
        
        cmin = str2double(get(hTFcmin,'String'));
        cmax = str2double(get(hTFcmax,'String'));
        
        tfMethod = get(hTFPlotMethod,'val');
    end

    function rescaleCaxis_Callback(~,~)
        cmin = str2double(get(hTFcmin,'String'));
        cmax = str2double(get(hTFcmax,'String'));
        
        
        set(hVaryParamTFPlot,'cLim',[cmin cmax]);
        
    end
        
    
    function resetHHTParams_Callback(~,~)
       
        disp('Resetting HHT parameters...');
        Nstd = str2double(get(hHHTnstd,'String'));
        NE = str2double(get(hHHTne,'String'));
        gaussFtr = str2double(get(hHHTgaussFtr,'String'));
        
    end
        

    %----main TF plotting function----------------
    function plotTFData_Callback(~,~)
        
        if strncmp(protocolName,'GRF',3)
            a=get(hAzimuth,'val');
            e=get(hElevation,'val');
            s=get(hSigma,'val');
            f=get(hSpatialFreq,'val');
            o=get(hOrientation,'val');
            c=get(hContrast,'val');
            t=get(hTemporalFreq,'val');
            r=[];
            p=[];
        elseif strncmp(protocolName,'CRS',3)
            a=cell2mat(get(hAzimuth,'val'));
            e=cell2mat(get(hElevation,'val'));
            s=cell2mat(get(hSigma,'val'));
            f=cell2mat(get(hSpatialFreq,'val'));
            o=cell2mat(get(hOrientation,'val'));
            c=cell2mat(get(hContrast,'val'));
            t=cell2mat(get(hTemporalFreq,'val'));
            r=cell2mat(get(hRadius,'val'));
            p=cell2mat(get(hSpatialPhase,'val'));
        end
        
        analysisType = get(hAnalysisType,'val');
        plotColor = colorNames(get(hChooseColor,'val'));
        BLMin = str2double(get(hBaselineMin,'String'));
        BLMax = str2double(get(hBaselineMax,'String'));
        STMin = str2double(get(hStimPeriodMin,'String'));
        STMax = str2double(get(hStimPeriodMax,'String'));
        holdOnState = get(hHoldOn,'val');
        notchData = get(hNotchData, 'val');
        
        plotLineWidth = linewidths(get(hChooseWidth,'val')); % Vinay - for plot line width
        
        
        % Load TF Parameters
        % MTM
        mtmParams.Fs = str2double(get(hMTMFs,'String'));
        mtmParams.tapers(1) = str2double(get(hMTMTapersTW,'String'));
        mtmParams.tapers(2) = str2double(get(hMTMTapersK,'String'));
        movingWin(1) = str2double(get(hMTMwLen,'String'));
        movingWin(2) = str2double(get(hMTMwStep,'String'));
        mtmParams.trialave = 1;
        mtmParams.err=0;
        mtmParams.pad=-1; % no padding
        
        % MP
        numAtomsMP = str2double(get(hMPnumAtoms,'String'));
        
        % HHT
        Nstd = str2double(get(hHHTnstd,'String'));
        NE = str2double(get(hHHTne,'String'));
        gaussFtr = str2double(get(hHHTgaussFtr,'String'));
        
        % Plot settings
        plotStyle = get(hTFPlotStyle,'val');
        spectrumType = get(hTFSpectrumType,'val');
        
        cmin = str2double(get(hTFcmin,'String'));
        cmax = str2double(get(hTFcmax,'String'));
        
        tfMethod = get(hTFPlotMethod,'val');
        
        
        
        if strncmp(protocolName,'CRS',3)
            hVaryParamTFPlot = getPlotHandles(length(getValsUnique(gabor7, param7)),length(getValsUnique(gabor8, param8)),varyGridTFPos,gapSmall);
        elseif strncmp(protocolName,'GRF',3)
            hVaryParamTFPlot = getPlotHandles(length(getValsUniqueGRF(param7)),length(getValsUniqueGRF(param8)),varyGridTFPos,gapSmall);
        end
        
        
        % Vinay - get analogChannel info freshly here from the first figure
        % selections. Otherwise the string array contains values from the
        % reference section if that has been previously selected
        [analogChannelsStored,~,~,analogInputNums] = loadlfpInfo(folderLFP);
        [~,analogChannelStringArray] = getAnalogStringFromValues(analogChannelsStored,analogInputNums);

        
        % Read analog channel/electrode 1
        analogChannelPos = get(hAnalogChannel,'val');
        analogChannelString = analogChannelStringArray{analogChannelPos};

        % [Vinay] - read the analog channel 2 string for bipolar case
        analogChannelPos2 = get(hAnalogChannel2,'val');
        analogChannelString2 = ('none');
        if analogChannelPos2 ~= 1
            analogChannelString2 = analogChannelStringArray{analogChannelPos2-1};
        end
        
        
        tfplotLFPDataVaryParameters1Channel(hVaryParamTFPlot,analogChannelString,analogChannelString2,a,e,s,f,o,c,t,r,p,gabor7, param7, gabor8, param8,folderLFP,...
            timeVals,plotColor,BLMin,BLMax,STMin,STMax,folderName, protocolNumber, notchData,useBipolar,...
            tfMethod,mtmParams,movingWin,numAtomsMP,plotStyle,spectrumType,cmin,cmax, holdOnState, saveMPFlag,loadProtocolNumber,plotLineWidth,...
            Nstd,NE,gaussFtr,saveHHTFlag,protocolName,useAllBadTrials);
        
    
    end


% -------------------------------------------------------------------------
% Reference scheme related callbacks

function diff_Callback(hObject,~,~)
   if (get(hObject,'val') ~= 1)
       useDiff = 1;
   else
       useDiff = 0;
   end
end

function resetRefType_Callback(~,~)
    
    refType = get(hRefType,'val');
       
    if refType==2 && existsBipolarData
        
        [electrodesStoredPair,electrodesStoredCz] = loadBipolarlfpInfo(folderBipolar);
        analogChannelsStored = [electrodesStoredPair,electrodesStoredCz];
    else
        [analogChannelsStored,~,~,~] = loadlfpInfo(folderLFP);
    end

    [analogChannelStringList,analogChannelStringArray] = getAnalogStringFromValues(analogChannelsStored,[]);
    
    
    set(hRefAnalogChannel1,'String',analogChannelStringList);
    set(hRefAnalogChannel2,'String',['none|' analogChannelStringList]);
    
    
end

% -------------------------------------------------------------------------

 function plotRefData_Callback(~,~)
        
        if strncmp(protocolName,'GRF',3)
            a=get(hAzimuth,'val');
            e=get(hElevation,'val');
            s=get(hSigma,'val');
            f=get(hSpatialFreq,'val');
            o=get(hOrientation,'val');
            c=get(hContrast,'val');
            t=get(hTemporalFreq,'val');
            r=[];
            p=[];
        elseif strncmp(protocolName,'CRS',3)
            a=cell2mat(get(hAzimuth,'val'));
            e=cell2mat(get(hElevation,'val'));
            s=cell2mat(get(hSigma,'val'));
            f=cell2mat(get(hSpatialFreq,'val'));
            o=cell2mat(get(hOrientation,'val'));
            c=cell2mat(get(hContrast,'val'));
            t=cell2mat(get(hTemporalFreq,'val'));
            r=cell2mat(get(hRadius,'val'));
            p=cell2mat(get(hSpatialPhase,'val'));
        end
        
        analysisType = get(hAnalysisType,'val');
        plotColor = colorNames(get(hChooseColor,'val'));
        BLMin = str2double(get(hBaselineMin,'String'));
        BLMax = str2double(get(hBaselineMax,'String'));
        STMin = str2double(get(hStimPeriodMin,'String'));
        STMax = str2double(get(hStimPeriodMax,'String'));
        holdOnState = get(hHoldOn,'val');
        notchData = get(hNotchData, 'val');
        
        plotLineWidth = linewidths(get(hChooseWidth,'val')); % Vinay - for plot line width
        
        % Load TF Parameters
        % MTM
        mtmParams.Fs = str2double(get(hMTMFs,'String'));
        mtmParams.tapers(1) = str2double(get(hMTMTapersTW,'String'));
        mtmParams.tapers(2) = str2double(get(hMTMTapersK,'String'));
        movingWin(1) = str2double(get(hMTMwLen,'String'));
        movingWin(2) = str2double(get(hMTMwStep,'String'));
        mtmParams.trialave = 1;
        mtmParams.err=0;
        mtmParams.pad=-1; % no padding
        
        % MP
        numAtomsMP = str2double(get(hMPnumAtoms,'String'));
        
        % HHT
        Nstd = str2double(get(hHHTnstd,'String'));
        NE = str2double(get(hHHTne,'String'));
        gaussFtr = str2double(get(hHHTgaussFtr,'String'));
        
        % Plot settings
        plotStyle = get(hTFPlotStyle,'val');
        spectrumType = get(hTFSpectrumType,'val');
        
        cmin = str2double(get(hTFcmin,'String'));
        cmax = str2double(get(hTFcmax,'String'));
        
        tfMethod = get(hTFPlotMethod,'val');
        
        
        if strncmp(protocolName,'CRS',3)
            hVaryParamTFPlot = getPlotHandles(length(getValsUnique(gabor7, param7)),length(getValsUnique(gabor8, param8)),varyGridTFPos,gapSmall);
        elseif strncmp(protocolName,'GRF',3)
            hVaryParamTFPlot = getPlotHandles(length(getValsUniqueGRF(param7)),length(getValsUniqueGRF(param8)),varyGridTFPos,gapSmall);
        end
        
        
        % Vinay - load analogChannel info freshly for reference options
        refType = get(hRefType,'val');
        if refType==2 && existsBipolarData
            [electrodesStoredPair,electrodesStoredCz] = loadBipolarlfpInfo(folderBipolar);
            analogChannelsStored = [electrodesStoredPair,electrodesStoredCz];
        else
            [analogChannelsStored,~,~,~] = loadlfpInfo(folderLFP);
        end

        [~,analogChannelStringArray] = getAnalogStringFromValues(analogChannelsStored,[]);
        
        
        % Load REF information
        refCombineType = get(hRefCombineChannel,'val');
        refDiffType = get(hRefDiffChannel,'val');
        numElectrodes = length(analogChannelsStored);
        
        analogChannelPos = get(hRefAnalogChannel1,'val');
        analogChannelString = analogChannelStringArray{analogChannelPos};

        % [Vinay] - read the analog channel 2 string for bipolar case
        analogChannelPos2 = get(hRefAnalogChannel2,'val');
        analogChannelString2 = ('none');
        if analogChannelPos2 ~= 1
            analogChannelString2 = analogChannelStringArray{analogChannelPos2-1};
        end

        
        reftfplotLFPDataVaryParameters1Channel(hVaryParamTFPlot,analogChannelString,analogChannelString2,a,e,s,f,o,c,t,r,p,gabor7, param7, gabor8, param8,folderLFP,...
            timeVals,plotColor,BLMin,BLMax,STMin,STMax,folderName, protocolNumber, notchData,useDiff,...
            tfMethod,mtmParams,movingWin,numAtomsMP,plotStyle,spectrumType,cmin,cmax, holdOnState, saveMPFlag,loadProtocolNumber,plotLineWidth,...
            Nstd,NE,gaussFtr,saveHHTFlag,protocolName,useAllBadTrials,...
            refType,refCombineType,refDiffType,numElectrodes,intersectTrials,...
            takeLogTrial,existsBipolarData,contralateralNeighbour);
        
    
    end


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotLFPDataVaryParameters1Channel(plotHandles,channelString,analogChannelString2,a,e,s,f,o,c,t,r,p,gaborNum1,paramNum1,gaborNum2,paramNum2,folderLFP,...
analysisType,timeVals,plotColor,BLMin,BLMax,STMin,STMax,folderName, protocolNumber, notchData,useBipolar,plotSEM,holdOnState,plotLineWidth,...
protocolName,useAllBadTrials,mtmParams,movingWin,fBandLow,fBandHigh)

folderExtract = fullfile(folderName,'extractedData');
folderSegment = fullfile(folderName,'segmentedData');

titleFontSize = 10;

timeForComputation = [40 100]/1000; % ms
freqForComputation = [30 80]; % Hz

[parameterCombinations,aValsUnique,eValsUnique,sValsUnique,...
    fValsUnique,oValsUnique,cValsUnique,tValsUnique,rValsUnique,pValsUnique] = loadParameterCombinations(folderExtract);

numRows = size(plotHandles,1);
numCols = size(plotHandles,2);

% colors for the conjunct plots
% colorsList = colormap(winter(numRows*numCols));
% colorsList = colormap(default(numRows*numCols));

% Get the data
clear signal analogData analogDataNotched
load(fullfile(folderLFP,channelString));

% Get bad trials
badTrialFile = fullfile(folderSegment,'badTrials.mat');
existsBadTrialFile = 0;
if ~exist(badTrialFile,'file')
    disp('Bad trial file does not exist...');
    badTrials=[];
else
    [allBadTrials, badTrials, nameElec] = loadBadTrials(badTrialFile);
    disp([num2str(length(badTrials)) ' bad trials']);
    existsBadTrialFile = 1;
end


% Choosing the good position trials
if strncmp(protocolName,'GRF',3)
    % [Vinay] - repeat the set of parameters for each gabor depending on the
    % number of columns to be drawn
    aList = repmat(a,numRows,numCols);
    eList = repmat(e,numRows,numCols);
    sList = repmat(s,numRows,numCols);
    fList = repmat(f,numRows,numCols);
    oList = repmat(o,numRows,numCols);
    cList = repmat(c,numRows,numCols);
    tList = repmat(t,numRows,numCols);


        % [Vinay] decide the row parameter based on paramNum1
            for row = 1:numRows
                switch (paramNum1-1)
                    case 1
                        aList(row,:) = row;                    
                        % Every new row takes a new value of param1.
                        % So assign the 'row' value in the list above
                        titleParam1 = 'Azi: ';
                        titleList1 = aValsUnique;
                    case 2
                        eList(row,:) = row; 
                        titleParam1 = 'Ele: ';
                        titleList1 = eValsUnique;
                    case 3
                        sList(row,:) = row; 
                        titleParam1 = 'Sigma: ';
                        titleList1 = sValsUnique;
                    case 4
                        fList(row,:) = row; 
                        titleParam1 = 'SF: ';
                        titleList1 = fValsUnique;
                    case 5
                        oList(row,:) = row; 
                        titleParam1 = 'Ori: ';
                        titleList1 = oValsUnique;
                    case 6
                        cList(row,:) = row; 
                        titleParam1 = 'Contr: ';
                        titleList1 = cValsUnique;
                    case 7
                        tList(row,:) = row; 
                        titleParam1 = 'TF: ';
                        titleList1 = tValsUnique;
                    otherwise
                        disp('No particular gabor or parameter selected');
                end
            end

            % [Vinay] decide the column parameter based on paramNum2
            for row = 1:numRows
                switch (paramNum2-1)
                    case 1
                        aList(row,:) = 1:length(aValsUnique);
                        % For every row we now have to assign indices
                        % incrementing along the columns.
                        % Basically the column entries go from 1 to 
                        % length{nValsUnique}
                        titleParam2 = 'Azi: ';
                        titleList2 = aValsUnique;
                    case 2
                        eList(row,:) = 1:length(eValsUnique); 
                        titleParam2 = 'Ele: ';
                        titleList2 = eValsUnique;
                    case 3
                        sList(row,:) = 1:length(sValsUnique);
                        titleParam2 = 'Sigma: ';
                        titleList2 = sValsUnique;
                    case 4
                        fList(row,:) = 1:length(fValsUnique);
                        titleParam2 = 'SF: ';
                        titleList2 = fValsUnique;
                    case 5
                        oList(row,:) = 1:length(oValsUnique);
                        titleParam2 = 'Ori: ';
                        titleList2 = oValsUnique;
                    case 6
                        cList(row,:) = 1:length(cValsUnique);
                        titleParam2 = 'Contr: ';
                        titleList2 = cValsUnique;
                    case 7
                        tList(row,:) = 1:length(tValsUnique);
                        titleParam2 = 'TF: ';
                        titleList2 = tValsUnique;
                    otherwise
                        disp('No particular gabor or parameter selected');
                end
            end


    % [Vinay] - Get the lengths of indices in parameterCombinations

    aLen = length(aValsUnique);
    eLen = length(eValsUnique);
    sLen = length(sValsUnique);
    fLen = length(fValsUnique);
    oLen = length(oValsUnique);
    cLen = length(cValsUnique);
    tLen = length(tValsUnique);

    % If more than one value, then length is one greater for all the values
    % together
    if (aLen> 1)           aLen=aLen+1;                    end
    if (eLen> 1)           eLen=eLen+1;                    end
    if (sLen> 1)           sLen=sLen+1;                    end
    if (fLen> 1)           fLen=fLen+1;                    end
    if (oLen> 1)           oLen=oLen+1;                    end
    if (cLen> 1)           cLen=cLen+1;                    end
    if (tLen> 1)           tLen=tLen+1;                    end

elseif strncmp(protocolName,'CRS',3)
    
    % [Vinay] - repeat the set of parameters for each gabor depending on the
    % number of columns to be drawn
    aList = repmat(a,numRows,numCols);
    eList = repmat(e,numRows,numCols);
    sList = repmat(s,numRows,numCols);
    fList = repmat(f,numRows,numCols);
    oList = repmat(o,numRows,numCols);
    cList = repmat(c,numRows,numCols);
    tList = repmat(t,numRows,numCols);
    rList = repmat(r,numRows,numCols);
    pList = repmat(p,numRows,numCols);

            % [Vinay] decide the row parameter based on paramNum1
            for row = 1:numRows
                switch (paramNum1-1)
                    case 1
                        aList(gaborNum1+((row-1)*3),:) = row;                    
                        % along the rows, we get three consecutive entries for the 3 
                        % gabors. Therefore to go to the next index of the same gabor,
                        % we go ahead by 3 indices
                        % For gabor1 every new row takes a new value of param1.
                        % So assign the 'row' value in the list above
                        titleParam1 = 'Azi: ';
                        titleList1 = aValsUnique{gaborNum1};
                    case 2
                        eList(gaborNum1+((row-1)*3),:) = row; 
                        titleParam1 = 'Ele: ';
                        titleList1 = eValsUnique{gaborNum1};
                    case 3
                        sList(gaborNum1+((row-1)*3),:) = row; 
                        titleParam1 = 'Sigma: ';
                        titleList1 = sValsUnique{gaborNum1};
                    case 4
                        fList(gaborNum1+((row-1)*3),:) = row; 
                        titleParam1 = 'SF: ';
                        titleList1 = fValsUnique{gaborNum1};
                    case 5
                        oList(gaborNum1+((row-1)*3),:) = row; 
                        titleParam1 = 'Ori: ';
                        titleList1 = oValsUnique{gaborNum1};
                    case 6
                        cList(gaborNum1+((row-1)*3),:) = row; 
                        titleParam1 = 'Contr: ';
                        titleList1 = cValsUnique{gaborNum1};
                    case 7
                        tList(gaborNum1+((row-1)*3),:) = row; 
                        titleParam1 = 'TF: ';
                        titleList1 = tValsUnique{gaborNum1};
                    case 8
                        rList(gaborNum1+((row-1)*3),:) = row; 
                        titleParam1 = 'Rad: ';
                        titleList1 = rValsUnique{gaborNum1};
                    case 9
                        pList(gaborNum1+((row-1)*3),:) = row; 
                        titleParam1 = 'SPhs: ';
                        titleList1 = pValsUnique{gaborNum1};
                    otherwise
                        disp('No particular gabor or parameter selected');
                end
            end


            % [Vinay] decide the column parameter based on paramNum2
            for row = 1:numRows
                switch (paramNum2-1)
                    case 1
                        aList(gaborNum2+((row-1)*3),:) = 1:length(aValsUnique{gaborNum2});
                        % along the rows, we get three consecutive entries for the 3 
                        % gabors. Therefore to go to the next index of the same gabor,
                        % we go ahead by 3 indices
                        % Besides for every row we now have to assign indices
                        % incrementing along the columns.
                        % Basically the column entries go from 1 to 
                        % length{nValsUnique}
                        titleParam2 = 'Azi: ';
                        titleList2 = aValsUnique{gaborNum2};
                    case 2
                        eList(gaborNum2+((row-1)*3),:) = 1:length(eValsUnique{gaborNum2}); 
                        titleParam2 = 'Ele: ';
                        titleList2 = eValsUnique{gaborNum2};
                    case 3
                        sList(gaborNum2+((row-1)*3),:) = 1:length(sValsUnique{gaborNum2});
                        titleParam2 = 'Sigma: ';
                        titleList2 = sValsUnique{gaborNum2};
                    case 4
                        fList(gaborNum2+((row-1)*3),:) = 1:length(fValsUnique{gaborNum2});
                        titleParam2 = 'SF: ';
                        titleList2 = fValsUnique{gaborNum2};
                    case 5
                        oList(gaborNum2+((row-1)*3),:) = 1:length(oValsUnique{gaborNum2});
                        titleParam2 = 'Ori: ';
                        titleList2 = oValsUnique{gaborNum2};
                    case 6
                        cList(gaborNum2+((row-1)*3),:) = 1:length(cValsUnique{gaborNum2});
                        titleParam2 = 'Contr: ';
                        titleList2 = cValsUnique{gaborNum2};
                    case 7
                        tList(gaborNum2+((row-1)*3),:) = 1:length(tValsUnique{gaborNum2});
                        titleParam2 = 'TF: ';
                        titleList2 = tValsUnique{gaborNum2};
                    case 8
                        rList(gaborNum2+((row-1)*3),:) = 1:length(rValsUnique{gaborNum2});
                        titleParam2 = 'Rad: ';
                        titleList2 = rValsUnique{gaborNum2};
                    case 9
                        pList(gaborNum2+((row-1)*3),:) = 1:length(pValsUnique{gaborNum2});
                        titleParam2 = 'SPhs: ';
                        titleList2 = pValsUnique{gaborNum2};
                    otherwise
                        disp('No particular gabor or parameter selected');
                end
            end


    % [Vinay] - Get the lengths of indices in parameterCombinations
    for i = 1:3
        aLen(i) = length(aValsUnique{i});
        eLen(i) = length(eValsUnique{i});
        sLen(i) = length(sValsUnique{i});
        fLen(i) = length(fValsUnique{i});
        oLen(i) = length(oValsUnique{i});
        cLen(i) = length(cValsUnique{i});
        tLen(i) = length(tValsUnique{i});
        pLen(i) = length(pValsUnique{i}); % Vinay - for CRS
        rLen(i) = length(rValsUnique{i}); % Vinay - for CRS

        % If more than one value, then length is one greater for all the values
        % together
        if (aLen(i)> 1)           aLen(i)=aLen(i)+1;                    end
        if (eLen(i)> 1)           eLen(i)=eLen(i)+1;                    end
        if (sLen(i)> 1)           sLen(i)=sLen(i)+1;                    end
        if (fLen(i)> 1)           fLen(i)=fLen(i)+1;                    end
        if (oLen(i)> 1)           oLen(i)=oLen(i)+1;                    end
        if (cLen(i)> 1)           cLen(i)=cLen(i)+1;                    end
        if (tLen(i)> 1)           tLen(i)=tLen(i)+1;                    end
        if (pLen(i)> 1)           pLen(i)=pLen(i)+1;                    end % Vinay - for CRS
        if (rLen(i)> 1)           rLen(i)=rLen(i)+1;                    end % Vinay - for CRS
    end
    
end

% Main loop
computationVals=zeros(1,numCols);
for k=1:numRows
    for j = 1:numCols
        
        if strncmp(protocolName,'GRF',3)
            clear goodPos
            goodPos = parameterCombinations{aList(k,j),eList(k,j),sList(k,j),fList(k,j),oList(k,j),cList(k,j),tList(k,j)};
    %         goodPos = setdiff(goodPos,badTrials);

            %------------------
            % Vinay - select good trials as per the electrode(s)
    %         useAllBadTrials = 1;
            if useAllBadTrials && existsBadTrialFile
                elecIndex1 = (strcmp(channelString,nameElec) == 1);

                if ~useBipolar
                    elecBadTrials = allBadTrials{elecIndex1};
                    disp(['No. of Bad trials for ' channelString ': ' num2str(length(elecBadTrials))]);
                else
                    elecIndex2 = (strcmp(analogChannelString2,nameElec) == 1);
                    elecBadTrials = unique([allBadTrials{elecIndex1} allBadTrials{elecIndex2}]);
                    disp(['No. of Bad trials for ' channelString ' and ' analogChannelString2 ': ' num2str(length(elecBadTrials))]);
                end
                
                goodPos = setdiff(goodPos,elecBadTrials);
            else
                goodPos = setdiff(goodPos,badTrials);
                
            end
            %-------------------

        elseif strncmp(protocolName,'CRS',3)
            
            clear goodPos
        %     goodPos = parameterCombinations{aList(j),eList(j),sList(j),fList(j),oList(j),cList(j),tList(j)};
        %     goodPos = setdiff(goodPos,badTrials);

            % [Vinay] - goodPos will be obtained by taking the goodPos for all
            % the three gabors and finding the common trials in them
            clear pos1 pos2 pos3
            switch protocolNumber
                       case 1 % Ring protocol
                           pos1 = parameterCombinations{aList(1+(k-1)*3,j),eList(1+(k-1)*3,j),sList(1+(k-1)*3,j),fList(1+(k-1)*3,j),oList(1+(k-1)*3,j),cList(1+(k-1)*3,j),tList(1+(k-1)*3,j),pList(1+(k-1)*3,j),rList(1+(k-1)*3,j),1}; 
                           % S params define C & S
                           % S i.e. kGabor0 corresponds to row number 1,4,7,..
                           % and so on. Hence it is 1+(k-1)*3
                           pos2 = parameterCombinations{aLen(2),eLen(2),sLen(2),fLen(2),oLen(2),cLen(2),tLen(2),pLen(2),rList(2+(k-1)*3,j),2}; % for R - only rad is imp
                           % R i.e. kGabor1 corresponds to row number 2,5,8,..
                           % and so on. Hence it is 2+(k-1)*3

                           % [Vinay] -  C is matched to S and therefore the entries
                           % for gabor3 i.e. can be ignored. Otherwise this can 
                           % give spurious results if the param values for C are 
                           % not set equal to those for S. Hence pos3 is the full
                           % set here
                           pos3 = parameterCombinations{aLen(3),eLen(3),sLen(3),fLen(3),oLen(3),cLen(3),tLen(3),pLen(3),rList(3+(k-1)*3,j),3}; % for C
                           % also keeping r variable here, just in case you
                           % have multiple r values

                       case 2 % Contrast Ring Protocol
                           pos1 = parameterCombinations{aList(1+(k-1)*3,j),eList(1+(k-1)*3,j),sList(1+(k-1)*3,j),fList(1+(k-1)*3,j),oList(1+(k-1)*3,j),cList(1+(k-1)*3,j),tList(1+(k-1)*3,j),pList(1+(k-1)*3,j),rList(1+(k-1)*3,j),1};
                           % S params define C & S
                           % S i.e. kGabor0 corresponds to row number 1,4,7,..
                           % and so on. Hence it is 1+(k-1)*3
                           pos2 = parameterCombinations{aLen(2),eLen(2),sLen(2),fLen(2),oLen(2),cList(2+(k-1)*3,j),tLen(2),pLen(2),rList(2+(k-1)*3,j),2}; % for R - only rad & cont
                           % R i.e. kGabor1 corresponds to row number 2,5,8,..
                           % and so on. Hence it is 2+(k-1)*3

                           % [Vinay] -  C is matched to S and therefore the entries
                           % for gabor3 i.e. can be ignored. Otherwise this can 
                           % give spurious results if the param values for C are 
                           % not set equal to those for S. Hence pos3 is the full
                           % set here
                           pos3 = parameterCombinations{aLen(3),eLen(3),sLen(3),fLen(3),oLen(3),cLen(3),tLen(3),pLen(3),rList(3+(k-1)*3,j),3}; % for C
                           % also keeping r variable here, just in case you
                           % have multiple r values

                       case 3 % Dual Contrast Protocol
                           pos1 = parameterCombinations{aLen(1),eLen(1),sLen(1),fLen(1),oLen(1),cLen(1),tLen(1),pLen(1),rLen(1),1}; % for C
                           % [Vinay] - S is hidden in this case and only C & R
                           % define the stimuli. pos1 is therefore the full set
                           % here
                           pos2 = parameterCombinations{aLen(2),eLen(2),sLen(2),fLen(2),oLen(2),cList(2+(k-1)*3,j),tLen(2),pLen(2),rList(2+(k-1)*3,j),2}; % for R - only rad & cont
                           % [Vinay] -  R params - cont & rad. Rest are matched to
                           % C's params
                           pos3 = parameterCombinations{aList(3+(k-1)*3,j),eList(3+(k-1)*3,j),sList(3+(k-1)*3,j),fList(3+(k-1)*3,j),oList(3+(k-1)*3,j),cList(3+(k-1)*3,j),tList(3+(k-1)*3,j),pList(3+(k-1)*3,j),rList(3+(k-1)*3,j),3};
                           % for C. C is the defining gabor here

                       case 4 % Dual Orientation Protocol
                           pos1 = parameterCombinations{aLen(1),eLen(1),sLen(1),fLen(1),oLen(1),cLen(1),tLen(1),pLen(1),rLen(1),1}; % for C
                           % [Vinay] - S is hidden in this case and only C & R
                           % define the stimuli. pos1 is therefore the full set
                           % here
                           pos2 = parameterCombinations{aLen(2),eLen(2),sLen(2),fLen(2),oList(2+(k-1)*3,j),cLen(2),tLen(2),pLen(2),rList(2+(k-1)*3,j),2}; % for R - only rad & ori
                           % [Vinay] -  R params - ori & rad. Rest are matched to
                           % C's params
                           pos3 = parameterCombinations{aList(3+(k-1)*3,j),eList(3+(k-1)*3,j),sList(3+(k-1)*3,j),fList(3+(k-1)*3,j),oList(3+(k-1)*3,j),cList(3+(k-1)*3,j),tList(3+(k-1)*3,j),pList(3+(k-1)*3,j),rList(3+(k-1)*3,j),3};
                           % for C. C is the defining gabor here

                       case 5 % Dual Phase Protocol
                           pos1 = parameterCombinations{aLen(1),eLen(1),sLen(1),fLen(1),oLen(1),cLen(1),tLen(1),pLen(1),rLen(1),1}; % for C
                           % [Vinay] - S is hidden in this case and only C & R
                           % define the stimuli. pos1 is therefore the full set
                           % here
                           pos2 = parameterCombinations{aLen(2),eLen(2),sLen(2),fLen(2),oLen(2),cLen(2),tLen(2),pList(2+(k-1)*3,j),rList(2+(k-1)*3,j),2}; % for R - only rad & p
                           % [Vinay] -  R params - p & rad. Rest are matched to
                           % C's params
                           pos3 = parameterCombinations{aList(3+(k-1)*3,j),eList(3+(k-1)*3,j),sList(3+(k-1)*3,j),fList(3+(k-1)*3,j),oList(3+(k-1)*3,j),cList(3+(k-1)*3,j),tList(3+(k-1)*3,j),pList(3+(k-1)*3,j),rList(3+(k-1)*3,j),3};
                           % for C. C is the defining gabor here

                       case 6 %Orientation Ring Protocol
                           pos1 = parameterCombinations{aList(1+(k-1)*3,j),eList(1+(k-1)*3,j),sList(1+(k-1)*3,j),fList(1+(k-1)*3,j),oList(1+(k-1)*3,j),cList(1+(k-1)*3,j),tList(1+(k-1)*3,j),pList(1+(k-1)*3,j),rList(1+(k-1)*3,j),1};
                           % S params define C & S
                           % S i.e. kGabor0 corresponds to row number 1,4,7,..
                           % and so on. Hence it is 1+(k-1)*3
                           pos2 = parameterCombinations{aLen(2),eLen(2),sLen(2),fLen(2),oList(2+(k-1)*3,j),cLen(2),tLen(2),pLen(2),rList(2+(k-1)*3,j),2}; % for R - only rad & ori
                           % [Vinay] -  C is matched to S and therefore the entries
                           % for gabor3 i.e. can be ignored. Otherwise this can 
                           % give spurious results if the param values for C are 
                           % not set equal to those for S. Hence pos3 is the full
                           % set here
                           pos3 = parameterCombinations{aLen(3),eLen(3),sLen(3),fLen(3),oLen(3),cLen(3),tLen(3),pLen(3),rList(3+(k-1)*3,j),3}; % for C
                           % also keeping r variable here, just in case you
                           % have multiple r values

                       case 7 %Phase Ring Protocol
                           pos1 = parameterCombinations{aList(1+(k-1)*3,j),eList(1+(k-1)*3,j),sList(1+(k-1)*3,j),fList(1+(k-1)*3,j),oList(1+(k-1)*3,j),cList(1+(k-1)*3,j),tList(1+(k-1)*3,j),pList(1+(k-1)*3,j),rList(1+(k-1)*3,j),1};
                           % S params define C & S
                           % S i.e. kGabor0 corresponds to row number 1,4,7,..
                           % and so on. Hence it is 1+(k-1)*3
                           pos2 = parameterCombinations{aLen(2),eLen(2),sLen(2),fLen(2),oLen(2),cLen(2),tLen(2),pList(2+(k-1)*3,j),rList(2+(k-1)*3,j),2}; % for R - only rad & p
                           % [Vinay] -  C is matched to S and therefore the entries
                           % for gabor3 i.e. can be ignored. Otherwise this can 
                           % give spurious results if the param values for C are 
                           % not set equal to those for S. Hence pos3 is the full
                           % set here
                           pos3 = parameterCombinations{aLen(3),eLen(3),sLen(3),fLen(3),oLen(3),cLen(3),tLen(3),pLen(3),rList(3+(k-1)*3,j),3}; % for C
                           % also keeping r variable here, just in case you
                           % have multiple r values

                       case 8 %Drifting Ring Protocol
                           pos1 = parameterCombinations{aList(1+(k-1)*3,j),eList(1+(k-1)*3,j),sList(1+(k-1)*3,j),fList(1+(k-1)*3,j),oList(1+(k-1)*3,j),cList(1+(k-1)*3,j),tList(1+(k-1)*3,j),pList(1+(k-1)*3,j),rList(1+(k-1)*3,j),1};
                           % S params define C & S
                           % S i.e. kGabor0 corresponds to row number 1,4,7,..
                           % and so on. Hence it is 1+(k-1)*3
                           pos2 = parameterCombinations{aLen(2),eLen(2),sLen(2),fLen(2),oLen(2),cLen(2),tList(2+(k-1)*3,j),pLen(2),rList(2+(k-1)*3,j),2}; % for R - only rad & t
                           % [Vinay] -  C is matched to S and therefore the entries
                           % for gabor3 i.e. can be ignored. Otherwise this can 
                           % give spurious results if the param values for C are 
                           % not set equal to those for S. Hence pos3 is the full
                           % set here
                           pos3 = parameterCombinations{aLen(3),eLen(3),sLen(3),fLen(3),oLen(3),cLen(3),tLen(3),pLen(3),rList(3+(k-1)*3,j),3}; % for C
                           % also keeping r variable here, just in case you
                           % have multiple r values

                       case 9 %Cross Orientation Protocol
                           pos1 = parameterCombinations{aLen(1),eLen(1),sLen(1),fLen(1),oLen(1),cLen(1),tLen(1),pLen(1),rLen(1),1}; % for C
                           % [Vinay] - S is hidden in this case and only C & R
                           % define the stimuli. pos1 is therefore the full set
                           % here
                           pos2 = parameterCombinations{aList(2+(k-1)*3,j),eList(2+(k-1)*3,j),sList(2+(k-1)*3,j),fList(2+(k-1)*3,j),oList(2+(k-1)*3,j),cList(2+(k-1)*3,j),tList(2+(k-1)*3,j),pList(2+(k-1)*3,j),rList(2+(k-1)*3,j),2};
                           % for R. R is one of the defining gabor here
                           pos3 = parameterCombinations{aList(3+(k-1)*3,j),eList(3+(k-1)*3,j),sList(3+(k-1)*3,j),fList(3+(k-1)*3,j),oList(3+(k-1)*3,j),cList(3+(k-1)*3,j),tList(3+(k-1)*3,j),pList(3+(k-1)*3,j),rList(3+(k-1)*3,j),3};
                           % for C. C is one of the defining gabor here

                       case 10 %Annulus Fixed Protocol
                           pos1 = parameterCombinations{aList(1+(k-1)*3,j),eList(1+(k-1)*3,j),sList(1+(k-1)*3,j),fList(1+(k-1)*3,j),oList(1+(k-1)*3,j),cList(1+(k-1)*3,j),tList(1+(k-1)*3,j),pList(1+(k-1)*3,j),rList(1+(k-1)*3,j),1};
                           % S params define C & S
                           % S i.e. kGabor0 corresponds to row number 1,4,7,..
                           % and so on. Hence it is 1+(k-1)*3
                           pos2 = parameterCombinations{aLen(2),eLen(2),sLen(2),fLen(2),oLen(2),cList(2+(k-1)*3,j),tLen(2),pLen(2),rList(2+(k-1)*3,j),2}; % for R - only rad & cont
                           % [Vinay] -  C is matched to S and therefore the entries
                           % for gabor3 i.e. can be ignored. Otherwise this can 
                           % give spurious results if the param values for C are 
                           % not set equal to those for S. However, the rad of C is
                           % still a relevant parameter and has to be taken care of
                           pos3 = parameterCombinations{aLen(3),eLen(3),sLen(3),fLen(3),oLen(3),cLen(3),tLen(3),pLen(3),rList(3+(k-1)*3,j),3}; % for C - only rad

                otherwise
                    pos1 = parameterCombinations{aList(1+(k-1)*3,j),eList(1+(k-1)*3,j),sList(1+(k-1)*3,j),fList(1+(k-1)*3,j),oList(1+(k-1)*3,j),cList(1+(k-1)*3,j),tList(1+(k-1)*3,j),pList(1+(k-1)*3,j),rList(1+(k-1)*3,j),1};
                    % good positions for gabor 1 i.e. S
                    pos2 = parameterCombinations{aList(2+(k-1)*3,j),eList(2+(k-1)*3,j),sList(2+(k-1)*3,j),fList(2+(k-1)*3,j),oList(2+(k-1)*3,j),cList(2+(k-1)*3,j),tList(2+(k-1)*3,j),pList(2+(k-1)*3,j),rList(2+(k-1)*3,j),2};
                    pos3 = parameterCombinations{aList(3+(k-1)*3,j),eList(3+(k-1)*3,j),sList(3+(k-1)*3,j),fList(3+(k-1)*3,j),oList(3+(k-1)*3,j),cList(3+(k-1)*3,j),tList(3+(k-1)*3,j),pList(3+(k-1)*3,j),rList(3+(k-1)*3,j),3};
            end

            goodPos = intersect(pos1,pos2);
            goodPos = intersect(goodPos,pos3);
    %         goodPos = setdiff(goodPos,badTrials);

            %------------------
            % Vinay - select good trials as per the electrode(s)
    %         useAllBadTrials = 1;
            if useAllBadTrials && existsBadTrialFile
                elecIndex1 = (strcmp(channelString,nameElec) == 1);

                if ~useBipolar
                    elecBadTrials = allBadTrials{elecIndex1};
                    disp(['No. of Bad trials for ' channelString ': ' num2str(length(elecBadTrials))]);
                else
                    elecIndex2 = (strcmp(analogChannelString2,nameElec) == 1);
                    elecBadTrials = unique([allBadTrials{elecIndex1} allBadTrials{elecIndex2}]);
                    disp(['No. of Bad trials for ' channelString ' and ' analogChannelString2 ': ' num2str(length(elecBadTrials))]);
                end
                
                goodPos = setdiff(goodPos,elecBadTrials);
            else
                goodPos = setdiff(goodPos,badTrials);
            end
            %-------------------
            
        end

            if isempty(goodPos)
                disp('No entries for this combination..')
            else
                disp(['pos=' num2str(k) ',' num2str(j) ',n=' num2str(length(goodPos))]);

                Fs = round(1/(timeVals(2)-timeVals(1)));
                BLRange = uint16((BLMax-BLMin)*Fs);
                STRange = uint16((STMax-STMin)*Fs);
                BLPos = find(timeVals>=BLMin,1)+ (1:BLRange);
                STPos = find(timeVals>=STMin,1)+ (1:STRange);

                xsBL = 0:1/(BLMax-BLMin):Fs-1/(BLMax-BLMin);
                xsST = 0:1/(STMax-STMin):Fs-1/(STMax-STMin);


                xsComputation = intersect(find(timeVals>=timeForComputation(1)),find(timeVals<timeForComputation(2)));
                freqComputation = intersect(find(xsST>=freqForComputation(1)),find(xsST<=freqForComputation(2)));

                % Vinay - added this notch data check
                if notchData
                    analogData = analogDataNotched;
                end

                if useBipolar
                    analogChannelString1 = channelString; % first electrode selected
                    clear analogData analogDataNotched
                    load(fullfile(folderLFP,analogChannelString1));

                    if notchData
                        analogData1 = analogDataNotched;
                    else
                        analogData1 = analogData; % Vinay - these are the 
                        % freshly loaded values for electrode 1
                    end
                    % Vinay - store the analogData 
                    % for electrode 1 in a separate variable 
                    % otherwise loading the data for 2nd electrode will
                    % overwrite analogData

                    clear analogData analogDataNotched
                    load(fullfile(folderLFP,analogChannelString2));

                    if notchData
                        analogData2 = analogDataNotched;
                    else
                        analogData2 = analogData; % Vinay - these are the 
                        % freshly loaded values for electrode 2
                    end
                    analogData = analogData1 - analogData2;
                end

                if analysisType == 1        % compute ERP
                    clear erp
                    erp = mean(analogData(goodPos,:),1);

                    plot(plotHandles(k,j),timeVals,erp,'color',plotColor,'Linewidth',plotLineWidth);

                    % Vinay - plotting SEM
                    if plotSEM
                        thisPlotColor = get(plot(plotHandles(k,j),timeVals,erp,'color',plotColor),'Color');
                        thisPlotColor(thisPlotColor==0) = 0.75;
                        erpSEM = std(analogData(goodPos,:),[],1)./sqrt(length(goodPos)); % Vinay - SEM = std/sqrt(n)

                        set(plotHandles(k,j),'Nextplot','add');
                        plot(plotHandles(k,j),timeVals,(erp+erpSEM),'color',thisPlotColor);
                        plot(plotHandles(k,j),timeVals,(erp-erpSEM),'color',thisPlotColor);
                        plot(plotHandles(k,j),timeVals,erp,'color',plotColor,'Linewidth',plotLineWidth);

                        if ~holdOnState
                            set(plotHandles(k,j),'Nextplot','replace');
                        end

                    end


                    if paramNum1==6 % Orientation tuning
                        computationVals(j) = abs(min(erp(xsComputation)));
                    end
                    
                    % ---- plot ERPs for all cases together----------
%                     set(hERP,'Nextplot','add');
%                     colorNumber = (k-1)*numCols + j; % there are numRows*numCols numbers and we go in an increasing order
%                     plot(hERP,timeVals,erp,'color',colorsList(colorNumber,:),'Linewidth',plotLineWidth);
%                     hLegERP{colorNumber} = [num2str(k) num2str(j)];


                else

                    fftBL = abs(fft(analogData(goodPos,BLPos),[],2));
                    fftST = abs(fft(analogData(goodPos,STPos),[],2));

                    if analysisType == 2
                        plot(plotHandles(k,j),xsBL,conv2Log(mean(fftBL)),'color',plotColor,'Linewidth',plotLineWidth,'LineStyle','--');
                        set(plotHandles(k,j),'Nextplot','add');
                        plot(plotHandles(k,j),xsST,conv2Log(mean(fftST)),'color',plotColor,'Linewidth',plotLineWidth);

                        if ~holdOnState
                            set(plotHandles(k,j),'Nextplot','replace');
                        end
                    end

                    if analysisType == 3
                        if xsBL == xsST %#ok<BDSCI>
                            plot(plotHandles(k,j),xsBL,conv2Log(mean(fftST))-conv2Log(mean(fftBL)),'color',plotColor,'Linewidth',plotLineWidth);
                        else
                            disp('Choose same baseline and stimulus periods..');
                        end
                    end
                    
                    if analysisType == 4
                        takeLogTrial = 0;
                        showMean = 1;
                        mtmParams.trialave=0;
                        [~,dS,t2,f2] = getSpectrum(analogData(goodPos,:),mtmParams,movingWin,takeLogTrial,BLMin,BLMax,timeVals);
                        
                        freqRange = (f2 >= fBandLow) & (f2 <= fBandHigh);
                        
                        dBandPower = sum(dS(:,freqRange),2);
                        
                        meandBandPower = mean(dS(:,freqRange),2);
                        semBandPower = std(dS(:,freqRange),[],2)./sqrt(size(dS,2));
                        
                        tST = (t2>=STMin) & (t2<=STMax); % stimulus time indices
                        totalGamma(k,j) = sum(dBandPower(tST,:))./(fBandHigh-fBandLow);
                        
                        if ~showMean
                            plot(plotHandles(k,j),t2,dBandPower,'color',plotColor,'Linewidth',plotLineWidth);
                        else
                            thisPlotColor = get(plot(plotHandles(k,j),t2,meandBandPower,'color',plotColor,'Linewidth',plotLineWidth),'Color');
                            thisPlotColor(thisPlotColor==0) = 0.75;
                            
                            set(plotHandles(k,j),'Nextplot','add');
                            plot(plotHandles(k,j),t2,(meandBandPower+semBandPower),'color',thisPlotColor);
                            plot(plotHandles(k,j),t2,(meandBandPower-semBandPower),'color',thisPlotColor);
                            plot(plotHandles(k,j),t2,meandBandPower,'color',plotColor,'Linewidth',plotLineWidth);
                        end
                        
                    end

                end
                
                if strncmp(protocolName,'CRS',3)
                    titleGabor1 = titleGabor(gaborNum1);
                    titleGabor2 = titleGabor(gaborNum2);
                elseif strncmp(protocolName,'GRF',3)
                    titleGabor1 = 'gabor';
                    titleGabor2 = 'gabor';
                end

                % Display title
                if (j==1 && k==1)
                    title(plotHandles(k,j),[titleGabor1 '-' titleParam1 'vs' titleGabor2 '-' titleParam2 ':' num2str(titleList1(k)) ',' num2str(titleList2(j))],'FontSize',titleFontSize);
                else
                    title(plotHandles(k,j),[num2str(titleList1(k)) ',' num2str(titleList2(j))],'FontSize',titleFontSize);
                end
                
                
                % show number of trials for this condition
                text(0.1,0.9,['n = ' num2str(length(goodPos))],'unit','normalized','fontsize',7,'Parent',plotHandles(k,j));
                
            end
            
    end
end

if analysisType == 4
    % plot the gamma trend
    figure;
    hTG = subplot(2,2,1);
    contourf(totalGamma);

    set(hTG,'XTick',1:numCols);
    set(hTG,'XTickLabel',titleList2);

    set(hTG,'YTick',1:numRows);
    set(hTG,'YTickLabel',titleList1);
end

end


function titleGaborN = titleGabor(gaborNum)
    if gaborNum == 1
        titleGaborN = 'S';
    elseif gaborNum == 2
        titleGaborN = 'R';
    elseif gaborNum == 3
        titleGaborN = 'C';
    else
        titleGaborN = 'Wrong gabor number';
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function yLims = getYLims(plotHandles)

[numRows,numCols] = size(plotHandles);
% Initialize
yMin = inf;
yMax = -inf;

for row=1:numRows
    for column=1:numCols
        % get positions
        axis(plotHandles(row,column),'tight');
        tmpAxisVals = axis(plotHandles(row,column));
        if tmpAxisVals(3) < yMin
            yMin = tmpAxisVals(3);
        end
        if tmpAxisVals(4) > yMax
            yMax = tmpAxisVals(4);
        end
    end
end

yLims=[yMin yMax];
end
function rescaleData(plotHandles,xMin,xMax,yLims)

[numRows,numCols] = size(plotHandles);
labelSize=12;
for i=1:numRows
    for j=1:numCols
        axis(plotHandles(i,j),[xMin xMax yLims]);
        if (i==numRows && rem(j,2)==1)
            if j~=1
                set(plotHandles(i,j),'YTickLabel',[],'fontSize',labelSize);
            end
        elseif (rem(i,2)==0 && j==1)
            set(plotHandles(i,j),'XTickLabel',[],'fontSize',labelSize);
        else
            set(plotHandles(i,j),'XTickLabel',[],'YTickLabel',[],'fontSize',labelSize);
        end
    end
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function outString = getStringFromValues(valsUnique,decimationFactor, gaborNum)

if length(valsUnique{gaborNum})==1
    outString = convertNumToStr(valsUnique{gaborNum}(1),decimationFactor);
else
    outString='';
    for i=1:length(valsUnique{gaborNum})
        outString = cat(2,outString,[convertNumToStr(valsUnique{gaborNum}(i),decimationFactor) '|']);
    end
    outString = [outString 'all'];
end

    function str = convertNumToStr(num,f)
        if num > 16384
            num=num-32768;
        end
        str = num2str(num/f);
    end
end

%--------------------------------------------------------------------------

function outString = getStringFromValuesGRF(valsUnique,decimationFactor)

if length(valsUnique)==1
    outString = convertNumToStr(valsUnique(1),decimationFactor);
else
    outString='';
    for i=1:length(valsUnique)
        outString = cat(2,outString,[convertNumToStr(valsUnique(i),decimationFactor) '|']);
    end
    outString = [outString 'all'];
end

    function str = convertNumToStr(num,f)
        if num > 16384
            num=num-32768;
        end
        str = num2str(num/f);
    end
end

%--------------------------------------------------------------------------


function [outString,outArray] = getAnalogStringFromValues(analogChannelsStored,analogInputNums)
outString='';
count=1;
for i=1:length(analogChannelsStored)
    outArray{count} = ['elec' num2str(analogChannelsStored(i))]; %#ok<AGROW>
    outString = cat(2,outString,[outArray{count} '|']);
    count=count+1;
end
if ~isempty(analogInputNums)
    for i=1:length(analogInputNums)
        outArray{count} = ['ainp' num2str(analogInputNums(i))]; %#ok<AGROW>
        outString = cat(2,outString,[outArray{count} '|']);
        count=count+1;
    end
end
end
function outString = getNeuralStringFromValues(neuralChannelsStored,SourceUnitIDs)
outString='';
for i=1:length(neuralChannelsStored)
    outString = cat(2,outString,[num2str(neuralChannelsStored(i)) ', SID ' num2str(SourceUnitIDs(i)) '|']);
end 
end
function [colorString, colorNames] = getColorString

colorNames = 'brkgcmy';
colorString = 'blue|red|black|green|cyan|magenta|yellow';

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%c%%%%%%%%%
%%%%%%%%%%%%c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load Data
function [analogChannelsStored,timeVals,goodStimPos,analogInputNums] = loadlfpInfo(folderLFP) %#ok<*STOUT>
load(fullfile(folderLFP,'lfpInfo'));
if ~exist('analogInputNums','var')
    analogInputNums=[];
end
end
function [neuralChannelsStored,SourceUnitID] = loadspikeInfo(folderSpikes)
fileName = [folderSpikes 'spikeInfo.mat'];
if exist(fileName,'file')
    load(fileName);
else
    neuralChannelsStored=[];
    SourceUnitID=[];
end
end
function [parameterCombinations,aValsUnique,eValsUnique,sValsUnique,...
    fValsUnique,oValsUnique,cValsUnique,tValsUnique, rValsUnique, pValsUnique] = loadParameterCombinations(folderExtract)
% [Vinay] - added rValsUnique and pValsUnique for radius and spatial phase
load(fullfile(folderExtract,'parameterCombinations.mat'));

if ~exist('rValsUnique','var')
    rValsUnique=[];
end

if ~exist('pValsUnique','var')
    pValsUnique=[];
end

if ~exist('sValsUnique','var')
    sValsUnique=rValsUnique;
end

if ~exist('cValsUnique','var')
    cValsUnique=[];
end

if ~exist('tValsUnique','var')
    tValsUnique=[];
end
end
% function stimResults = loadStimResults(folderExtract)
% load ([folderExtract 'stimResults']);
% end
function [allBadTrials, badTrials, nameElec] = loadBadTrials(badTrialFile)
load(badTrialFile);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function protocolNumber = getProtocolNumber(folderExtract)
load (fullfile(folderExtract,'stimResults'));
protocolNumber = stimResults.protocolNumber;
end

function gaborsDisplayed = getGaborsDisplayed(folderExtract)
load (fullfile(folderExtract,'stimResults'));
gaborsDisplayed = stimResults.side;
end

function [electrodesStoredPair,electrodesStoredCz] = loadBipolarlfpInfo(folderBipolar) %#ok<*STOUT>
load(fullfile(folderBipolar,'lfpInfo'));
end

function tfplotLFPDataVaryParameters1Channel(tfplotHandles,channelString,analogChannelString2,a,e,s,f,o,c,t,r,p,gaborNum1,paramNum1,gaborNum2,paramNum2,folderLFP,...
            timeVals,plotColor,BLMin,BLMax,STMin,STMax,folderName, protocolNumber, notchData,useBipolar,...
            tfMethod,mtmParams,movingWin,numAtomsMP,plotStyle,spectrumType,cmin,cmax,holdOnState,saveMPFlag,loadProtocolNumber,plotLineWidth,...
            Nstd,NE,gaussFtr,saveHHTFlag,protocolName,useAllBadTrials)
        
        
        folderExtract = fullfile(folderName,'extractedData');
        folderSegment = fullfile(folderName,'segmentedData');
        mpFolder = fullfile(folderName,'mpAnalysis'); % for using stored MP data
        hhtFolder = fullfile(folderName,'hhtAnalysis');
        
        if ispc
            mpFolderTemp = 'C:\Documents\MP\data\'; % for online calculations
        else
            mpFolderTemp = '/home/vinay/Documents/MP/data/'; % for online calculations
        end


        titleFontSize = 9;

%         timeForComputation = [40 100]/1000; % ms
%         freqForComputation = [40 60]; % Hz

        [parameterCombinations,aValsUnique,eValsUnique,sValsUnique,...
            fValsUnique,oValsUnique,cValsUnique,tValsUnique,rValsUnique,pValsUnique] = loadParameterCombinations(folderExtract);
        
        numRows = size(tfplotHandles,1);
        numCols = size(tfplotHandles,2);


        % Get the data
        clear signal analogData
        load(fullfile(folderLFP,channelString));

        % Get bad trials
        badTrialFile = fullfile(folderSegment,'badTrials.mat');
        existsBadTrialFile = 0;
        if ~exist(badTrialFile,'file')
            disp('Bad trial file does not exist...');
            badTrials=[];
        else
            [allBadTrials, badTrials, nameElec] = loadBadTrials(badTrialFile);
            disp([num2str(length(badTrials)) ' bad trials']);
            existsBadTrialFile = 1;
        end
        
        
        % Choosing the good position trials
if strncmp(protocolName,'GRF',3)
    % [Vinay] - repeat the set of parameters for each gabor depending on the
    % number of columns to be drawn
    aList = repmat(a,numRows,numCols);
    eList = repmat(e,numRows,numCols);
    sList = repmat(s,numRows,numCols);
    fList = repmat(f,numRows,numCols);
    oList = repmat(o,numRows,numCols);
    cList = repmat(c,numRows,numCols);
    tList = repmat(t,numRows,numCols);


        % [Vinay] decide the row parameter based on paramNum1
            for row = 1:numRows
                switch (paramNum1-1)
                    case 1
                        aList(row,:) = row;                    
                        % Every new row takes a new value of param1.
                        % So assign the 'row' value in the list above
                        titleParam1 = 'Azi: ';
                        titleList1 = aValsUnique;
                    case 2
                        eList(row,:) = row; 
                        titleParam1 = 'Ele: ';
                        titleList1 = eValsUnique;
                    case 3
                        sList(row,:) = row; 
                        titleParam1 = 'Sigma: ';
                        titleList1 = sValsUnique;
                    case 4
                        fList(row,:) = row; 
                        titleParam1 = 'SF: ';
                        titleList1 = fValsUnique;
                    case 5
                        oList(row,:) = row; 
                        titleParam1 = 'Ori: ';
                        titleList1 = oValsUnique;
                    case 6
                        cList(row,:) = row; 
                        titleParam1 = 'Contr: ';
                        titleList1 = cValsUnique;
                    case 7
                        tList(row,:) = row; 
                        titleParam1 = 'TF: ';
                        titleList1 = tValsUnique;
                    otherwise
                        disp('No particular gabor or parameter selected');
                end
            end

            % [Vinay] decide the column parameter based on paramNum2
            for row = 1:numRows
                switch (paramNum2-1)
                    case 1
                        aList(row,:) = 1:length(aValsUnique);
                        % For every row we now have to assign indices
                        % incrementing along the columns.
                        % Basically the column entries go from 1 to 
                        % length{nValsUnique}
                        titleParam2 = 'Azi: ';
                        titleList2 = aValsUnique;
                    case 2
                        eList(row,:) = 1:length(eValsUnique); 
                        titleParam2 = 'Ele: ';
                        titleList2 = eValsUnique;
                    case 3
                        sList(row,:) = 1:length(sValsUnique);
                        titleParam2 = 'Sigma: ';
                        titleList2 = sValsUnique;
                    case 4
                        fList(row,:) = 1:length(fValsUnique);
                        titleParam2 = 'SF: ';
                        titleList2 = fValsUnique;
                    case 5
                        oList(row,:) = 1:length(oValsUnique);
                        titleParam2 = 'Ori: ';
                        titleList2 = oValsUnique;
                    case 6
                        cList(row,:) = 1:length(cValsUnique);
                        titleParam2 = 'Contr: ';
                        titleList2 = cValsUnique;
                    case 7
                        tList(row,:) = 1:length(tValsUnique);
                        titleParam2 = 'TF: ';
                        titleList2 = tValsUnique;
                    otherwise
                        disp('No particular gabor or parameter selected');
                end
            end


    % [Vinay] - Get the lengths of indices in parameterCombinations

    aLen = length(aValsUnique);
    eLen = length(eValsUnique);
    sLen = length(sValsUnique);
    fLen = length(fValsUnique);
    oLen = length(oValsUnique);
    cLen = length(cValsUnique);
    tLen = length(tValsUnique);

    % If more than one value, then length is one greater for all the values
    % together
    if (aLen> 1)           aLen=aLen+1;                    end
    if (eLen> 1)           eLen=eLen+1;                    end
    if (sLen> 1)           sLen=sLen+1;                    end
    if (fLen> 1)           fLen=fLen+1;                    end
    if (oLen> 1)           oLen=oLen+1;                    end
    if (cLen> 1)           cLen=cLen+1;                    end
    if (tLen> 1)           tLen=tLen+1;                    end


elseif strncmp(protocolName,'CRS',3)

        % [Vinay] - repeat the set of parameters for each gabor depending on the
        % number of rows and columns to be drawn
        aList = repmat(a,numRows,numCols);
        eList = repmat(e,numRows,numCols);
        sList = repmat(s,numRows,numCols);
        fList = repmat(f,numRows,numCols);
        oList = repmat(o,numRows,numCols);
        cList = repmat(c,numRows,numCols);
        tList = repmat(t,numRows,numCols);
        rList = repmat(r,numRows,numCols);
        pList = repmat(p,numRows,numCols);

        % [Vinay] decide the row parameter based on paramNum1
        for row = 1:numRows
            switch (paramNum1-1)
                case 1
                    aList(gaborNum1+((row-1)*3),:) = row;                    
                    % along the rows, we get three consecutive entries for the 3 
                    % gabors. Therefore to go to the next index of the same gabor,
                    % we go ahead by 3 indices
                    % For gabor1 every new row takes a new value of param1.
                    % So assign the 'row' value in the list above
                    titleParam1 = 'Azi: ';
                    titleList1 = aValsUnique{gaborNum1};
                case 2
                    eList(gaborNum1+((row-1)*3),:) = row; 
                    titleParam1 = 'Ele: ';
                    titleList1 = eValsUnique{gaborNum1};
                case 3
                    sList(gaborNum1+((row-1)*3),:) = row; 
                    titleParam1 = 'Sigma: ';
                    titleList1 = sValsUnique{gaborNum1};
                case 4
                    fList(gaborNum1+((row-1)*3),:) = row; 
                    titleParam1 = 'SF: ';
                    titleList1 = fValsUnique{gaborNum1};
                case 5
                    oList(gaborNum1+((row-1)*3),:) = row; 
                    titleParam1 = 'Ori: ';
                    titleList1 = oValsUnique{gaborNum1};
                case 6
                    cList(gaborNum1+((row-1)*3),:) = row; 
                    titleParam1 = 'Contr: ';
                    titleList1 = cValsUnique{gaborNum1};
                case 7
                    tList(gaborNum1+((row-1)*3),:) = row; 
                    titleParam1 = 'TF: ';
                    titleList1 = tValsUnique{gaborNum1};
                case 8
                    rList(gaborNum1+((row-1)*3),:) = row; 
                    titleParam1 = 'Rad: ';
                    titleList1 = rValsUnique{gaborNum1};
                case 9
                    pList(gaborNum1+((row-1)*3),:) = row; 
                    titleParam1 = 'SPhs: ';
                    titleList1 = pValsUnique{gaborNum1};
                otherwise
                    disp('No particular gabor or parameter selected');
            end
        end


        % [Vinay] decide the column parameter based on paramNum2
        for row = 1:numRows
            switch (paramNum2-1)
                case 1
                    aList(gaborNum2+((row-1)*3),:) = 1:length(aValsUnique{gaborNum2});
                    % along the rows, we get three consecutive entries for the 3 
                    % gabors. Therefore to go to the next index of the same gabor,
                    % we go ahead by 3 indices
                    % Besides for every row we now have to assign indices
                    % incrementing along the columns.
                    % Basically the column entries go from 1 to 
                    % length{nValsUnique}
                    titleParam2 = 'Azi: ';
                    titleList2 = aValsUnique{gaborNum2};
                case 2
                    eList(gaborNum2+((row-1)*3),:) = 1:length(eValsUnique{gaborNum2}); 
                    titleParam2 = 'Ele: ';
                    titleList2 = eValsUnique{gaborNum2};
                case 3
                    sList(gaborNum2+((row-1)*3),:) = 1:length(sValsUnique{gaborNum2});
                    titleParam2 = 'Sigma: ';
                    titleList2 = sValsUnique{gaborNum2};
                case 4
                    fList(gaborNum2+((row-1)*3),:) = 1:length(fValsUnique{gaborNum2});
                    titleParam2 = 'SF: ';
                    titleList2 = fValsUnique{gaborNum2};
                case 5
                    oList(gaborNum2+((row-1)*3),:) = 1:length(oValsUnique{gaborNum2});
                    titleParam2 = 'Ori: ';
                    titleList2 = oValsUnique{gaborNum2};
                case 6
                    cList(gaborNum2+((row-1)*3),:) = 1:length(cValsUnique{gaborNum2});
                    titleParam2 = 'Contr: ';
                    titleList2 = cValsUnique{gaborNum2};
                case 7
                    tList(gaborNum2+((row-1)*3),:) = 1:length(tValsUnique{gaborNum2});
                    titleParam2 = 'TF: ';
                    titleList2 = tValsUnique{gaborNum2};
                case 8
                    rList(gaborNum2+((row-1)*3),:) = 1:length(rValsUnique{gaborNum2});
                    titleParam2 = 'Rad: ';
                    titleList2 = rValsUnique{gaborNum2};
                case 9
                    pList(gaborNum2+((row-1)*3),:) = 1:length(pValsUnique{gaborNum2});
                    titleParam2 = 'SPhs: ';
                    titleList2 = pValsUnique{gaborNum2};
                otherwise
                    disp('No particular gabor or parameter selected');
            end
        end


        % [Vinay] - Get the lengths of indices in parameterCombinations
        for i = 1:3
            aLen(i) = length(aValsUnique{i});
            eLen(i) = length(eValsUnique{i});
            sLen(i) = length(sValsUnique{i});
            fLen(i) = length(fValsUnique{i});
            oLen(i) = length(oValsUnique{i});
            cLen(i) = length(cValsUnique{i});
            tLen(i) = length(tValsUnique{i});
            pLen(i) = length(pValsUnique{i}); % Vinay - for CRS
            rLen(i) = length(rValsUnique{i}); % Vinay - for CRS

            % If more than one value, then length is one greater for all the values
            % together
            if (aLen(i)> 1)           aLen(i)=aLen(i)+1;                    end
            if (eLen(i)> 1)           eLen(i)=eLen(i)+1;                    end
            if (sLen(i)> 1)           sLen(i)=sLen(i)+1;                    end
            if (fLen(i)> 1)           fLen(i)=fLen(i)+1;                    end
            if (oLen(i)> 1)           oLen(i)=oLen(i)+1;                    end
            if (cLen(i)> 1)           cLen(i)=cLen(i)+1;                    end
            if (tLen(i)> 1)           tLen(i)=tLen(i)+1;                    end
            if (pLen(i)> 1)           pLen(i)=pLen(i)+1;                    end % Vinay - for CRS
            if (rLen(i)> 1)           rLen(i)=rLen(i)+1;                    end % Vinay - for CRS
        end
        
end

        % Main loop
%         computationVals=zeros(1,numCols);
        for k=1:numRows
            for j = 1:numCols
                
                if strncmp(protocolName,'GRF',3)
                    clear goodPos
                    goodPos = parameterCombinations{aList(k,j),eList(k,j),sList(k,j),fList(k,j),oList(k,j),cList(k,j),tList(k,j)};
            %         goodPos = setdiff(goodPos,badTrials);
            
                    tagPos1 = [num2str(aList(k,j)) num2str(eList(k,j)) num2str(sList(k,j))...
                               num2str(fList(k,j)) num2str(oList(k,j)) num2str(cList(k,j))...
                               num2str(tList(k,j))];
                    
                    %------------------
                    % Vinay - select good trials as per the electrode(s)
            %         useAllBadTrials = 1;
                    if useAllBadTrials && existsBadTrialFile
                        elecIndex1 = (strcmp(channelString,nameElec) == 1);

                        if ~useBipolar
                            elecBadTrials = allBadTrials{elecIndex1};
                            disp(['No. of Bad trials for ' channelString ': ' num2str(length(elecBadTrials))]);
                        else
                            elecIndex2 = (strcmp(analogChannelString2,nameElec) == 1);
                            elecBadTrials = unique([allBadTrials{elecIndex1} allBadTrials{elecIndex2}]);
                            disp(['No. of Bad trials for ' channelString ' and ' analogChannelString2 ': ' num2str(length(elecBadTrials))]);
                        end
                        goodPos = setdiff(goodPos,elecBadTrials);
                    else
                        goodPos = setdiff(goodPos,badTrials);    
                    end
                    %-------------------

                elseif strncmp(protocolName,'CRS',3)
                
                clear goodPos
            %     goodPos = parameterCombinations{aList(j),eList(j),sList(j),fList(j),oList(j),cList(j),tList(j)};
            %     goodPos = setdiff(goodPos,badTrials);

                % [Vinay] - goodPos will be obtained by taking the goodPos for all
                % the three gabors and finding the common trials in them
                clear pos1 pos2 pos3
                switch protocolNumber
                           case 1 % Ring protocol
                               pos1 = parameterCombinations{aList(1+(k-1)*3,j),eList(1+(k-1)*3,j),sList(1+(k-1)*3,j),fList(1+(k-1)*3,j),oList(1+(k-1)*3,j),cList(1+(k-1)*3,j),tList(1+(k-1)*3,j),pList(1+(k-1)*3,j),rList(1+(k-1)*3,j),1}; 
                               % S params define C & S
                               % S i.e. kGabor0 corresponds to row number 1,4,7,..
                               % and so on. Hence it is 1+(k-1)*3
                               pos2 = parameterCombinations{aLen(2),eLen(2),sLen(2),fLen(2),oLen(2),cLen(2),tLen(2),pLen(2),rList(2+(k-1)*3,j),2}; % for R - only rad is imp
                               % R i.e. kGabor1 corresponds to row number 2,5,8,..
                               % and so on. Hence it is 2+(k-1)*3

                               % [Vinay] -  C is matched to S and therefore the entries
                               % for gabor3 i.e. can be ignored. Otherwise this can 
                               % give spurious results if the param values for C are 
                               % not set equal to those for S. Hence pos3 is the full
                               % set here
                               pos3 = parameterCombinations{aLen(3),eLen(3),sLen(3),fLen(3),oLen(3),cLen(3),tLen(3),pLen(3),rList(3+(k-1)*3,j),3}; % for C
                               % also keeping r variable here, just in case you
                               % have multiple r values
                               
                               % keep the position tags for storing MP data
                               tagPos1 = [num2str(aList(1+(k-1)*3,j)) num2str(eList(1+(k-1)*3,j)) num2str(sList(1+(k-1)*3,j))...
                                   num2str(fList(1+(k-1)*3,j)) num2str(oList(1+(k-1)*3,j)) num2str(cList(1+(k-1)*3,j))...
                                   num2str(tList(1+(k-1)*3,j)) num2str(pList(1+(k-1)*3,j)) num2str(rList(1+(k-1)*3,j))];
                               tagPos2 = [num2str(aLen(2)) num2str(eLen(2)) num2str(sLen(2)) num2str(fLen(2))...
                                   num2str(oLen(2)) num2str(cLen(2)) num2str(tLen(2)) num2str(pLen(2)) num2str(rList(2+(k-1)*3,j))];
                               tagPos3 = [num2str(aLen(3)) num2str(eLen(3)) num2str(sLen(3)) num2str(fLen(3))...
                                   num2str(oLen(3)) num2str(cLen(3)) num2str(tLen(3)) num2str(pLen(3)) num2str(rList(3+(k-1)*3,j))]; % for C

                           case 2 % Contrast Ring Protocol
                               pos1 = parameterCombinations{aList(1+(k-1)*3,j),eList(1+(k-1)*3,j),sList(1+(k-1)*3,j),fList(1+(k-1)*3,j),oList(1+(k-1)*3,j),cList(1+(k-1)*3,j),tList(1+(k-1)*3,j),pList(1+(k-1)*3,j),rList(1+(k-1)*3,j),1};
                               % S params define C & S
                               % S i.e. kGabor0 corresponds to row number 1,4,7,..
                               % and so on. Hence it is 1+(k-1)*3
                               pos2 = parameterCombinations{aLen(2),eLen(2),sLen(2),fLen(2),oLen(2),cList(2+(k-1)*3,j),tLen(2),pLen(2),rList(2+(k-1)*3,j),2}; % for R - only rad & cont
                               % R i.e. kGabor1 corresponds to row number 2,5,8,..
                               % and so on. Hence it is 2+(k-1)*3

                               % [Vinay] -  C is matched to S and therefore the entries
                               % for gabor3 i.e. can be ignored. Otherwise this can 
                               % give spurious results if the param values for C are 
                               % not set equal to those for S. Hence pos3 is the full
                               % set here
                               pos3 = parameterCombinations{aLen(3),eLen(3),sLen(3),fLen(3),oLen(3),cLen(3),tLen(3),pLen(3),rList(3+(k-1)*3,j),3}; % for C
                               % also keeping r variable here, just in case you
                               % have multiple r values
                               
                               tagPos1 = [num2str(aList(1+(k-1)*3,j)) num2str(eList(1+(k-1)*3,j)) num2str(sList(1+(k-1)*3,j))...
                                   num2str(fList(1+(k-1)*3,j)) num2str(oList(1+(k-1)*3,j)) num2str(cList(1+(k-1)*3,j))...
                                   num2str(tList(1+(k-1)*3,j)) num2str(pList(1+(k-1)*3,j)) num2str(rList(1+(k-1)*3,j))];
                               tagPos2 = [num2str(aLen(2)) num2str(eLen(2)) num2str(sLen(2)) num2str(fLen(2))...
                                   num2str(oLen(2)) num2str(cList(2+(k-1)*3,j)) num2str(tLen(2)) num2str(pLen(2)) num2str(rList(2+(k-1)*3,j))];
                               tagPos3 = [num2str(aLen(3)) num2str(eLen(3)) num2str(sLen(3)) num2str(fLen(3))...
                                   num2str(oLen(3)) num2str(cLen(3)) num2str(tLen(3)) num2str(pLen(3)) num2str(rList(3+(k-1)*3,j))];

                           case 3 % Dual Contrast Protocol
                               pos1 = parameterCombinations{aLen(1),eLen(1),sLen(1),fLen(1),oLen(1),cLen(1),tLen(1),pLen(1),rLen(1),1}; % for C
                               % [Vinay] - S is hidden in this case and only C & R
                               % define the stimuli. pos1 is therefore the full set
                               % here
                               pos2 = parameterCombinations{aLen(2),eLen(2),sLen(2),fLen(2),oLen(2),cList(2+(k-1)*3,j),tLen(2),pLen(2),rList(2+(k-1)*3,j),2}; % for R - only rad & cont
                               % [Vinay] -  R params - cont & rad. Rest are matched to
                               % C's params
                               pos3 = parameterCombinations{aList(3+(k-1)*3,j),eList(3+(k-1)*3,j),sList(3+(k-1)*3,j),fList(3+(k-1)*3,j),oList(3+(k-1)*3,j),cList(3+(k-1)*3,j),tList(3+(k-1)*3,j),pList(3+(k-1)*3,j),rList(3+(k-1)*3,j),3};
                               % for C. C is the defining gabor here
                               
                               tagPos1 = [num2str(aLen(1)) num2str(eLen(1)) num2str(sLen(1))...
                                   num2str(fLen(1)) num2str(oLen(1)) num2str(cLen(1))...
                                   num2str(tLen(1)) num2str(pLen(1)) num2str(rLen(1))];
                               tagPos2 = [num2str(aLen(2)) num2str(eLen(2)) num2str(sLen(2)) num2str(fLen(2))...
                                   num2str(oLen(2)) num2str(cList(2+(k-1)*3,j)) num2str(tLen(2)) num2str(pLen(2)) num2str(rList(2+(k-1)*3,j))];
                               tagPos3 = [num2str(aList(3+(k-1)*3,j)) num2str(eList(3+(k-1)*3,j)) num2str(sList(3+(k-1)*3,j))...
                                   num2str(fList(3+(k-1)*3,j)) num2str(oList(3+(k-1)*3,j)) num2str(cList(3+(k-1)*3,j))...
                                   num2str(tList(3+(k-1)*3,j)) num2str(pList(3+(k-1)*3,j)) num2str(rList(3+(k-1)*3,j))];
                               

                           case 4 % Dual Orientation Protocol
                               pos1 = parameterCombinations{aLen(1),eLen(1),sLen(1),fLen(1),oLen(1),cLen(1),tLen(1),pLen(1),rLen(1),1}; % for C
                               % [Vinay] - S is hidden in this case and only C & R
                               % define the stimuli. pos1 is therefore the full set
                               % here
                               pos2 = parameterCombinations{aLen(2),eLen(2),sLen(2),fLen(2),oList(2+(k-1)*3,j),cLen(2),tLen(2),pLen(2),rList(2+(k-1)*3,j),2}; % for R - only rad & ori
                               % [Vinay] -  R params - ori & rad. Rest are matched to
                               % C's params
                               pos3 = parameterCombinations{aList(3+(k-1)*3,j),eList(3+(k-1)*3,j),sList(3+(k-1)*3,j),fList(3+(k-1)*3,j),oList(3+(k-1)*3,j),cList(3+(k-1)*3,j),tList(3+(k-1)*3,j),pList(3+(k-1)*3,j),rList(3+(k-1)*3,j),3};
                               % for C. C is the defining gabor here
                               
                               tagPos1 = [num2str(aLen(1)) num2str(eLen(1)) num2str(sLen(1))...
                                   num2str(fLen(1)) num2str(oLen(1)) num2str(cLen(1))...
                                   num2str(tLen(1)) num2str(pLen(1)) num2str(rLen(1))];
                               tagPos2 = [num2str(aLen(2)) num2str(eLen(2)) num2str(sLen(2)) num2str(fLen(2))...
                                   num2str(oList(2+(k-1)*3,j)) num2str(cLen(2)) num2str(tLen(2)) num2str(pLen(2)) num2str(rList(2+(k-1)*3,j))];
                               tagPos3 = [num2str(aList(3+(k-1)*3,j)) num2str(eList(3+(k-1)*3,j)) num2str(sList(3+(k-1)*3,j))...
                                   num2str(fList(3+(k-1)*3,j)) num2str(oList(3+(k-1)*3,j)) num2str(cList(3+(k-1)*3,j))...
                                   num2str(tList(3+(k-1)*3,j)) num2str(pList(3+(k-1)*3,j)) num2str(rList(3+(k-1)*3,j))];

                           case 5 % Dual Phase Protocol
                               pos1 = parameterCombinations{aLen(1),eLen(1),sLen(1),fLen(1),oLen(1),cLen(1),tLen(1),pLen(1),rLen(1),1}; % for C
                               % [Vinay] - S is hidden in this case and only C & R
                               % define the stimuli. pos1 is therefore the full set
                               % here
                               pos2 = parameterCombinations{aLen(2),eLen(2),sLen(2),fLen(2),oLen(2),cLen(2),tLen(2),pList(2+(k-1)*3,j),rList(2+(k-1)*3,j),2}; % for R - only rad & p
                               % [Vinay] -  R params - p & rad. Rest are matched to
                               % C's params
                               pos3 = parameterCombinations{aList(3+(k-1)*3,j),eList(3+(k-1)*3,j),sList(3+(k-1)*3,j),fList(3+(k-1)*3,j),oList(3+(k-1)*3,j),cList(3+(k-1)*3,j),tList(3+(k-1)*3,j),pList(3+(k-1)*3,j),rList(3+(k-1)*3,j),3};
                               % for C. C is the defining gabor here
                               
                               tagPos1 = [num2str(aLen(1)) num2str(eLen(1)) num2str(sLen(1))...
                                   num2str(fLen(1)) num2str(oLen(1)) num2str(cLen(1))...
                                   num2str(tLen(1)) num2str(pLen(1)) num2str(rLen(1))];
                               tagPos2 = [num2str(aLen(2)) num2str(eLen(2)) num2str(sLen(2)) num2str(fLen(2))...
                                   num2str(oLen(2)) num2str(cLen(2)) num2str(tLen(2)) num2str(pList(2+(k-1)*3,j)) num2str(rList(2+(k-1)*3,j))];
                               tagPos3 = [num2str(aList(3+(k-1)*3,j)) num2str(eList(3+(k-1)*3,j)) num2str(sList(3+(k-1)*3,j))...
                                   num2str(fList(3+(k-1)*3,j)) num2str(oList(3+(k-1)*3,j)) num2str(cList(3+(k-1)*3,j))...
                                   num2str(tList(3+(k-1)*3,j)) num2str(pList(3+(k-1)*3,j)) num2str(rList(3+(k-1)*3,j))];

                           case 6 %Orientation Ring Protocol
                               pos1 = parameterCombinations{aList(1+(k-1)*3,j),eList(1+(k-1)*3,j),sList(1+(k-1)*3,j),fList(1+(k-1)*3,j),oList(1+(k-1)*3,j),cList(1+(k-1)*3,j),tList(1+(k-1)*3,j),pList(1+(k-1)*3,j),rList(1+(k-1)*3,j),1};
                               % S params define C & S
                               % S i.e. kGabor0 corresponds to row number 1,4,7,..
                               % and so on. Hence it is 1+(k-1)*3
                               pos2 = parameterCombinations{aLen(2),eLen(2),sLen(2),fLen(2),oList(2+(k-1)*3,j),cLen(2),tLen(2),pLen(2),rList(2+(k-1)*3,j),2}; % for R - only rad & ori
                               % [Vinay] -  C is matched to S and therefore the entries
                               % for gabor3 i.e. can be ignored. Otherwise this can 
                               % give spurious results if the param values for C are 
                               % not set equal to those for S. Hence pos3 is the full
                               % set here
                               pos3 = parameterCombinations{aLen(3),eLen(3),sLen(3),fLen(3),oLen(3),cLen(3),tLen(3),pLen(3),rList(3+(k-1)*3,j),3}; % for C
                               % also keeping r variable here, just in case you
                               % have multiple r values
                               
                               tagPos1 = [num2str(aList(1+(k-1)*3,j)) num2str(eList(1+(k-1)*3,j)) num2str(sList(1+(k-1)*3,j))...
                                   num2str(fList(1+(k-1)*3,j)) num2str(oList(1+(k-1)*3,j)) num2str(cList(1+(k-1)*3,j))...
                                   num2str(tList(1+(k-1)*3,j)) num2str(pList(1+(k-1)*3,j)) num2str(rList(1+(k-1)*3,j))];
                               tagPos2 = [num2str(aLen(2)) num2str(eLen(2)) num2str(sLen(2)) num2str(fLen(2))...
                                   num2str(oList(2+(k-1)*3,j)) num2str(cLen(2)) num2str(tLen(2)) num2str(pLen(2)) num2str(rList(2+(k-1)*3,j))];
                               tagPos3 = [num2str(aLen(3)) num2str(eLen(3)) num2str(sLen(3)) num2str(fLen(3))...
                                   num2str(oLen(3)) num2str(cLen(3)) num2str(tLen(3)) num2str(pLen(3)) num2str(rList(3+(k-1)*3,j))];

                           case 7 %Phase Ring Protocol
                               pos1 = parameterCombinations{aList(1+(k-1)*3,j),eList(1+(k-1)*3,j),sList(1+(k-1)*3,j),fList(1+(k-1)*3,j),oList(1+(k-1)*3,j),cList(1+(k-1)*3,j),tList(1+(k-1)*3,j),pList(1+(k-1)*3,j),rList(1+(k-1)*3,j),1};
                               % S params define C & S
                               % S i.e. kGabor0 corresponds to row number 1,4,7,..
                               % and so on. Hence it is 1+(k-1)*3
                               pos2 = parameterCombinations{aLen(2),eLen(2),sLen(2),fLen(2),oLen(2),cLen(2),tLen(2),pList(2+(k-1)*3,j),rList(2+(k-1)*3,j),2}; % for R - only rad & p
                               % [Vinay] -  C is matched to S and therefore the entries
                               % for gabor3 i.e. can be ignored. Otherwise this can 
                               % give spurious results if the param values for C are 
                               % not set equal to those for S. Hence pos3 is the full
                               % set here
                               pos3 = parameterCombinations{aLen(3),eLen(3),sLen(3),fLen(3),oLen(3),cLen(3),tLen(3),pLen(3),rList(3+(k-1)*3,j),3}; % for C
                               % also keeping r variable here, just in case you
                               % have multiple r values
                               
                               tagPos1 = [num2str(aList(1+(k-1)*3,j)) num2str(eList(1+(k-1)*3,j)) num2str(sList(1+(k-1)*3,j))...
                                   num2str(fList(1+(k-1)*3,j)) num2str(oList(1+(k-1)*3,j)) num2str(cList(1+(k-1)*3,j))...
                                   num2str(tList(1+(k-1)*3,j)) num2str(pList(1+(k-1)*3,j)) num2str(rList(1+(k-1)*3,j))];
                               tagPos2 = [num2str(aLen(2)) num2str(eLen(2)) num2str(sLen(2)) num2str(fLen(2))...
                                   num2str(oLen(2)) num2str(cLen(2)) num2str(tLen(2)) num2str(pList(2+(k-1)*3,j)) num2str(rList(2+(k-1)*3,j))];
                               tagPos3 = [num2str(aLen(3)) num2str(eLen(3)) num2str(sLen(3)) num2str(fLen(3))...
                                   num2str(oLen(3)) num2str(cLen(3)) num2str(tLen(3)) num2str(pLen(3)) num2str(rList(3+(k-1)*3,j))];

                           case 8 %Drifting Ring Protocol
                               pos1 = parameterCombinations{aList(1+(k-1)*3,j),eList(1+(k-1)*3,j),sList(1+(k-1)*3,j),fList(1+(k-1)*3,j),oList(1+(k-1)*3,j),cList(1+(k-1)*3,j),tList(1+(k-1)*3,j),pList(1+(k-1)*3,j),rList(1+(k-1)*3,j),1};
                               % S params define C & S
                               % S i.e. kGabor0 corresponds to row number 1,4,7,..
                               % and so on. Hence it is 1+(k-1)*3
                               pos2 = parameterCombinations{aLen(2),eLen(2),sLen(2),fLen(2),oLen(2),cLen(2),tList(2+(k-1)*3,j),pLen(2),rList(2+(k-1)*3,j),2}; % for R - only rad & t
                               % [Vinay] -  C is matched to S and therefore the entries
                               % for gabor3 i.e. can be ignored. Otherwise this can 
                               % give spurious results if the param values for C are 
                               % not set equal to those for S. Hence pos3 is the full
                               % set here
                               pos3 = parameterCombinations{aLen(3),eLen(3),sLen(3),fLen(3),oLen(3),cLen(3),tLen(3),pLen(3),rList(3+(k-1)*3,j),3}; % for C
                               % also keeping r variable here, just in case you
                               % have multiple r values
                               
                               tagPos1 = [num2str(aList(1+(k-1)*3,j)) num2str(eList(1+(k-1)*3,j)) num2str(sList(1+(k-1)*3,j))...
                                   num2str(fList(1+(k-1)*3,j)) num2str(oList(1+(k-1)*3,j)) num2str(cList(1+(k-1)*3,j))...
                                   num2str(tList(1+(k-1)*3,j)) num2str(pList(1+(k-1)*3,j)) num2str(rList(1+(k-1)*3,j))];
                               tagPos2 = [num2str(aLen(2)) num2str(eLen(2)) num2str(sLen(2)) num2str(fLen(2))...
                                   num2str(oLen(2)) num2str(cLen(2)) num2str(tList(2+(k-1)*3,j)) num2str(pLen(2)) num2str(rList(2+(k-1)*3,j))];
                               tagPos3 = [num2str(aLen(3)) num2str(eLen(3)) num2str(sLen(3)) num2str(fLen(3))...
                                   num2str(oLen(3)) num2str(cLen(3)) num2str(tLen(3)) num2str(pLen(3)) num2str(rList(3+(k-1)*3,j))];

                           case 9 %Cross Orientation Protocol
                               pos1 = parameterCombinations{aLen(1),eLen(1),sLen(1),fLen(1),oLen(1),cLen(1),tLen(1),pLen(1),rLen(1),1}; % for C
                               % [Vinay] - S is hidden in this case and only C & R
                               % define the stimuli. pos1 is therefore the full set
                               % here
                               pos2 = parameterCombinations{aList(2+(k-1)*3,j),eList(2+(k-1)*3,j),sList(2+(k-1)*3,j),fList(2+(k-1)*3,j),oList(2+(k-1)*3,j),cList(2+(k-1)*3,j),tList(2+(k-1)*3,j),pList(2+(k-1)*3,j),rList(2+(k-1)*3,j),2};
                               % for R. R is one of the defining gabor here
                               pos3 = parameterCombinations{aList(3+(k-1)*3,j),eList(3+(k-1)*3,j),sList(3+(k-1)*3,j),fList(3+(k-1)*3,j),oList(3+(k-1)*3,j),cList(3+(k-1)*3,j),tList(3+(k-1)*3,j),pList(3+(k-1)*3,j),rList(3+(k-1)*3,j),3};
                               % for C. C is one of the defining gabor here
                               
                               tagPos1 = [num2str(aLen(1)) num2str(eLen(1)) num2str(sLen(1))...
                                   num2str(fLen(1)) num2str(oLen(1)) num2str(cLen(1))...
                                   num2str(tLen(1)) num2str(pLen(1)) num2str(rLen(1))];
                               tagPos2 = [num2str(aList(2+(k-1)*3,j)) num2str(eList(2+(k-1)*3,j)) num2str(sList(2+(k-1)*3,j))...
                                   num2str(fList(2+(k-1)*3,j)) num2str(oList(2+(k-1)*3,j)) num2str(cList(2+(k-1)*3,j))...
                                   num2str(tList(2+(k-1)*3,j)) num2str(pList(2+(k-1)*3,j)) num2str(rList(2+(k-1)*3,j))];
                               tagPos3 = [num2str(aList(3+(k-1)*3,j)) num2str(eList(3+(k-1)*3,j)) num2str(sList(3+(k-1)*3,j))...
                                   num2str(fList(3+(k-1)*3,j)) num2str(oList(3+(k-1)*3,j)) num2str(cList(3+(k-1)*3,j))...
                                   num2str(tList(3+(k-1)*3,j)) num2str(pList(3+(k-1)*3,j)) num2str(rList(3+(k-1)*3,j))];

                           case 10 %Annulus Fixed Protocol
                               pos1 = parameterCombinations{aList(1+(k-1)*3,j),eList(1+(k-1)*3,j),sList(1+(k-1)*3,j),fList(1+(k-1)*3,j),oList(1+(k-1)*3,j),cList(1+(k-1)*3,j),tList(1+(k-1)*3,j),pList(1+(k-1)*3,j),rList(1+(k-1)*3,j),1};
                               % S params define C & S
                               % S i.e. kGabor0 corresponds to row number 1,4,7,..
                               % and so on. Hence it is 1+(k-1)*3
                               pos2 = parameterCombinations{aLen(2),eLen(2),sLen(2),fLen(2),oLen(2),cList(2+(k-1)*3,j),tLen(2),pLen(2),rList(2+(k-1)*3,j),2}; % for R - only rad & cont
                               % [Vinay] -  C is matched to S and therefore the entries
                               % for gabor3 i.e. can be ignored. Otherwise this can 
                               % give spurious results if the param values for C are 
                               % not set equal to those for S. However, the rad of C is
                               % still a relevant parameter and has to be taken care of
                               pos3 = parameterCombinations{aLen(3),eLen(3),sLen(3),fLen(3),oLen(3),cLen(3),tLen(3),pLen(3),rList(3+(k-1)*3,j),3}; % for C - only rad
                               
                               tagPos1 = [num2str(aList(1+(k-1)*3,j)) num2str(eList(1+(k-1)*3,j)) num2str(sList(1+(k-1)*3,j))...
                                   num2str(fList(1+(k-1)*3,j)) num2str(oList(1+(k-1)*3,j)) num2str(cList(1+(k-1)*3,j))...
                                   num2str(tList(1+(k-1)*3,j)) num2str(pList(1+(k-1)*3,j)) num2str(rList(1+(k-1)*3,j))];
                               tagPos2 = [num2str(aLen(2)) num2str(eLen(2)) num2str(sLen(2)) num2str(fLen(2))...
                                   num2str(oLen(2)) num2str(cList(2+(k-1)*3,j)) num2str(tLen(2)) num2str(pLen(2)) num2str(rList(2+(k-1)*3,j))];
                               tagPos3 = [num2str(aLen(3)) num2str(eLen(3)) num2str(sLen(3)) num2str(fLen(3))...
                                   num2str(oLen(3)) num2str(cLen(3)) num2str(tLen(3)) num2str(pLen(3)) num2str(rList(3+(k-1)*3,j))];

                    otherwise
                        pos1 = parameterCombinations{aList(1+(k-1)*3,j),eList(1+(k-1)*3,j),sList(1+(k-1)*3,j),fList(1+(k-1)*3,j),oList(1+(k-1)*3,j),cList(1+(k-1)*3,j),tList(1+(k-1)*3,j),pList(1+(k-1)*3,j),rList(1+(k-1)*3,j),1};
                        % good positions for gabor 1 i.e. S
                        pos2 = parameterCombinations{aList(2+(k-1)*3,j),eList(2+(k-1)*3,j),sList(2+(k-1)*3,j),fList(2+(k-1)*3,j),oList(2+(k-1)*3,j),cList(2+(k-1)*3,j),tList(2+(k-1)*3,j),pList(2+(k-1)*3,j),rList(2+(k-1)*3,j),2};
                        pos3 = parameterCombinations{aList(3+(k-1)*3,j),eList(3+(k-1)*3,j),sList(3+(k-1)*3,j),fList(3+(k-1)*3,j),oList(3+(k-1)*3,j),cList(3+(k-1)*3,j),tList(3+(k-1)*3,j),pList(3+(k-1)*3,j),rList(3+(k-1)*3,j),3};
                        
                           tagPos1 = [num2str(aList(1+(k-1)*3,j)) num2str(eList(1+(k-1)*3,j)) num2str(sList(1+(k-1)*3,j))...
                               num2str(fList(1+(k-1)*3,j)) num2str(oList(1+(k-1)*3,j)) num2str(cList(1+(k-1)*3,j))...
                               num2str(tList(1+(k-1)*3,j)) num2str(pList(1+(k-1)*3,j)) num2str(rList(1+(k-1)*3,j))];
                           tagPos2 = [num2str(aList(2+(k-1)*3,j)) num2str(eList(2+(k-1)*3,j)) num2str(sList(2+(k-1)*3,j))...
                               num2str(fList(2+(k-1)*3,j)) num2str(oList(2+(k-1)*3,j)) num2str(cList(2+(k-1)*3,j))...
                               num2str(tList(2+(k-1)*3,j)) num2str(pList(2+(k-1)*3,j)) num2str(rList(2+(k-1)*3,j))];
                           tagPos3 = [num2str(aList(3+(k-1)*3,j)) num2str(eList(3+(k-1)*3,j)) num2str(sList(3+(k-1)*3,j))...
                               num2str(fList(3+(k-1)*3,j)) num2str(oList(3+(k-1)*3,j)) num2str(cList(3+(k-1)*3,j))...
                               num2str(tList(3+(k-1)*3,j)) num2str(pList(3+(k-1)*3,j)) num2str(rList(3+(k-1)*3,j))];
                end


            goodPos = intersect(pos1,pos2);
            goodPos = intersect(goodPos,pos3);
            
            %------------------
            % Vinay - select good trials as per the electrode(s)
    %         useAllBadTrials = 1;
            if useAllBadTrials && existsBadTrialFile
                elecIndex1 = (strcmp(channelString,nameElec) == 1);

                if ~useBipolar
                    elecBadTrials = allBadTrials{elecIndex1};
                    disp(['No. of Bad trials for ' channelString ': ' num2str(length(elecBadTrials))]);
                else
                    elecIndex2 = (strcmp(analogChannelString2,nameElec) == 1);
                    elecBadTrials = unique([allBadTrials{elecIndex1} allBadTrials{elecIndex2}]);
                    disp(['No. of Bad trials for ' channelString ' and ' analogChannelString2 ': ' num2str(length(elecBadTrials))]);
                end
                goodPos = setdiff(goodPos,elecBadTrials);
            else
                goodPos = setdiff(goodPos,badTrials);
            end
            %-------------------
            
        end

            if isempty(goodPos)
                disp('No entries for this combination..')
            else
                disp(['pos=' num2str(k) ',' num2str(j) ',n=' num2str(length(goodPos))]);

                Fs = round(1/(timeVals(2)-timeVals(1)));
                BLRange = uint16((BLMax-BLMin)*Fs);
                STRange = uint16((STMax-STMin)*Fs);
                BLPos = find(timeVals>=BLMin,1)+ (1:BLRange);
                STPos = find(timeVals>=STMin,1)+ (1:STRange);

                xsBL = 0:1/(BLMax-BLMin):Fs-1/(BLMax-BLMin);
                xsST = 0:1/(STMax-STMin):Fs-1/(STMax-STMin);


                % Vinay - added this notch data check
                if notchData
                    analogData = analogDataNotched;
                end

                if useBipolar
                    analogChannelString1 = channelString; % first electrode selected
                    clear analogData analogDataNotched
                    load(fullfile(folderLFP,analogChannelString1));
                
                    if notchData
                        analogData1 = analogDataNotched;
                    else
                        analogData1 = analogData; % Vinay - these are the 
                        % freshly loaded values for electrode 1
                    end
                    % Vinay - store the analogData 
                    % for electrode 1 in a separate variable 
                    % otherwise loading the data for 2nd electrode will
                    % overwrite analogData
                    clear analogData analogDataNotched
                    load(fullfile(folderLFP,analogChannelString2));
                    if notchData
                        analogData2 = analogDataNotched;
                    else
                        analogData2 = analogData; % Vinay - these are the 
                        % freshly loaded values for electrode 2
                    end
                    analogData = analogData1 - analogData2;
                end
                
                
                %--------------------------------------------------------------
                % Implement the method set by tfMethod varibale

                if (tfMethod == 1) % MTM
                    %----------------------------------------------------------
                    if (plotStyle == 3) % line plot (power vs f)
                        if (spectrumType == 1) % raw
                            [SBL,f]=mtspectrumc((analogData(goodPos,BLPos))',mtmParams); % baseline period
                            [SST,f]=mtspectrumc((analogData(goodPos,STPos))',mtmParams); % stimulus period
                            plot(tfplotHandles(k,j),f,conv2Log(SST),'color',plotColor,'Linewidth',plotLineWidth);
                            set(tfplotHandles(k,j),'Nextplot','add');
                            plot(tfplotHandles(k,j),f,conv2Log(SBL),'color','g','Linewidth',plotLineWidth);
                            if holdOnState
                                set(tfplotHandles(k,j),'Nextplot','add');
                            else
                                set(tfplotHandles(k,j),'Nextplot','replace');
                            end

                        elseif (spectrumType == 2) % difference from baseline
                            if xsST == xsBL 
                                [SBL,f]=mtspectrumc((analogData(goodPos,BLPos))',mtmParams);
                                [SST,f]=mtspectrumc((analogData(goodPos,STPos))',mtmParams);
                                plot(tfplotHandles(k,j),f,10*(conv2Log(SST)-conv2Log(SBL)),'color',plotColor,'Linewidth',plotLineWidth);
                            else
                                disp('Choose same baseline and stimulus periods..');
                            end
                            if holdOnState
                                set(tfplotHandles(k,j),'Nextplot','add');
                            else
                                set(tfplotHandles(k,j),'Nextplot','replace');
                            end
                        end

                    %----------------------------------------------------------    
                    elseif (plotStyle == 1) % pcolor plot

                        if (spectrumType == 1) % raw
                            [S1,t,f]=mtspecgramc((analogData(goodPos,:))',movingWin,mtmParams);
                            t = t + timeVals(1); % shift the t values to the actual time
                            pcolor(tfplotHandles(k,j),t,f,conv2Log(S1')); shading(tfplotHandles(k,j),'interp');
                            caxis([cmin cmax]);

                        elseif (spectrumType == 2) % difference from baseline
                            [S1,t,f]=mtspecgramc((analogData(goodPos,:))',movingWin,mtmParams);
                            t = t + timeVals(1); % shift the t values to the actual time

                            % baseline period calculation
                            tBL = (t>=BLMin) & (t<=BLMax); % baseline time indices
                            S1BL = S1(tBL,:); % part of spectrum corresponding to the baseline period
                            mlogS1BL = mean(conv2Log(S1BL),1); % mean log power 
                            % across these time points at every frequency

                            % difference spectrum calculation
                            dS1 = 10*(conv2Log(S1) - repmat(mlogS1BL,size(S1,1),1)); % in dB 
                            % subtract the mean baseline power at each
                            % frequency from the raw power at that frequency at
                            % every time point

                            % plot the difference spectrum
                            pcolor(tfplotHandles(k,j),t,f,dS1'); shading(tfplotHandles(k,j),'interp');
                            caxis([cmin cmax]);
                        end

                    %----------------------------------------------------------    
                    elseif (plotStyle == 2) % imagesc plot

                        if (spectrumType == 1) % raw
                            [S1,t,f]=mtspecgramc((analogData(goodPos,:))',movingWin,mtmParams);
                            t = t + timeVals(1); % shift the t values to the actual time
                            imagesc(t,f,conv2Log(S1'),'Parent',tfplotHandles(k,j)); 
                            axis(tfplotHandles(k,j),'xy'); % for imagesc every 
                            % parameter has to be set using the plot handles 
                            % explicitly 
                            set(tfplotHandles(k,j),'cLim',[cmin cmax]);
                            
                            

                        elseif (spectrumType == 2) % difference from baseline
                            [S1,t,f]=mtspecgramc((analogData(goodPos,:))',movingWin,mtmParams);
                            t = t + timeVals(1); % shift the t values to the actual time

                            % baseline period calculation
                            tBL = (t>=BLMin) & (t<=BLMax); % baseline time indices
                            S1BL = S1(tBL,:); % part of spectrum corresponding to the baseline period
                            mlogS1BL = mean(conv2Log(S1BL),1); % mean log power 
                            % across these time points at every frequency

                            % difference spectrum calculation
                            dS1 = 10*(conv2Log(S1) - repmat(mlogS1BL,size(S1,1),1)); % in dB 
                            % subtract the mean baseline power at each
                            % frequency from the raw power at that frequency at
                            % every time point

                            % plot the difference spectrum                       
                            imagesc(t,f,dS1','Parent',tfplotHandles(k,j)); 
                            axis(tfplotHandles(k,j),'xy'); % for imagesc every 
                            % parameter has to be set using the plot handles 
                            % explicitly 
                            set(tfplotHandles(k,j),'cLim',[cmin cmax]);

                        end
                    end
                    
                    %======================================================

                elseif (tfMethod == 2) % MP Algorithm
                    
                    disp(['number of atoms for reconstruction:' num2str(numAtomsMP)]);
                    
                    % Vinay - create a directory (if not already present)
                    % to store the saved mp energy data
                    folderSaveMP = fullfile(folderName,'mpSavedSpectrum');
                    makeDirectory(folderSaveMP);                        
                    
%                     tagPos = (['Ca' num2str(a(3)) 'e' num2str(e(3)) 's' num2str(s(3)) 'f' num2str(f(3)) 'o' num2str(o(3)) 'c' num2str(c(3)) 't' num2str(t(3)) 'p' num2str(p(3)) 'r' num2str(r(3)) ...
%                                    '_Ra' num2str(a(2)) 'e' num2str(e(2)) 's' num2str(s(2)) 'f' num2str(f(2)) 'o' num2str(o(2)) 'c' num2str(c(2)) 't' num2str(t(2)) 'p' num2str(p(2)) 'r' num2str(r(2)) ...
%                                    '_Sa' num2str(a(1)) 'e' num2str(e(1)) 's' num2str(s(1)) 'f' num2str(f(1)) 'o' num2str(o(1)) 'c' num2str(c(1)) 't' num2str(t(1)) 'p' num2str(p(1)) 'r' num2str(r(1))]);
                    
                    if strncmp(protocolName,'GRF',3)
                        tagPos = ['L' tagPos1];
                    elseif strncmp(protocolName,'CRS',3)
                        tagPos = ['C' tagPos3 'R' tagPos2 'S' tagPos1];
                    end

                                                                          
                    if useBipolar
                        
                        folderSaveMP = fullfile(folderSaveMP,'bipolar');
                        makeDirectory(folderSaveMP);
                        
                        if strncmp(protocolName,'GRF',3)
                            mpSavedSpectrumName = (['mpSpectrum' num2str(numAtomsMP) 'atoms' '_notch' num2str(notchData)...
                                'pos' tagPos '_protoGRF'...
                                '_bipolar' num2str(channelString) num2str(analogChannelString2)...
                                '_vary' num2str(paramNum1) num2str(paramNum2)]);
                        elseif strncmp(protocolName,'CRS',3)
                            mpSavedSpectrumName = (['mpSpectrum' num2str(numAtomsMP) 'atoms' '_notch' num2str(notchData)...
                                'pos' tagPos '_proto' num2str(loadProtocolNumber)...
                                '_bipolar' num2str(channelString) num2str(analogChannelString2)...
                                '_vary' num2str(gaborNum1) num2str(paramNum1) num2str(gaborNum2) num2str(paramNum2)]);
                        end
                        
                        
                        % Vinay - save or load saved MP data only if the
                        % saveMPFlag is set, else just do all the
                        % computations
                        if (exist(fullfile(folderSaveMP,[mpSavedSpectrumName '.mat']),'file') && saveMPFlag)
                            disp('loading saved spectrum data...');
                            load(fullfile(folderSaveMP,[mpSavedSpectrumName '.mat']));
                            L = size(analogData,2);
                        else     
                            clear tag gaborInfo X
                            tag = 'temp/';

                            disp('Running MP for the bipolar case...');

                            % Import the data
                            X = analogData';
                            L = size(analogData,2);
                            signalRange = [1 L]; % full range
                            importData(X,mpFolderTemp,tag,signalRange,Fs);

                            % perform Gabor decomposition
                            Numb_points = L; % length of the signal
                            Max_iterations = 100; % number of iterations

                            disp(['MP decomposition with max iterations =' num2str(Max_iterations)]);

                            runGabor(mpFolderTemp,tag,Numb_points, Max_iterations);
                            disp('MP decomposition done');

                            gaborInfo = getGaborData(mpFolderTemp,tag,1);

                            gaborInfoGoodPos = gaborInfo(goodPos);

                            wrap = [];
                            atomList = (1:numAtomsMP);
                            mpSpectrum = [];
                            disp(['Reconstructing Energy from:' num2str(numAtomsMP) 'atoms, and'  num2str(length(goodPos)) 'trials']);
                            for m=1:length(goodPos)
                                disp(['trial number:' num2str(m) '(actual trial -)' num2str(goodPos(m))]);
                                rEnergy = reconstructEnergyFromAtomsMPP(gaborInfoGoodPos{m}.gaborData,L,wrap,atomList);
                                if m == 1
                                    mpSpectrum = rEnergy;
                                else
                                    mpSpectrum = mpSpectrum + rEnergy;
                                end
                            end
                            mpSpectrum = mpSpectrum/length(goodPos);
                        
                            if saveMPFlag
                                disp('saving spectrum data...');
                                save(fullfile(folderSaveMP,mpSavedSpectrumName), 'mpSpectrum');
                            end
                        end
                        
                    else
                        
                        folderSaveMP = fullfile(folderSaveMP,num2str(channelString));
                        makeDirectory(folderSaveMP);
                        
                        
                        if strncmp(protocolName,'GRF',3)
                            mpSavedSpectrumName = (['mpSpectrum' num2str(numAtomsMP) 'atoms' '_notch' num2str(notchData)...
                                'pos_' tagPos '_protoGRF'...
                                '_vary' num2str(paramNum1) num2str(paramNum2)]);
                        elseif strncmp(protocolName,'CRS',3)
                            mpSavedSpectrumName = (['mpSpectrum' num2str(numAtomsMP) 'atoms' '_notch' num2str(notchData)...
                                'pos_' tagPos '_proto' num2str(loadProtocolNumber)...
                                '_vary' num2str(gaborNum1) num2str(paramNum1) num2str(gaborNum2) num2str(paramNum2)]);
                        end
                        
                        
                        if (exist(fullfile(folderSaveMP,[mpSavedSpectrumName '.mat']),'file') && saveMPFlag)
                            disp('loading saved spectrum data...');
                            load(fullfile(folderSaveMP,[mpSavedSpectrumName '.mat']));
                            L = size(analogData,2);
                        else
                            clear tag gaborInfo
                            tag = channelString;
                            load(fullfile(mpFolder,tag,'gaborInfo.mat'));
                            gaborInfoGoodPos = gaborInfo(goodPos);

                            L = size(analogData,2);
                            wrap = [];
                            atomList = (1:numAtomsMP);
                            mpSpectrum = [];
                            disp(['Reconstructing Energy from:' num2str(numAtomsMP) 'atoms, and'  num2str(length(goodPos)) 'trials']);
                            for m=1:length(goodPos)
                                disp(['trial number:' num2str(m) '(actual trial -)' num2str(goodPos(m))]);
                                rEnergy = reconstructEnergyFromAtomsMPP(gaborInfoGoodPos{m}.gaborData,L,wrap,atomList);
                                if m == 1
                                    mpSpectrum = rEnergy;
                                else
                                    mpSpectrum = mpSpectrum + rEnergy;
                                end
                            end
                            mpSpectrum = mpSpectrum/length(goodPos);
                            
                            if saveMPFlag
                                disp('saving spectrum data...');
                                save(fullfile(folderSaveMP,mpSavedSpectrumName), 'mpSpectrum');
                            end
                        
                        end
                    end
                    
                    mpSpectrum = mpSpectrum'; % transpose so that first index corresponds to 'time' and second to 'freq'
                    
                    t = timeVals;
                    f = 0:Fs/L:Fs/2; 
                    
                    % plot the MP Spectrum
                    if plotStyle == 3 % line plot
                        disp('No line plot for MP method');
                        
                    elseif plotStyle == 1 % pcolor plot
                        
                        if (spectrumType == 1) % raw
                            pcolor(tfplotHandles(k,j),t,f,conv2Log(mpSpectrum')); 
                            shading(tfplotHandles(k,j),'interp');
                            caxis([cmin cmax]);
                            
                        elseif(spectrumType == 2) % difference from baseline                            
                            % baseline period calculation
                            tBL = (t>=BLMin) & (t<=BLMax); % baseline time indices
                            mpSpectrumBL = mpSpectrum(tBL,:); % part of spectrum corresponding to the baseline period
                            mlogmpSpectrumBL = mean(conv2Log(mpSpectrumBL),1); % mean log power 
                            % across these time points at every frequency

                            % difference spectrum calculation
                            dmpSpectrum = 10*(conv2Log(mpSpectrum) - repmat(mlogmpSpectrumBL,size(mpSpectrum,1),1)); % in dB 
                            % subtract the mean baseline power at each
                            % frequency from the raw power at that frequency at
                            % every time point

                            % plot the difference spectrum
                            pcolor(tfplotHandles(k,j),t,f,dmpSpectrum'); 
                            shading(tfplotHandles(k,j),'interp');
                            caxis([cmin cmax]);
                        end
                        
                    elseif plotStyle == 2 % imagesc plot
                        
                        if (spectrumType == 1) % raw
                            imagesc(t,f,conv2Log(mpSpectrum'),'Parent',tfplotHandles(k,j)); 
                            axis(tfplotHandles(k,j),'xy'); % for imagesc every 
                            % parameter has to be set using the plot handles 
                            % explicitly 
                            set(tfplotHandles(k,j),'cLim',[cmin cmax]);
                            
                        elseif(spectrumType == 2) % difference from baseline                            
                            % baseline period calculation
                            tBL = (t>=BLMin) & (t<=BLMax); % baseline time indices
                            mpSpectrumBL = mpSpectrum(tBL,:); % part of spectrum corresponding to the baseline period
                            mlogmpSpectrumBL = mean(conv2Log(mpSpectrumBL),1); % mean log power 
                            % across these time points at every frequency

                            % difference spectrum calculation
                            dmpSpectrum = 10*(conv2Log(mpSpectrum) - repmat(mlogmpSpectrumBL,size(mpSpectrum,1),1)); 
                            % subtract the mean baseline power at each
                            % frequency from the raw power at that frequency at
                            % every time point

                            % plot the difference spectrum
                            imagesc(t,f,dmpSpectrum','Parent',tfplotHandles(k,j)); 
                            axis(tfplotHandles(k,j),'xy'); % for imagesc every 
                            % parameter has to be set using the plot handles 
                            % explicitly 
                            set(tfplotHandles(k,j),'cLim',[cmin cmax]);
                        end
                        
                    end
                    
%------------------------------------------------------------------------                  
                    elseif (tfMethod == 3) % HHT
                    
                    disp(['HHT calculation...Noise ratio: ' num2str(Nstd) ', Ensemble number: ' num2str(NE)]);
                    
                    % Vinay - create a directory (if not already present)
                    % to store the saved mp energy data
                    folderSaveHHT = fullfile(folderName,'hhtSavedSpectrum');
                    makeDirectory(folderSaveHHT);                        
                    
                    % Required parameters for implementation
                    t0=timeVals(1);         % starting time value of the signal
                    t1=timeVals(end);       % last time value of the signal
                    fw0=0;                  % minimum frequency
                    fw1=250;                % maximum frequency
                    fres=251;               % frequency resolution i.e. no. of frequencies between the maximum and minimum frequencies
                    tres=length(timeVals);  % time resolution i.e. no. of time points between the maximum and minimum time values
                    tw0=t0;                 % minimum time for zooming
                    tw1=t1;                 % maximum time for zooming
                    lscale=[];              % => linear axis
                    nfilter=0;              % number of points to be median filtered
                    numIMFs = 9;            % number of IMFs
                    ifmethod='hilbert';     % instantaneous frequency method - default
                    normmethod='hilbert';   % normalization method - default
                    fastHHT = 1;            % Use FEEMD if set, EEMD if 0
                    
%                     tagPos = (['Ca' num2str(a(3)) 'e' num2str(e(3)) 's' num2str(s(3)) 'f' num2str(f(3)) 'o' num2str(o(3)) 'c' num2str(c(3)) 't' num2str(t(3)) 'p' num2str(p(3)) 'r' num2str(r(3)) ...
%                                    '_Ra' num2str(a(2)) 'e' num2str(e(2)) 's' num2str(s(2)) 'f' num2str(f(2)) 'o' num2str(o(2)) 'c' num2str(c(2)) 't' num2str(t(2)) 'p' num2str(p(2)) 'r' num2str(r(2)) ...
%                                    '_Sa' num2str(a(1)) 'e' num2str(e(1)) 's' num2str(s(1)) 'f' num2str(f(1)) 'o' num2str(o(1)) 'c' num2str(c(1)) 't' num2str(t(1)) 'p' num2str(p(1)) 'r' num2str(r(1))]);
                    
                    if strncmp(protocolName,'GRF',3)
                        tagPos = ['L' tagPos1];
                    elseif strncmp(protocolName,'CRS',3)
                        tagPos = ['C' tagPos3 'R' tagPos2 'S' tagPos1];
                    end
                    
                    clear E
                                                                          
                    if useBipolar
                        
                        folderSaveHHT = fullfile(folderSaveHHT,'bipolar');
                        makeDirectory(folderSaveHHT);
                        
                        
                        if strncmp(protocolName,'GRF',3)
                            hhtSavedSpectrumName = (['hhtSpectrum_Nstd' num2str(Nstd) '_NE' num2str(NE) '_numIMFs' num2str(numIMFs) '_notch' num2str(notchData)...
                                'pos' tagPos '_protoGRF'...
                                '_bipolar' num2str(channelString) num2str(analogChannelString2)...
                                '_vary' num2str(paramNum1) num2str(paramNum2)]);
                        elseif strncmp(protocolName,'CRS',3)
                            hhtSavedSpectrumName = (['hhtSpectrum_Nstd' num2str(Nstd) '_NE' num2str(NE) '_numIMFs' num2str(numIMFs) '_notch' num2str(notchData)...
                                'pos' tagPos '_proto' num2str(protocolNumber)...
                                '_bipolar' num2str(channelString) num2str(analogChannelString2)...
                                '_vary' num2str(gaborNum1) num2str(paramNum1) num2str(gaborNum2) num2str(paramNum2)]);
                        end
                        
                        
                        % Vinay - save or load saved HHT data only if the
                        % saveHHTFlag is set, else just do all the
                        % computations
                        if (exist(fullfile(folderSaveHHT,[hhtSavedSpectrumName '.mat']),'file') && saveHHTFlag)
                            disp('loading saved spectrum data...');
                            load(fullfile(folderSaveHHT,[hhtSavedSpectrumName '.mat']));
                        
                        else
                            disp(['Computing HHT spectrum...no. of trials: ' num2str(numTrials)]);
                            % Initialize the energy distribution matrices
                            E = zeros(fres,tres);
                            numTrials = length(goodPos);
                            for hi=1:numTrials
                                disp(['trial number:' num2str(hi) '(actual trial -)' num2str(goodPos(hi))]);
                                if ~fastHHT
                                    imfList = eemd(analogData(goodPos(hi),:),Nstd,NE);
                                else
                                    imfList = rcada_eemd(analogData(goodPos(hi),:),Nstd,NE,numIMFs);
                                end
                                [P,t,f] = nnspe(imfList(:,2:(numIMFs+1)),t0,t1,fres,tres,fw0,fw1,tw0,tw1,ifmethod,normmethod,nfilter,lscale);
                                E=E+P;
                            end
                            
                            E = E/numTrials; % averaged across trials

                            if saveHHTFlag
                                save(fullfile(folderSaveHHT,[hhtSavedSpectrumName '.mat']),'E','t','f'); % save the HHT spectrum matrices
                            end
                            
                        end
                        
                    else
                        
                        folderSaveHHT = fullfile(folderSaveHHT,num2str(channelString));
                        makeDirectory(folderSaveHHT);
                        
                        
                        if strncmp(protocolName,'GRF',3)
                            hhtSavedSpectrumName = (['hhtSpectrum_Nstd' num2str(Nstd) '_NE' num2str(NE) '_numIMFs' num2str(numIMFs) '_notch' num2str(notchData)...
                                'pos_' tagPos '_protoGRF'...
                                '_vary' num2str(paramNum1) num2str(paramNum2)]);
                        elseif strncmp(protocolName,'CRS',3)
                            hhtSavedSpectrumName = (['hhtSpectrum_Nstd' num2str(Nstd) '_NE' num2str(NE) '_numIMFs' num2str(numIMFs) '_notch' num2str(notchData)...
                                'pos_' tagPos '_proto' num2str(protocolNumber)...
                                '_vary' num2str(gaborNum1) num2str(paramNum1) num2str(gaborNum2) num2str(paramNum2)]);
                        end
                        
                        
                        if (exist(fullfile(folderSaveHHT,[hhtSavedSpectrumName '.mat']),'file') && saveHHTFlag)
                            disp('loading saved spectrum data...');
                            load(fullfile(folderSaveHHT,[hhtSavedSpectrumName '.mat']));
                        else
                            
                            % Initialize the energy distribution matrices
                            E = zeros(fres,tres);
                            numTrials = length(goodPos);
                            disp(['Computing HHT spectrum...no. of trials: ' num2str(numTrials)]);
                            numIMFs = 7;
                            for hi=1:numTrials
                                disp(['trial number:' num2str(hi) '(actual trial -)' num2str(goodPos(hi))]);
                                if ~fastHHT
                                    imfList = eemd(analogData(goodPos(hi),:),Nstd,NE);
                                    [P,t,f] = nnspe(imfList(:,2:(numIMFs+1)),t0,t1,fres,tres,fw0,fw1,tw0,tw1,ifmethod,normmethod,nfilter,lscale);
                                else
                                    imfList = rcada_eemd(analogData(goodPos(hi),:),Nstd,NE,numIMFs);
                                    imfList = imfList';
                                    [P,t,f] = nnspe(imfList(:,1:numIMFs),t0,t1,fres,tres,fw0,fw1,tw0,tw1,ifmethod,normmethod,nfilter,lscale);
                                end
                                
                                E=E+P;
                            end
                            
                            E = E/numTrials; % averaged across trials

                            if saveHHTFlag
                                save(fullfile(folderSaveHHT,[hhtSavedSpectrumName '.mat']),'E','t','f'); % save the HHT spectrum matrices
                            end
                        
                        end
                    end
                    
                    
                    if gaussFtr~=0
                        %%%%%%%%%%%%%%%%%____Smoothing of the HHT spectra____%%%%%%%%%%%%%%%%%%%%%%

                        gaussFilter = fspecial('gaussian', gaussFtr, 1); % a Gaussian low-pass filter

                        % filter the energy distribution E
                        E = filter2(gaussFilter, E);
                    end
                    
                    E = E'; % transpose for plotting

                    
                    % plot the HHT Spectrum
                    if plotStyle == 3 % line plot
%                         analogData = mean(analogData(goodPos,:),1); %
%                         Vinay - can't do this because this changes
%                         analogData itself and we are not reloading it
%                         everytime. Instead do as follows
                        analogDataMean = mean(analogData(goodPos,:),1);
                        disp('Computing IMFs and displaying them...');
                        if ~fastHHT
                            imfList = eemd(analogDataMean,Nstd,NE);
                        else
                            imfList = rcada_eemd(analogDataMean,Nstd,NE,numIMFs);
                        end
%                         imfList = imfList';
                        
%                         figure(3);
                        nimf = 4;
                        IMFcolors = {'r','c','b','k'};
                        for m=1:nimf
%                           subplot(nimf,1,m); 
                          plot(tfplotHandles(k,j),t,imfList(m,:),'color',IMFcolors{m},'LineWidth',1.5);
                          set(tfplotHandles(k,j),'Nextplot','add');
                          pause;
                        end
%                         cFigC = get(cFig,'all');
%                         tfC = get(tfplotHandles(k,j),'Children');
%                         copyobj(cFig,tfplotHandles(k,j));
%                         subplot(nimf,1,1); title('IMFs');

                            if holdOnState
                                set(tfplotHandles(k,j),'Nextplot','add');
                            else
                                set(tfplotHandles(k,j),'Nextplot','replace');
                            end
                    
                        %--------------------------------------------------
                    
                    elseif plotStyle == 1 % pcolor plot
                        
                        if (spectrumType == 1) % raw
                            pcolor(tfplotHandles(k,j),t,f,conv2Log(E)); 
                            shading(tfplotHandles(k,j),'interp');
                            caxis([cmin cmax]);
                            
                        elseif(spectrumType == 2) % difference from baseline                            
                            % baseline period calculation
                            tBL = (t>=BLMin) & (t<=BLMax); % baseline time indices
                            EBL = E(tBL,:); % part of spectrum corresponding to the baseline period
                            mlogEBL = mean(conv2Log(EBL),1); % mean log power 
                            % across these time points at every frequency

                            % difference spectrum calculation
                            dhhtSpectrum = 10*(conv2Log(E) - repmat(mlogEBL,size(E,1),1)); % in dB 
                            % subtract the mean baseline power at each
                            % frequency from the raw power at that frequency at
                            % every time point

                            % plot the difference spectrum
                            pcolor(tfplotHandles(k,j),t,f,dhhtSpectrum); 
                            shading(tfplotHandles(k,j),'interp');
                            caxis([cmin cmax]);
                        end
                        
                        %--------------------------------------------------
                        
                    elseif plotStyle == 2 % imagesc plot
                        
                        if (spectrumType == 1) % raw
                            imagesc(t,f,conv2Log(E'),'Parent',tfplotHandles(k,j)); 
                            axis(tfplotHandles(k,j),'xy'); % for imagesc every 
                            % parameter has to be set using the plot handles 
                            % explicitly 
                            set(tfplotHandles(k,j),'cLim',[cmin cmax]);
                            
                        elseif(spectrumType == 2) % difference from baseline                            
                            % baseline period calculation
                            tBL = (t>=BLMin) & (t<=BLMax); % baseline time indices
                            EBL = E(tBL,:); % part of spectrum corresponding to the baseline period
                            mlogEBL = mean(conv2Log(EBL),1); % mean log power 
                            % across these time points at every frequency

                            % difference spectrum calculation
                            dhhtSpectrum = 10*(conv2Log(E) - repmat(mlogEBL,size(E,1),1)); 
                            % subtract the mean baseline power at each
                            % frequency from the raw power at that frequency at
                            % every time point

                            % plot the difference spectrum
                            imagesc(t,f,dhhtSpectrum','Parent',tfplotHandles(k,j)); 
                            axis(tfplotHandles(k,j),'xy'); % for imagesc every 
                            % parameter has to be set using the plot handles 
                            % explicitly 
                            set(tfplotHandles(k,j),'cLim',[cmin cmax]);
                        end
                        
                    end
                   
                end
                
                if strncmp(protocolName,'CRS',3)
                    titleGabor1 = titleGabor(gaborNum1);
                    titleGabor2 = titleGabor(gaborNum2);
                elseif strncmp(protocolName,'GRF',3)
                    titleGabor1 = 'gabor';
                    titleGabor2 = 'gabor';
                end

                % Display title
                if (j==1 && k==1)
                    title(tfplotHandles(k,j),[titleGabor1 '-' titleParam1 'vs' titleGabor2 '-' titleParam2 ':' num2str(titleList1(k)) ',' num2str(titleList2(j))],'FontSize',titleFontSize);
                else
                    title(tfplotHandles(k,j),[num2str(titleList1(k)) ',' num2str(titleList2(j))],'FontSize',titleFontSize);
                end
                
                
                % show number of trials for this condition
                text(0.1,0.9,['n = ' num2str(length(goodPos))],'unit','normalized','fontsize',7,'Parent',tfplotHandles(k,j));
                
            end  
            end
        end
        
end


%--------------------------------------------------------------------------
% Reference scheme functions
%--------------------------------------------------------------------------

function elecBadTrials = getBipolarBadTrials(channelString,allBadTrials)
    switch channelString
        case 'elec12'
            elecBadTrials = unique([allBadTrials{1} allBadTrials{2}]);
        case 'elec34'
            elecBadTrials = unique([allBadTrials{3} allBadTrials{4}]);
        case 'elec56'
            elecBadTrials = unique([allBadTrials{5} allBadTrials{6}]);
        case 'elec78'
            elecBadTrials = unique([allBadTrials{7} allBadTrials{8}]);
        case 'elec910'
            elecBadTrials = unique([allBadTrials{9} allBadTrials{10}]);
        case 'elec1112'
            elecBadTrials = unique([allBadTrials{11} allBadTrials{12}]);
        case 'elec1314'
            elecBadTrials = unique([allBadTrials{13} allBadTrials{14}]);
        case 'elec1516'
            elecBadTrials = unique([allBadTrials{15} allBadTrials{16}]);
        case 'elec1719'
            elecBadTrials = unique([allBadTrials{17} allBadTrials{19}]);
        case 'elec2021'
            elecBadTrials = unique([allBadTrials{20} allBadTrials{21}]);
        case 'elec118'
            elecBadTrials = unique([allBadTrials{1} allBadTrials{18}]);
        case 'elec218'
            elecBadTrials = unique([allBadTrials{2} allBadTrials{18}]);
        case 'elec318'
            elecBadTrials = unique([allBadTrials{3} allBadTrials{18}]);
        case 'elec418'
            elecBadTrials = unique([allBadTrials{4} allBadTrials{18}]);
        case 'elec518'
            elecBadTrials = unique([allBadTrials{5} allBadTrials{18}]);
        case 'elec618'
            elecBadTrials = unique([allBadTrials{6} allBadTrials{18}]);
        case 'elec718'
            elecBadTrials = unique([allBadTrials{7} allBadTrials{18}]);
        case 'elec818'
            elecBadTrials = unique([allBadTrials{8} allBadTrials{18}]);
        case 'elec918'
            elecBadTrials = unique([allBadTrials{9} allBadTrials{18}]);
        case 'elec1018'
            elecBadTrials = unique([allBadTrials{10} allBadTrials{18}]);
        case 'elec1118'
            elecBadTrials = unique([allBadTrials{11} allBadTrials{18}]);
        case 'elec1218'
            elecBadTrials = unique([allBadTrials{12} allBadTrials{18}]);
        case 'elec1318'
            elecBadTrials = unique([allBadTrials{13} allBadTrials{18}]);
        case 'elec1418'
            elecBadTrials = unique([allBadTrials{14} allBadTrials{18}]);
        case 'elec1518'
            elecBadTrials = unique([allBadTrials{15} allBadTrials{18}]);
        case 'elec1618'
            elecBadTrials = unique([allBadTrials{16} allBadTrials{18}]);
        case 'elec1718'
            elecBadTrials = unique([allBadTrials{17} allBadTrials{18}]);
        case 'elec1818'
            elecBadTrials = unique([allBadTrials{18} allBadTrials{18}]);
        case 'elec1918'
            elecBadTrials = unique([allBadTrials{19} allBadTrials{18}]);
        case 'elec2018'
            elecBadTrials = unique([allBadTrials{20} allBadTrials{18}]);
        case 'elec2118'
            elecBadTrials = unique([allBadTrials{21} allBadTrials{18}]);
    end
end

%--------------------------------------------------------------------------

function reftfplotLFPDataVaryParameters1Channel(tfplotHandles,channelString,analogChannelString2,a,e,s,f,o,c,t,r,p,gaborNum1,paramNum1,gaborNum2,paramNum2,folderLFP,...
            timeVals,plotColor,BLMin,BLMax,STMin,STMax,folderName, protocolNumber, notchData,useDiff,...
            tfMethod,mtmParams,movingWin,numAtomsMP,plotStyle,spectrumType,cmin,cmax,holdOnState,saveMPFlag,loadProtocolNumber,plotLineWidth,...
            Nstd,NE,gaussFtr,saveHHTFlag,protocolName,useAllBadTrials,...
            refType,refCombineType,refDiffType,numElectrodes,intersectTrials,...
            takeLogTrial,existsBipolarData,contralateralNeighbour)
        
        folderExtract = fullfile(folderName,'extractedData');
        folderSegment = fullfile(folderName,'segmentedData');
        mpFolder = fullfile(folderName,'mpAnalysis'); % for using stored MP data    
        folderBipolar = fullfile(folderLFP,'bipolar');
        folderAverage = folderLFP;
        folderCSD = fullfile(folderLFP,'csd');

        mpFolderBipolar = fullfile(mpFolder,'bipolar');
        mpFolderAverage = fullfile(mpFolder,'average');
        mpFolderCSD = fullfile(mpFolder,'csd');
        
        if ispc
            mpFolderTemp = 'C:\Documents\MP\data\'; % for online calculations
        else
            mpFolderTemp = '/home/vinay/Documents/MP/data/'; % for online calculations
        end


        titleFontSize = 9;

%         timeForComputation = [40 100]/1000; % ms
%         freqForComputation = [40 60]; % Hz

        [parameterCombinations,aValsUnique,eValsUnique,sValsUnique,...
            fValsUnique,oValsUnique,cValsUnique,tValsUnique,rValsUnique,pValsUnique] = loadParameterCombinations(folderExtract);
        
        numRows = size(tfplotHandles,1);
        numCols = size(tfplotHandles,2);
        
        % Get the data
        clear signal analogData
%         load([folderLFP channelString]);
        

        % Get bad trials
        badTrialFile = fullfile(folderSegment,'badTrials.mat');
        existsBadTrialFile = 0;
        if ~exist(badTrialFile,'file')
            disp('Bad trial file does not exist...');
            badTrials=[];
        else
            [allBadTrials, badTrials, nameElec] = loadBadTrials(badTrialFile);
            disp([num2str(length(badTrials)) ' bad trials']);
            existsBadTrialFile = 1;
        end
        

        % Choosing the good position trials
if strncmp(protocolName,'GRF',3)
    % [Vinay] - repeat the set of parameters for each gabor depending on the
    % number of columns to be drawn
    aList = repmat(a,numRows,numCols);
    eList = repmat(e,numRows,numCols);
    sList = repmat(s,numRows,numCols);
    fList = repmat(f,numRows,numCols);
    oList = repmat(o,numRows,numCols);
    cList = repmat(c,numRows,numCols);
    tList = repmat(t,numRows,numCols);


        % [Vinay] decide the row parameter based on paramNum1
            for row = 1:numRows
                switch (paramNum1-1)
                    case 1
                        aList(row,:) = row;                    
                        % Every new row takes a new value of param1.
                        % So assign the 'row' value in the list above
                        titleParam1 = 'Azi: ';
                        titleList1 = aValsUnique;
                    case 2
                        eList(row,:) = row; 
                        titleParam1 = 'Ele: ';
                        titleList1 = eValsUnique;
                    case 3
                        sList(row,:) = row; 
                        titleParam1 = 'Sigma: ';
                        titleList1 = sValsUnique;
                    case 4
                        fList(row,:) = row; 
                        titleParam1 = 'SF: ';
                        titleList1 = fValsUnique;
                    case 5
                        oList(row,:) = row; 
                        titleParam1 = 'Ori: ';
                        titleList1 = oValsUnique;
                    case 6
                        cList(row,:) = row; 
                        titleParam1 = 'Contr: ';
                        titleList1 = cValsUnique;
                    case 7
                        tList(row,:) = row; 
                        titleParam1 = 'TF: ';
                        titleList1 = tValsUnique;
                    otherwise
                        disp('No particular gabor or parameter selected');
                end
            end

            % [Vinay] decide the column parameter based on paramNum2
            for row = 1:numRows
                switch (paramNum2-1)
                    case 1
                        aList(row,:) = 1:length(aValsUnique);
                        % For every row we now have to assign indices
                        % incrementing along the columns.
                        % Basically the column entries go from 1 to 
                        % length{nValsUnique}
                        titleParam2 = 'Azi: ';
                        titleList2 = aValsUnique;
                    case 2
                        eList(row,:) = 1:length(eValsUnique); 
                        titleParam2 = 'Ele: ';
                        titleList2 = eValsUnique;
                    case 3
                        sList(row,:) = 1:length(sValsUnique);
                        titleParam2 = 'Sigma: ';
                        titleList2 = sValsUnique;
                    case 4
                        fList(row,:) = 1:length(fValsUnique);
                        titleParam2 = 'SF: ';
                        titleList2 = fValsUnique;
                    case 5
                        oList(row,:) = 1:length(oValsUnique);
                        titleParam2 = 'Ori: ';
                        titleList2 = oValsUnique;
                    case 6
                        cList(row,:) = 1:length(cValsUnique);
                        titleParam2 = 'Contr: ';
                        titleList2 = cValsUnique;
                    case 7
                        tList(row,:) = 1:length(tValsUnique);
                        titleParam2 = 'TF: ';
                        titleList2 = tValsUnique;
                    otherwise
                        disp('No particular gabor or parameter selected');
                end
            end


    % [Vinay] - Get the lengths of indices in parameterCombinations

    aLen = length(aValsUnique);
    eLen = length(eValsUnique);
    sLen = length(sValsUnique);
    fLen = length(fValsUnique);
    oLen = length(oValsUnique);
    cLen = length(cValsUnique);
    tLen = length(tValsUnique);

    % If more than one value, then length is one greater for all the values
    % together
    if (aLen> 1)           aLen=aLen+1;                    end
    if (eLen> 1)           eLen=eLen+1;                    end
    if (sLen> 1)           sLen=sLen+1;                    end
    if (fLen> 1)           fLen=fLen+1;                    end
    if (oLen> 1)           oLen=oLen+1;                    end
    if (cLen> 1)           cLen=cLen+1;                    end
    if (tLen> 1)           tLen=tLen+1;                    end


elseif strncmp(protocolName,'CRS',3)

        % [Vinay] - repeat the set of parameters for each gabor depending on the
        % number of rows and columns to be drawn
        aList = repmat(a,numRows,numCols);
        eList = repmat(e,numRows,numCols);
        sList = repmat(s,numRows,numCols);
        fList = repmat(f,numRows,numCols);
        oList = repmat(o,numRows,numCols);
        cList = repmat(c,numRows,numCols);
        tList = repmat(t,numRows,numCols);
        rList = repmat(r,numRows,numCols);
        pList = repmat(p,numRows,numCols);

        % [Vinay] decide the row parameter based on paramNum1
        for row = 1:numRows
            switch (paramNum1-1)
                case 1
                    aList(gaborNum1+((row-1)*3),:) = row;                    
                    % along the rows, we get three consecutive entries for the 3 
                    % gabors. Therefore to go to the next index of the same gabor,
                    % we go ahead by 3 indices
                    % For gabor1 every new row takes a new value of param1.
                    % So assign the 'row' value in the list above
                    titleParam1 = 'Azi: ';
                    titleList1 = aValsUnique{gaborNum1};
                case 2
                    eList(gaborNum1+((row-1)*3),:) = row; 
                    titleParam1 = 'Ele: ';
                    titleList1 = eValsUnique{gaborNum1};
                case 3
                    sList(gaborNum1+((row-1)*3),:) = row; 
                    titleParam1 = 'Sigma: ';
                    titleList1 = sValsUnique{gaborNum1};
                case 4
                    fList(gaborNum1+((row-1)*3),:) = row; 
                    titleParam1 = 'SF: ';
                    titleList1 = fValsUnique{gaborNum1};
                case 5
                    oList(gaborNum1+((row-1)*3),:) = row; 
                    titleParam1 = 'Ori: ';
                    titleList1 = oValsUnique{gaborNum1};
                case 6
                    cList(gaborNum1+((row-1)*3),:) = row; 
                    titleParam1 = 'Contr: ';
                    titleList1 = cValsUnique{gaborNum1};
                case 7
                    tList(gaborNum1+((row-1)*3),:) = row; 
                    titleParam1 = 'TF: ';
                    titleList1 = tValsUnique{gaborNum1};
                case 8
                    rList(gaborNum1+((row-1)*3),:) = row; 
                    titleParam1 = 'Rad: ';
                    titleList1 = rValsUnique{gaborNum1};
                case 9
                    pList(gaborNum1+((row-1)*3),:) = row; 
                    titleParam1 = 'SPhs: ';
                    titleList1 = pValsUnique{gaborNum1};
                otherwise
                    disp('No particular gabor or parameter selected');
            end
        end


        % [Vinay] decide the column parameter based on paramNum2
        for row = 1:numRows
            switch (paramNum2-1)
                case 1
                    aList(gaborNum2+((row-1)*3),:) = 1:length(aValsUnique{gaborNum2});
                    % along the rows, we get three consecutive entries for the 3 
                    % gabors. Therefore to go to the next index of the same gabor,
                    % we go ahead by 3 indices
                    % Besides for every row we now have to assign indices
                    % incrementing along the columns.
                    % Basically the column entries go from 1 to 
                    % length{nValsUnique}
                    titleParam2 = 'Azi: ';
                    titleList2 = aValsUnique{gaborNum2};
                case 2
                    eList(gaborNum2+((row-1)*3),:) = 1:length(eValsUnique{gaborNum2}); 
                    titleParam2 = 'Ele: ';
                    titleList2 = eValsUnique{gaborNum2};
                case 3
                    sList(gaborNum2+((row-1)*3),:) = 1:length(sValsUnique{gaborNum2});
                    titleParam2 = 'Sigma: ';
                    titleList2 = sValsUnique{gaborNum2};
                case 4
                    fList(gaborNum2+((row-1)*3),:) = 1:length(fValsUnique{gaborNum2});
                    titleParam2 = 'SF: ';
                    titleList2 = fValsUnique{gaborNum2};
                case 5
                    oList(gaborNum2+((row-1)*3),:) = 1:length(oValsUnique{gaborNum2});
                    titleParam2 = 'Ori: ';
                    titleList2 = oValsUnique{gaborNum2};
                case 6
                    cList(gaborNum2+((row-1)*3),:) = 1:length(cValsUnique{gaborNum2});
                    titleParam2 = 'Contr: ';
                    titleList2 = cValsUnique{gaborNum2};
                case 7
                    tList(gaborNum2+((row-1)*3),:) = 1:length(tValsUnique{gaborNum2});
                    titleParam2 = 'TF: ';
                    titleList2 = tValsUnique{gaborNum2};
                case 8
                    rList(gaborNum2+((row-1)*3),:) = 1:length(rValsUnique{gaborNum2});
                    titleParam2 = 'Rad: ';
                    titleList2 = rValsUnique{gaborNum2};
                case 9
                    pList(gaborNum2+((row-1)*3),:) = 1:length(pValsUnique{gaborNum2});
                    titleParam2 = 'SPhs: ';
                    titleList2 = pValsUnique{gaborNum2};
                otherwise
                    disp('No particular gabor or parameter selected');
            end
        end

        % [Vinay] - Get the lengths of indices in parameterCombinations
        for i = 1:3
            aLen(i) = length(aValsUnique{i});
            eLen(i) = length(eValsUnique{i});
            sLen(i) = length(sValsUnique{i});
            fLen(i) = length(fValsUnique{i});
            oLen(i) = length(oValsUnique{i});
            cLen(i) = length(cValsUnique{i});
            tLen(i) = length(tValsUnique{i});
            pLen(i) = length(pValsUnique{i}); % Vinay - for CRS
            rLen(i) = length(rValsUnique{i}); % Vinay - for CRS

            % If more than one value, then length is one greater for all the values
            % together
            if (aLen(i)> 1)           aLen(i)=aLen(i)+1;                    end
            if (eLen(i)> 1)           eLen(i)=eLen(i)+1;                    end
            if (sLen(i)> 1)           sLen(i)=sLen(i)+1;                    end
            if (fLen(i)> 1)           fLen(i)=fLen(i)+1;                    end
            if (oLen(i)> 1)           oLen(i)=oLen(i)+1;                    end
            if (cLen(i)> 1)           cLen(i)=cLen(i)+1;                    end
            if (tLen(i)> 1)           tLen(i)=tLen(i)+1;                    end
            if (pLen(i)> 1)           pLen(i)=pLen(i)+1;                    end % Vinay - for CRS
            if (rLen(i)> 1)           rLen(i)=rLen(i)+1;                    end % Vinay - for CRS
        end
        
end

        % Main loop
        for k=1:numRows
            for j = 1:numCols
                
                if strncmp(protocolName,'GRF',3)
                    clear goodPos
                    goodPos = parameterCombinations{aList(k,j),eList(k,j),sList(k,j),fList(k,j),oList(k,j),cList(k,j),tList(k,j)};
            
                    tagPos1 = [num2str(aList(k,j)) num2str(eList(k,j)) num2str(sList(k,j))...
                               num2str(fList(k,j)) num2str(oList(k,j)) num2str(cList(k,j))...
                               num2str(tList(k,j))];
                           

                elseif strncmp(protocolName,'CRS',3)
                
                clear goodPos

                % [Vinay] - goodPos will be obtained by taking the goodPos for all
                % the three gabors and finding the common trials in them
                clear pos1 pos2 pos3
                switch protocolNumber
                           case 1 % Ring protocol
                               pos1 = parameterCombinations{aList(1+(k-1)*3,j),eList(1+(k-1)*3,j),sList(1+(k-1)*3,j),fList(1+(k-1)*3,j),oList(1+(k-1)*3,j),cList(1+(k-1)*3,j),tList(1+(k-1)*3,j),pList(1+(k-1)*3,j),rList(1+(k-1)*3,j),1}; 
                               % S params define C & S
                               % S i.e. kGabor0 corresponds to row number 1,4,7,..
                               % and so on. Hence it is 1+(k-1)*3
                               pos2 = parameterCombinations{aLen(2),eLen(2),sLen(2),fLen(2),oLen(2),cLen(2),tLen(2),pLen(2),rList(2+(k-1)*3,j),2}; % for R - only rad is imp
                               % R i.e. kGabor1 corresponds to row number 2,5,8,..
                               % and so on. Hence it is 2+(k-1)*3

                               % [Vinay] -  C is matched to S and therefore the entries
                               % for gabor3 i.e. can be ignored. Otherwise this can 
                               % give spurious results if the param values for C are 
                               % not set equal to those for S. Hence pos3 is the full
                               % set here
                               pos3 = parameterCombinations{aLen(3),eLen(3),sLen(3),fLen(3),oLen(3),cLen(3),tLen(3),pLen(3),rList(3+(k-1)*3,j),3}; % for C
                               % also keeping r variable here, just in case you
                               % have multiple r values
                               
                               % keep the position tags for storing MP data
                               tagPos1 = [num2str(aList(1+(k-1)*3,j)) num2str(eList(1+(k-1)*3,j)) num2str(sList(1+(k-1)*3,j))...
                                   num2str(fList(1+(k-1)*3,j)) num2str(oList(1+(k-1)*3,j)) num2str(cList(1+(k-1)*3,j))...
                                   num2str(tList(1+(k-1)*3,j)) num2str(pList(1+(k-1)*3,j)) num2str(rList(1+(k-1)*3,j))];
                               tagPos2 = [num2str(aLen(2)) num2str(eLen(2)) num2str(sLen(2)) num2str(fLen(2))...
                                   num2str(oLen(2)) num2str(cLen(2)) num2str(tLen(2)) num2str(pLen(2)) num2str(rList(2+(k-1)*3,j))];
                               tagPos3 = [num2str(aLen(3)) num2str(eLen(3)) num2str(sLen(3)) num2str(fLen(3))...
                                   num2str(oLen(3)) num2str(cLen(3)) num2str(tLen(3)) num2str(pLen(3)) num2str(rList(3+(k-1)*3,j))]; % for C

                           case 2 % Contrast Ring Protocol
                               pos1 = parameterCombinations{aList(1+(k-1)*3,j),eList(1+(k-1)*3,j),sList(1+(k-1)*3,j),fList(1+(k-1)*3,j),oList(1+(k-1)*3,j),cList(1+(k-1)*3,j),tList(1+(k-1)*3,j),pList(1+(k-1)*3,j),rList(1+(k-1)*3,j),1};
                               % S params define C & S
                               % S i.e. kGabor0 corresponds to row number 1,4,7,..
                               % and so on. Hence it is 1+(k-1)*3
                               pos2 = parameterCombinations{aLen(2),eLen(2),sLen(2),fLen(2),oLen(2),cList(2+(k-1)*3,j),tLen(2),pLen(2),rList(2+(k-1)*3,j),2}; % for R - only rad & cont
                               % R i.e. kGabor1 corresponds to row number 2,5,8,..
                               % and so on. Hence it is 2+(k-1)*3

                               % [Vinay] -  C is matched to S and therefore the entries
                               % for gabor3 i.e. can be ignored. Otherwise this can 
                               % give spurious results if the param values for C are 
                               % not set equal to those for S. Hence pos3 is the full
                               % set here
                               pos3 = parameterCombinations{aLen(3),eLen(3),sLen(3),fLen(3),oLen(3),cLen(3),tLen(3),pLen(3),rList(3+(k-1)*3,j),3}; % for C
                               % also keeping r variable here, just in case you
                               % have multiple r values
                               
                               tagPos1 = [num2str(aList(1+(k-1)*3,j)) num2str(eList(1+(k-1)*3,j)) num2str(sList(1+(k-1)*3,j))...
                                   num2str(fList(1+(k-1)*3,j)) num2str(oList(1+(k-1)*3,j)) num2str(cList(1+(k-1)*3,j))...
                                   num2str(tList(1+(k-1)*3,j)) num2str(pList(1+(k-1)*3,j)) num2str(rList(1+(k-1)*3,j))];
                               tagPos2 = [num2str(aLen(2)) num2str(eLen(2)) num2str(sLen(2)) num2str(fLen(2))...
                                   num2str(oLen(2)) num2str(cList(2+(k-1)*3,j)) num2str(tLen(2)) num2str(pLen(2)) num2str(rList(2+(k-1)*3,j))];
                               tagPos3 = [num2str(aLen(3)) num2str(eLen(3)) num2str(sLen(3)) num2str(fLen(3))...
                                   num2str(oLen(3)) num2str(cLen(3)) num2str(tLen(3)) num2str(pLen(3)) num2str(rList(3+(k-1)*3,j))];

                           case 3 % Dual Contrast Protocol
                               pos1 = parameterCombinations{aLen(1),eLen(1),sLen(1),fLen(1),oLen(1),cLen(1),tLen(1),pLen(1),rLen(1),1}; % for C
                               % [Vinay] - S is hidden in this case and only C & R
                               % define the stimuli. pos1 is therefore the full set
                               % here
                               pos2 = parameterCombinations{aLen(2),eLen(2),sLen(2),fLen(2),oLen(2),cList(2+(k-1)*3,j),tLen(2),pLen(2),rList(2+(k-1)*3,j),2}; % for R - only rad & cont
                               % [Vinay] -  R params - cont & rad. Rest are matched to
                               % C's params
                               pos3 = parameterCombinations{aList(3+(k-1)*3,j),eList(3+(k-1)*3,j),sList(3+(k-1)*3,j),fList(3+(k-1)*3,j),oList(3+(k-1)*3,j),cList(3+(k-1)*3,j),tList(3+(k-1)*3,j),pList(3+(k-1)*3,j),rList(3+(k-1)*3,j),3};
                               % for C. C is the defining gabor here
                               
                               tagPos1 = [num2str(aLen(1)) num2str(eLen(1)) num2str(sLen(1))...
                                   num2str(fLen(1)) num2str(oLen(1)) num2str(cLen(1))...
                                   num2str(tLen(1)) num2str(pLen(1)) num2str(rLen(1))];
                               tagPos2 = [num2str(aLen(2)) num2str(eLen(2)) num2str(sLen(2)) num2str(fLen(2))...
                                   num2str(oLen(2)) num2str(cList(2+(k-1)*3,j)) num2str(tLen(2)) num2str(pLen(2)) num2str(rList(2+(k-1)*3,j))];
                               tagPos3 = [num2str(aList(3+(k-1)*3,j)) num2str(eList(3+(k-1)*3,j)) num2str(sList(3+(k-1)*3,j))...
                                   num2str(fList(3+(k-1)*3,j)) num2str(oList(3+(k-1)*3,j)) num2str(cList(3+(k-1)*3,j))...
                                   num2str(tList(3+(k-1)*3,j)) num2str(pList(3+(k-1)*3,j)) num2str(rList(3+(k-1)*3,j))];
                               

                           case 4 % Dual Orientation Protocol
                               pos1 = parameterCombinations{aLen(1),eLen(1),sLen(1),fLen(1),oLen(1),cLen(1),tLen(1),pLen(1),rLen(1),1}; % for C
                               % [Vinay] - S is hidden in this case and only C & R
                               % define the stimuli. pos1 is therefore the full set
                               % here
                               pos2 = parameterCombinations{aLen(2),eLen(2),sLen(2),fLen(2),oList(2+(k-1)*3,j),cLen(2),tLen(2),pLen(2),rList(2+(k-1)*3,j),2}; % for R - only rad & ori
                               % [Vinay] -  R params - ori & rad. Rest are matched to
                               % C's params
                               pos3 = parameterCombinations{aList(3+(k-1)*3,j),eList(3+(k-1)*3,j),sList(3+(k-1)*3,j),fList(3+(k-1)*3,j),oList(3+(k-1)*3,j),cList(3+(k-1)*3,j),tList(3+(k-1)*3,j),pList(3+(k-1)*3,j),rList(3+(k-1)*3,j),3};
                               % for C. C is the defining gabor here
                               
                               tagPos1 = [num2str(aLen(1)) num2str(eLen(1)) num2str(sLen(1))...
                                   num2str(fLen(1)) num2str(oLen(1)) num2str(cLen(1))...
                                   num2str(tLen(1)) num2str(pLen(1)) num2str(rLen(1))];
                               tagPos2 = [num2str(aLen(2)) num2str(eLen(2)) num2str(sLen(2)) num2str(fLen(2))...
                                   num2str(oList(2+(k-1)*3,j)) num2str(cLen(2)) num2str(tLen(2)) num2str(pLen(2)) num2str(rList(2+(k-1)*3,j))];
                               tagPos3 = [num2str(aList(3+(k-1)*3,j)) num2str(eList(3+(k-1)*3,j)) num2str(sList(3+(k-1)*3,j))...
                                   num2str(fList(3+(k-1)*3,j)) num2str(oList(3+(k-1)*3,j)) num2str(cList(3+(k-1)*3,j))...
                                   num2str(tList(3+(k-1)*3,j)) num2str(pList(3+(k-1)*3,j)) num2str(rList(3+(k-1)*3,j))];

                           case 5 % Dual Phase Protocol
                               pos1 = parameterCombinations{aLen(1),eLen(1),sLen(1),fLen(1),oLen(1),cLen(1),tLen(1),pLen(1),rLen(1),1}; % for C
                               % [Vinay] - S is hidden in this case and only C & R
                               % define the stimuli. pos1 is therefore the full set
                               % here
                               pos2 = parameterCombinations{aLen(2),eLen(2),sLen(2),fLen(2),oLen(2),cLen(2),tLen(2),pList(2+(k-1)*3,j),rList(2+(k-1)*3,j),2}; % for R - only rad & p
                               % [Vinay] -  R params - p & rad. Rest are matched to
                               % C's params
                               pos3 = parameterCombinations{aList(3+(k-1)*3,j),eList(3+(k-1)*3,j),sList(3+(k-1)*3,j),fList(3+(k-1)*3,j),oList(3+(k-1)*3,j),cList(3+(k-1)*3,j),tList(3+(k-1)*3,j),pList(3+(k-1)*3,j),rList(3+(k-1)*3,j),3};
                               % for C. C is the defining gabor here
                               
                               tagPos1 = [num2str(aLen(1)) num2str(eLen(1)) num2str(sLen(1))...
                                   num2str(fLen(1)) num2str(oLen(1)) num2str(cLen(1))...
                                   num2str(tLen(1)) num2str(pLen(1)) num2str(rLen(1))];
                               tagPos2 = [num2str(aLen(2)) num2str(eLen(2)) num2str(sLen(2)) num2str(fLen(2))...
                                   num2str(oLen(2)) num2str(cLen(2)) num2str(tLen(2)) num2str(pList(2+(k-1)*3,j)) num2str(rList(2+(k-1)*3,j))];
                               tagPos3 = [num2str(aList(3+(k-1)*3,j)) num2str(eList(3+(k-1)*3,j)) num2str(sList(3+(k-1)*3,j))...
                                   num2str(fList(3+(k-1)*3,j)) num2str(oList(3+(k-1)*3,j)) num2str(cList(3+(k-1)*3,j))...
                                   num2str(tList(3+(k-1)*3,j)) num2str(pList(3+(k-1)*3,j)) num2str(rList(3+(k-1)*3,j))];

                           case 6 %Orientation Ring Protocol
                               pos1 = parameterCombinations{aList(1+(k-1)*3,j),eList(1+(k-1)*3,j),sList(1+(k-1)*3,j),fList(1+(k-1)*3,j),oList(1+(k-1)*3,j),cList(1+(k-1)*3,j),tList(1+(k-1)*3,j),pList(1+(k-1)*3,j),rList(1+(k-1)*3,j),1};
                               % S params define C & S
                               % S i.e. kGabor0 corresponds to row number 1,4,7,..
                               % and so on. Hence it is 1+(k-1)*3
                               pos2 = parameterCombinations{aLen(2),eLen(2),sLen(2),fLen(2),oList(2+(k-1)*3,j),cLen(2),tLen(2),pLen(2),rList(2+(k-1)*3,j),2}; % for R - only rad & ori
                               % [Vinay] -  C is matched to S and therefore the entries
                               % for gabor3 i.e. can be ignored. Otherwise this can 
                               % give spurious results if the param values for C are 
                               % not set equal to those for S. Hence pos3 is the full
                               % set here
                               pos3 = parameterCombinations{aLen(3),eLen(3),sLen(3),fLen(3),oLen(3),cLen(3),tLen(3),pLen(3),rList(3+(k-1)*3,j),3}; % for C
                               % also keeping r variable here, just in case you
                               % have multiple r values
                               
                               tagPos1 = [num2str(aList(1+(k-1)*3,j)) num2str(eList(1+(k-1)*3,j)) num2str(sList(1+(k-1)*3,j))...
                                   num2str(fList(1+(k-1)*3,j)) num2str(oList(1+(k-1)*3,j)) num2str(cList(1+(k-1)*3,j))...
                                   num2str(tList(1+(k-1)*3,j)) num2str(pList(1+(k-1)*3,j)) num2str(rList(1+(k-1)*3,j))];
                               tagPos2 = [num2str(aLen(2)) num2str(eLen(2)) num2str(sLen(2)) num2str(fLen(2))...
                                   num2str(oList(2+(k-1)*3,j)) num2str(cLen(2)) num2str(tLen(2)) num2str(pLen(2)) num2str(rList(2+(k-1)*3,j))];
                               tagPos3 = [num2str(aLen(3)) num2str(eLen(3)) num2str(sLen(3)) num2str(fLen(3))...
                                   num2str(oLen(3)) num2str(cLen(3)) num2str(tLen(3)) num2str(pLen(3)) num2str(rList(3+(k-1)*3,j))];

                           case 7 %Phase Ring Protocol
                               pos1 = parameterCombinations{aList(1+(k-1)*3,j),eList(1+(k-1)*3,j),sList(1+(k-1)*3,j),fList(1+(k-1)*3,j),oList(1+(k-1)*3,j),cList(1+(k-1)*3,j),tList(1+(k-1)*3,j),pList(1+(k-1)*3,j),rList(1+(k-1)*3,j),1};
                               % S params define C & S
                               % S i.e. kGabor0 corresponds to row number 1,4,7,..
                               % and so on. Hence it is 1+(k-1)*3
                               pos2 = parameterCombinations{aLen(2),eLen(2),sLen(2),fLen(2),oLen(2),cLen(2),tLen(2),pList(2+(k-1)*3,j),rList(2+(k-1)*3,j),2}; % for R - only rad & p
                               % [Vinay] -  C is matched to S and therefore the entries
                               % for gabor3 i.e. can be ignored. Otherwise this can 
                               % give spurious results if the param values for C are 
                               % not set equal to those for S. Hence pos3 is the full
                               % set here
                               pos3 = parameterCombinations{aLen(3),eLen(3),sLen(3),fLen(3),oLen(3),cLen(3),tLen(3),pLen(3),rList(3+(k-1)*3,j),3}; % for C
                               % also keeping r variable here, just in case you
                               % have multiple r values
                               
                               tagPos1 = [num2str(aList(1+(k-1)*3,j)) num2str(eList(1+(k-1)*3,j)) num2str(sList(1+(k-1)*3,j))...
                                   num2str(fList(1+(k-1)*3,j)) num2str(oList(1+(k-1)*3,j)) num2str(cList(1+(k-1)*3,j))...
                                   num2str(tList(1+(k-1)*3,j)) num2str(pList(1+(k-1)*3,j)) num2str(rList(1+(k-1)*3,j))];
                               tagPos2 = [num2str(aLen(2)) num2str(eLen(2)) num2str(sLen(2)) num2str(fLen(2))...
                                   num2str(oLen(2)) num2str(cLen(2)) num2str(tLen(2)) num2str(pList(2+(k-1)*3,j)) num2str(rList(2+(k-1)*3,j))];
                               tagPos3 = [num2str(aLen(3)) num2str(eLen(3)) num2str(sLen(3)) num2str(fLen(3))...
                                   num2str(oLen(3)) num2str(cLen(3)) num2str(tLen(3)) num2str(pLen(3)) num2str(rList(3+(k-1)*3,j))];

                           case 8 %Drifting Ring Protocol
                               pos1 = parameterCombinations{aList(1+(k-1)*3,j),eList(1+(k-1)*3,j),sList(1+(k-1)*3,j),fList(1+(k-1)*3,j),oList(1+(k-1)*3,j),cList(1+(k-1)*3,j),tList(1+(k-1)*3,j),pList(1+(k-1)*3,j),rList(1+(k-1)*3,j),1};
                               % S params define C & S
                               % S i.e. kGabor0 corresponds to row number 1,4,7,..
                               % and so on. Hence it is 1+(k-1)*3
                               pos2 = parameterCombinations{aLen(2),eLen(2),sLen(2),fLen(2),oLen(2),cLen(2),tList(2+(k-1)*3,j),pLen(2),rList(2+(k-1)*3,j),2}; % for R - only rad & t
                               % [Vinay] -  C is matched to S and therefore the entries
                               % for gabor3 i.e. can be ignored. Otherwise this can 
                               % give spurious results if the param values for C are 
                               % not set equal to those for S. Hence pos3 is the full
                               % set here
                               pos3 = parameterCombinations{aLen(3),eLen(3),sLen(3),fLen(3),oLen(3),cLen(3),tLen(3),pLen(3),rList(3+(k-1)*3,j),3}; % for C
                               % also keeping r variable here, just in case you
                               % have multiple r values
                               
                               tagPos1 = [num2str(aList(1+(k-1)*3,j)) num2str(eList(1+(k-1)*3,j)) num2str(sList(1+(k-1)*3,j))...
                                   num2str(fList(1+(k-1)*3,j)) num2str(oList(1+(k-1)*3,j)) num2str(cList(1+(k-1)*3,j))...
                                   num2str(tList(1+(k-1)*3,j)) num2str(pList(1+(k-1)*3,j)) num2str(rList(1+(k-1)*3,j))];
                               tagPos2 = [num2str(aLen(2)) num2str(eLen(2)) num2str(sLen(2)) num2str(fLen(2))...
                                   num2str(oLen(2)) num2str(cLen(2)) num2str(tList(2+(k-1)*3,j)) num2str(pLen(2)) num2str(rList(2+(k-1)*3,j))];
                               tagPos3 = [num2str(aLen(3)) num2str(eLen(3)) num2str(sLen(3)) num2str(fLen(3))...
                                   num2str(oLen(3)) num2str(cLen(3)) num2str(tLen(3)) num2str(pLen(3)) num2str(rList(3+(k-1)*3,j))];

                           case 9 %Cross Orientation Protocol
                               pos1 = parameterCombinations{aLen(1),eLen(1),sLen(1),fLen(1),oLen(1),cLen(1),tLen(1),pLen(1),rLen(1),1}; % for C
                               % [Vinay] - S is hidden in this case and only C & R
                               % define the stimuli. pos1 is therefore the full set
                               % here
                               pos2 = parameterCombinations{aList(2+(k-1)*3,j),eList(2+(k-1)*3,j),sList(2+(k-1)*3,j),fList(2+(k-1)*3,j),oList(2+(k-1)*3,j),cList(2+(k-1)*3,j),tList(2+(k-1)*3,j),pList(2+(k-1)*3,j),rList(2+(k-1)*3,j),2};
                               % for R. R is one of the defining gabor here
                               pos3 = parameterCombinations{aList(3+(k-1)*3,j),eList(3+(k-1)*3,j),sList(3+(k-1)*3,j),fList(3+(k-1)*3,j),oList(3+(k-1)*3,j),cList(3+(k-1)*3,j),tList(3+(k-1)*3,j),pList(3+(k-1)*3,j),rList(3+(k-1)*3,j),3};
                               % for C. C is one of the defining gabor here
                               
                               tagPos1 = [num2str(aLen(1)) num2str(eLen(1)) num2str(sLen(1))...
                                   num2str(fLen(1)) num2str(oLen(1)) num2str(cLen(1))...
                                   num2str(tLen(1)) num2str(pLen(1)) num2str(rLen(1))];
                               tagPos2 = [num2str(aList(2+(k-1)*3,j)) num2str(eList(2+(k-1)*3,j)) num2str(sList(2+(k-1)*3,j))...
                                   num2str(fList(2+(k-1)*3,j)) num2str(oList(2+(k-1)*3,j)) num2str(cList(2+(k-1)*3,j))...
                                   num2str(tList(2+(k-1)*3,j)) num2str(pList(2+(k-1)*3,j)) num2str(rList(2+(k-1)*3,j))];
                               tagPos3 = [num2str(aList(3+(k-1)*3,j)) num2str(eList(3+(k-1)*3,j)) num2str(sList(3+(k-1)*3,j))...
                                   num2str(fList(3+(k-1)*3,j)) num2str(oList(3+(k-1)*3,j)) num2str(cList(3+(k-1)*3,j))...
                                   num2str(tList(3+(k-1)*3,j)) num2str(pList(3+(k-1)*3,j)) num2str(rList(3+(k-1)*3,j))];

                           case 10 %Annulus Fixed Protocol
                               pos1 = parameterCombinations{aList(1+(k-1)*3,j),eList(1+(k-1)*3,j),sList(1+(k-1)*3,j),fList(1+(k-1)*3,j),oList(1+(k-1)*3,j),cList(1+(k-1)*3,j),tList(1+(k-1)*3,j),pList(1+(k-1)*3,j),rList(1+(k-1)*3,j),1};
                               % S params define C & S
                               % S i.e. kGabor0 corresponds to row number 1,4,7,..
                               % and so on. Hence it is 1+(k-1)*3
                               pos2 = parameterCombinations{aLen(2),eLen(2),sLen(2),fLen(2),oLen(2),cList(2+(k-1)*3,j),tLen(2),pLen(2),rList(2+(k-1)*3,j),2}; % for R - only rad & cont
                               % [Vinay] -  C is matched to S and therefore the entries
                               % for gabor3 i.e. can be ignored. Otherwise this can 
                               % give spurious results if the param values for C are 
                               % not set equal to those for S. However, the rad of C is
                               % still a relevant parameter and has to be taken care of
                               pos3 = parameterCombinations{aLen(3),eLen(3),sLen(3),fLen(3),oLen(3),cLen(3),tLen(3),pLen(3),rList(3+(k-1)*3,j),3}; % for C - only rad
                               
                               tagPos1 = [num2str(aList(1+(k-1)*3,j)) num2str(eList(1+(k-1)*3,j)) num2str(sList(1+(k-1)*3,j))...
                                   num2str(fList(1+(k-1)*3,j)) num2str(oList(1+(k-1)*3,j)) num2str(cList(1+(k-1)*3,j))...
                                   num2str(tList(1+(k-1)*3,j)) num2str(pList(1+(k-1)*3,j)) num2str(rList(1+(k-1)*3,j))];
                               tagPos2 = [num2str(aLen(2)) num2str(eLen(2)) num2str(sLen(2)) num2str(fLen(2))...
                                   num2str(oLen(2)) num2str(cList(2+(k-1)*3,j)) num2str(tLen(2)) num2str(pLen(2)) num2str(rList(2+(k-1)*3,j))];
                               tagPos3 = [num2str(aLen(3)) num2str(eLen(3)) num2str(sLen(3)) num2str(fLen(3))...
                                   num2str(oLen(3)) num2str(cLen(3)) num2str(tLen(3)) num2str(pLen(3)) num2str(rList(3+(k-1)*3,j))];

                    otherwise
                        pos1 = parameterCombinations{aList(1+(k-1)*3,j),eList(1+(k-1)*3,j),sList(1+(k-1)*3,j),fList(1+(k-1)*3,j),oList(1+(k-1)*3,j),cList(1+(k-1)*3,j),tList(1+(k-1)*3,j),pList(1+(k-1)*3,j),rList(1+(k-1)*3,j),1};
                        % good positions for gabor 1 i.e. S
                        pos2 = parameterCombinations{aList(2+(k-1)*3,j),eList(2+(k-1)*3,j),sList(2+(k-1)*3,j),fList(2+(k-1)*3,j),oList(2+(k-1)*3,j),cList(2+(k-1)*3,j),tList(2+(k-1)*3,j),pList(2+(k-1)*3,j),rList(2+(k-1)*3,j),2};
                        pos3 = parameterCombinations{aList(3+(k-1)*3,j),eList(3+(k-1)*3,j),sList(3+(k-1)*3,j),fList(3+(k-1)*3,j),oList(3+(k-1)*3,j),cList(3+(k-1)*3,j),tList(3+(k-1)*3,j),pList(3+(k-1)*3,j),rList(3+(k-1)*3,j),3};
                        
                           tagPos1 = [num2str(aList(1+(k-1)*3,j)) num2str(eList(1+(k-1)*3,j)) num2str(sList(1+(k-1)*3,j))...
                               num2str(fList(1+(k-1)*3,j)) num2str(oList(1+(k-1)*3,j)) num2str(cList(1+(k-1)*3,j))...
                               num2str(tList(1+(k-1)*3,j)) num2str(pList(1+(k-1)*3,j)) num2str(rList(1+(k-1)*3,j))];
                           tagPos2 = [num2str(aList(2+(k-1)*3,j)) num2str(eList(2+(k-1)*3,j)) num2str(sList(2+(k-1)*3,j))...
                               num2str(fList(2+(k-1)*3,j)) num2str(oList(2+(k-1)*3,j)) num2str(cList(2+(k-1)*3,j))...
                               num2str(tList(2+(k-1)*3,j)) num2str(pList(2+(k-1)*3,j)) num2str(rList(2+(k-1)*3,j))];
                           tagPos3 = [num2str(aList(3+(k-1)*3,j)) num2str(eList(3+(k-1)*3,j)) num2str(sList(3+(k-1)*3,j))...
                               num2str(fList(3+(k-1)*3,j)) num2str(oList(3+(k-1)*3,j)) num2str(cList(3+(k-1)*3,j))...
                               num2str(tList(3+(k-1)*3,j)) num2str(pList(3+(k-1)*3,j)) num2str(rList(3+(k-1)*3,j))];
                end


                    goodPos = intersect(pos1,pos2);
                    goodPos = intersect(goodPos,pos3);
        
                end
            
            %------------------
            % Vinay - select good trials as per the electrode(s)
    %         useAllBadTrials = 1;
            if refType ~= 2 || ~existsBipolarData % not bipolar -> electrodes remain the same
                
                if useAllBadTrials && existsBadTrialFile
                    elecIndex1 = (strcmp(channelString,nameElec) == 1);

                    if ~useDiff
                        elecBadTrials = allBadTrials{elecIndex1};
                        disp(['No. of Bad trials for ' channelString ': ' num2str(length(elecBadTrials))]);
                    else
                        elecIndex2 = (strcmp(analogChannelString2,nameElec) == 1);
                        elecBadTrials = unique([allBadTrials{elecIndex1} allBadTrials{elecIndex2}]);
                        disp(['No. of Bad trials for ' channelString ' and ' analogChannelString2 ': ' num2str(length(elecBadTrials))]);
                    end
                    
                    goodPos = setdiff(goodPos,elecBadTrials);
                    
                else

                    goodPos = setdiff(goodPos,badTrials);
                end
            else
                
                if useAllBadTrials && existsBadTrialFile
                    if ~useDiff 
                        elecBadTrials = getBipolarBadTrials(channelString,allBadTrials);
                        disp(['No. of Bad trials for ' channelString ': ' num2str(length(elecBadTrials))]);
                    else
                        elecBadTrials1 = getBipolarBadTrials(channelString,allBadTrials);
                        elecBadTrials2 = getBipolarBadTrials(analogChannelString2,allBadTrials);
                        elecBadTrials = unique([elecBadTrials1 elecBadTrials2]);
                        disp(['No. of Bad trials for ' channelString ' and ' analogChannelString2 ': ' num2str(length(elecBadTrials))]);
                    end
                    goodPos = setdiff(goodPos,elecBadTrials);
                    
                else
                    goodPos = setdiff(goodPos,badTrials);
                end
            end
            %-------------------
            
            
            if refType == 2
                % Determine combination schemes
                % none, all, occ, par, temp, front
                switch refCombineType
                    case 2
                        combineElec = (1:numElectrodes);
                    case 3
                        if contralateralNeighbour
                            combineElec = [910,4546,5960,6364]; % occ
                        else % nearest neighbour
                            combineElec = [910,4559,4660,6364]; % occ
                        end
                    case 4
                        if contralateralNeighbour
                            combineElec = [78,1516,3738,5152,1963]; % par
                        else % nearest neighbour
                            combineElec = [1551,737,838,1652,1963]; % par
                        end
                    case 5
                        if contralateralNeighbour
                            combineElec = [1314,3132,2930,5556,5758]; % tempo
                        else % nearest neighbour
                            combineElec = [1349,1450,3157,3258,2955,3056]; % tempo
                        end
                    case 6
                        if contralateralNeighbour
                            combineElec = [12,3940,5354,1112,4748,34,3334,1761,2526,4142,2122]; % front
                        else % nearest neighbour
                            combineElec = [12,1761,3953,4054,1147,333,1248,434,2541,2122,2642]; % front
                        end
                    case 7
                        if contralateralNeighbour
                            combineElec = [549,1862,3536,650,2743,2324,2844]; % central
                        else % nearest neighbour
                            combineElec = [56,1862,3536,4950,2728,2324,4344]; % central
                        end
                end

                switch refDiffType
                    case 2
                        diffElec = (1:numElectrodes);
                    case 3
                        if contralateralNeighbour
                            diffElec = [910,4546,5960,6364]; % occ
                        else % nearest neighbour
                            diffElec = [910,4559,4660,6364]; % occ
                        end
                    case 4
                        if contralateralNeighbour
                            diffElec = [78,1516,3738,5152,1963]; % par
                        else % nearest neighbour
                            diffElec = [1551,737,838,1652,1963]; % par
                        end
                    case 5
                        if contralateralNeighbour
                            diffElec = [1314,3132,2930,5556,5758]; % tempo
                        else % nearest neighbour
                            diffElec = [1349,1450,3157,3258,2955,3056]; % tempo
                        end
                    case 6
                        if contralateralNeighbour
                            diffElec = [12,3940,5354,1112,4748,34,3334,1761,2526,4142,2122]; % front
                        else % nearest neighbour
                            diffElec = [12,1761,3953,4054,1147,333,1248,434,2541,2122,2642]; % front
                        end
                    case 7
                        if contralateralNeighbour
                            diffElec = [549,1862,3536,650,2743,2324,2844]; % central
                        else % nearest neighbour
                            diffElec = [56,1862,3536,4950,2728,2324,4344]; % central
                        end
                end
                
            else
                % Determine combination schemes
                % none, all, occ, par, temp, front, central
                switch refCombineType
                    case 2
                        combineElec = (1:numElectrodes);
                    case 3
                        combineElec = [9,10,45,46,59,60,63,64]; % occ
                    case 4
                        combineElec = [15,51,7,37,19,38,8,52,16]; % par
                    case 5
                        combineElec = [29,55,13,57,31,30,56,14,58,32]; % tempo
                    case 6
                        combineElec = [1,2,61,53,39,40,54,11,47,3,33,17,34,4,48,12,25,41,21,22,42,26]; % front
                    case 7
                        combineElec = [49,5,35,18,36,6,50,27,43,23,62,24,44,28]; % centre
                end

                switch refDiffType
                    case 2
                        diffElec = (1:numElectrodes);
                    case 3
                        diffElec = [9,10,45,46,59,60,63,64]; % occ
                    case 4
                        diffElec = [15,51,7,37,19,38,8,52,16]; % par
                    case 5
                        diffElec = [29,55,13,57,31,30,56,14,58,32]; % tempo
                    case 6
                        diffElec = [1,2,61,53,39,40,54,11,47,3,33,17,34,4,48,12,25,41,21,22,42,26]; % front
                    case 7
                        diffElec = [49,5,35,18,36,6,50,27,43,23,62,24,44,28]; % centre
                end
            end
            
            % Determine goodPos for the combination cases
            
            if refCombineType ~= 1
                numCombineElec = length(combineElec);
                
                elecBadTrials = []; elecBadTrials1 = [];
                for c = 1:numCombineElec
                    combineElecString = ['elec' num2str(combineElec(c))];
                    if refType == 2 && existsBipolarData
                        elecBadTrials1{c} = getBipolarBadTrials(combineElecString,allBadTrials);
                    elseif refType == 2 && ~existsBipolarData
                        combineElecStr = num2str(combineElec(c));
                        if length(combineElecStr) == 2
                            combineElecString1 = ['elec' combineElecStr(1)];
                            combineElecString2 = ['elec' combineElecStr(2)];
                        elseif length(combineElecStr) == 3
                            combineElecString1 = ['elec' combineElecStr(1)];
                            combineElecString2 = ['elec' combineElecStr(2:3)];
                        elseif length(combineElecStr) == 4
                            combineElecString1 = ['elec' combineElecStr(1:2)];
                            combineElecString2 = ['elec' combineElecStr(3:4)];
                        end
                        elecIndex1 = (strcmp(combineElecString1,nameElec) == 1);
                        elecBadTrials1{c} = allBadTrials{elecIndex1};
                        
                        elecIndex2 = (strcmp(combineElecString2,nameElec) == 1);
                        elecBadTrials2{c} = allBadTrials{elecIndex2};
                        
                        elecBadTrials1{c} = unique([elecBadTrials1{c} elecBadTrials2{c}]);
                        
                    else
                        elecIndex1 = (strcmp(combineElecString,nameElec) == 1);
                        elecBadTrials1{c} = allBadTrials{elecIndex1}; 
                    end
                    disp('COMBINATION ELECTRODES');
                    disp(['No. of Bad trials for ' combineElecString ': ' num2str(length(elecBadTrials1{c}))]);
                    elecBadTrials = unique([elecBadTrials elecBadTrials1{c}]);
                end
            end
            
            % Determine goodPos for the difference cases
            
            if refDiffType ~= 1
                numDiffElec = length(diffElec);
                
                elecBadTrialsD = []; elecBadTrials1D = [];
                for c = 1:numDiffElec
                    diffElecString = ['elec' num2str(diffElec(c))];
                    if refType == 2 && existsBipolarData
                        elecBadTrials1D{c} = getBipolarBadTrials(diffElecString,allBadTrials);
                    elseif refType == 2 && ~existsBipolarData
                        diffElecStr = num2str(diffElec(c));
                        if length(diffElecStr) == 2
                            diffElecString1 = ['elec' diffElecStr(1)];
                            diffElecString2 = ['elec' diffElecStr(2)];
                        elseif length(diffElecStr) == 3
                            diffElecString1 = ['elec' diffElecStr(1)];
                            diffElecString2 = ['elec' diffElecStr(2:3)];
                        elseif length(diffElecStr) == 4
                            diffElecString1 = ['elec' diffElecStr(1:2)];
                            diffElecString2 = ['elec' diffElecStr(3:4)];
                        end
                        elecIndex1 = (strcmp(diffElecString1,nameElec) == 1);
                        elecBadTrials1{c} = allBadTrials{elecIndex1};
                        
                        elecIndex2 = (strcmp(diffElecString2,nameElec) == 1);
                        elecBadTrials2{c} = allBadTrials{elecIndex2};
                        
                        elecBadTrials1D{c} = unique([elecBadTrials1{c} elecBadTrials2{c}]);
                    else
                        elecIndex1 = (strcmp(diffElecString,nameElec) == 1);
                        elecBadTrials1D{c} = allBadTrials{elecIndex1}; 
                    end
                    disp('DIFFERENCE ELECTRODES');
                    disp(['No. of Bad trials for ' diffElecString ': ' num2str(length(elecBadTrials1D{c}))]);
                    elecBadTrialsD = unique([elecBadTrialsD elecBadTrials1D{c}]);
                end
            end
            
            
            %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            % Main computation loop
            %--------------------------------------------------------------

            if isempty(goodPos)
                disp('No entries for this combination..')
            else
                disp(['pos=' num2str(k) ',' num2str(j) ',n=' num2str(length(goodPos))]);

                Fs = round(1/(timeVals(2)-timeVals(1)));
                BLRange = uint16((BLMax-BLMin)*Fs);
                STRange = uint16((STMax-STMin)*Fs);
                BLPos = find(timeVals>=BLMin,1)+ (1:BLRange);
                STPos = find(timeVals>=STMin,1)+ (1:STRange);

                xsBL = 0:1/(BLMax-BLMin):Fs-1/(BLMax-BLMin);
                xsST = 0:1/(STMax-STMin):Fs-1/(STMax-STMin);

            % Get the data
            clear signal analogData analogDataNotched
            if (refCombineType == 1 && refDiffType == 1) % no combination or difference
                switch refType
                    case 1
                        load(fullfile(folderLFP,channelString));
                    case 2
                        bipolarSignalFile = fullfile(folderBipolar,channelString);
                        if exist(bipolarSignalFile,'file')
                            load(bipolarSignalFile);
                        else
                            analogChannelString1 = channelString;
                            load(fullfile(folderLFP,analogChannelString1));
                            if notchData
                                analogData1 = analogDataNotched;
                            else
                                analogData1 = analogData; % Vinay - these are the 
                                % freshly loaded values for electrode 1
                            end
                            
                            clear analogData analogDataNotched
                            load(fullfile(folderAverage,analogChannelString2));
                            if notchData
                                analogData2 = analogDataNotched;
                            else
                                analogData2 = analogData; % Vinay - these are the 
                                % freshly loaded values for electrode 2
                            end

                            analogData = analogData1 - analogData2;
                        end
                            
                    case 3
                        load(fullfile(folderAverage,channelString));
                        aD1 = analogData;
                        if notchData
                            aDN1 = analogDataNotched;
                        end
                        clear analogData analogDataNotched
                        load(fullfile(folderAverage,'avgRef.mat'));
                        analogData = aD1 - analogData;
                        if notchData
                            analogDataNotched = aDN1 - analogDataNotched;
                        end
                    case 4
                        load(fullfile(folderCSD,channelString));
                    otherwise
                        load(fullfile(folderLFP,channelString));
                end

                % Vinay - added this notch data check
                if notchData
                    analogData = analogDataNotched;
                end

                if useDiff
                    analogChannelString1 = channelString; % first electrode selected
                    clear analogData analogDataNotched
                    switch refType
                        case 1
                            load(fullfile(folderLFP,analogChannelString1));
                        case 2
                            if existsBipolarData
                                load(fullfile(folderBipolar,analogChannelString1));
                            else
                                load(fullfile(folderLFP,analogChannelString1));
                            end
                        case 3
                            load(fullfile(folderAverage,analogChannelString1));
                            aD1 = analogData;
                            if notchData
                                aDN1 = analogDataNotched;
                            end
                            clear analogData analogDataNotched
                            load(fullfile(folderAverage,'avgRef.mat'));
                            analogData = aD1 - analogData;
                            if notchData
                                analogDataNotched = aDN1 - analogDataNotched;
                            end
                        case 4
                            load(fullfile(folderCSD,analogChannelString1));
                        otherwise
                            load(fullfile(folderLFP,analogChannelString1));
                    end
                
                    if notchData
                        analogData1 = analogDataNotched;
                    else
                        analogData1 = analogData; % Vinay - these are the 
                        % freshly loaded values for electrode 1
                    end
                    % Vinay - store the analogData 
                    % for electrode 1 in a separate variable 
                    % otherwise loading the data for 2nd electrode will
                    % overwrite analogData
                    clear analogData analogDataNotched
                    switch refType
                        case 1
                            load(fullfile(folderLFP,analogChannelString2));
                        case 2
                            if existsBipolarData
                                load(fullfile(folderBipolar,analogChannelString2));
                            else
                                load(fullfile(folderLFP,analogChannelString2));
                            end
                        case 3
                            load(fullfile(folderAverage,analogChannelString2));
                            aD1 = analogData;
                            if notchData
                                aDN1 = analogDataNotched;
                            end
                            clear analogData analogDataNotched
                            load(fullfile(folderAverage,'avgRef.mat'));
                            analogData = aD1 - analogData;
                            if notchData
                                analogDataNotched = aDN1 - analogDataNotched;
                            end
                        case 4
                            load(fullfile(folderCSD,analogChannelString2));
                        otherwise
                            load(fullfile(folderLFP,analogChannelString2));
                    end
                    
                    if notchData
                        analogData2 = analogDataNotched;
                    else
                        analogData2 = analogData; % Vinay - these are the 
                        % freshly loaded values for electrode 2
                    end
                    
                    
                    analogData = analogData1 - analogData2;
                end
                
                analogData = analogData(goodPos,:); % gives the relevant signal
            end
            
            % For the combination cases, directly extract the goodPos part
            % of analogData (for each electrode) and combine those.
            % Post-extaction isn't as straightforward as in the isolated
            % case (used above)
            if refCombineType ~= 1
                numCombineElec = length(combineElec);
                totalNumStim = 0;
                for c = 1:numCombineElec
                    clear analogData analogDataNotched
                    combineElecString = ['elec' num2str(combineElec(c)) '.mat'];
                    switch refType
                    case 1
                        load(fullfile(folderLFP,combineElecString));                  
                    case 2
                        combineElecStr = num2str(combineElec(c));
                        if length(combineElecStr) == 2
                            combineElecString1 = ['elec' combineElecStr(1) '.mat'];
                            combineElecString2 = ['elec' combineElecStr(2) '.mat'];
                        elseif length(combineElecStr) == 3
                            combineElecString1 = ['elec' combineElecStr(1) '.mat'];
                            combineElecString2 = ['elec' combineElecStr(2:3) '.mat'];
                        elseif length(combineElecStr) == 4
                            combineElecString1 = ['elec' combineElecStr(1:2) '.mat'];
                            combineElecString2 = ['elec' combineElecStr(3:4) '.mat'];
                        end
                        clear analogData analogDataNotched
                        load(fullfile(folderLFP,combineElecString1));
                        aD1 = analogData;
                        if notchData
                            aDN1 = analogDataNotched;
                        end
                        clear analogData analogDataNotched
                        load(fullfile(folderLFP,combineElecString2));
                        analogData = aD1 - analogData;
                        if notchData
                            analogDataNotched = aDN1 - analogDataNotched;
                        end
                    case 3
                        load(fullfile(folderAverage,combineElecString));
                        aD1 = analogData;
                        if notchData
                            aDN1 = analogDataNotched;
                        end
                        clear analogData analogDataNotched
                        load(fullfile(folderAverage,'avgRef.mat'));
                        analogData = aD1 - analogData;
                        if notchData
                            analogDataNotched = aDN1 - analogDataNotched;
                        end
                    case 4
                        load(fullfile(folderCSD,combineElecString));
                    otherwise
                        load(fullfile(folderLFP,combineElecString));
                    end
                    
                    % Get the relevant trials here
                    if intersectTrials
                        goodPosElec = setdiff(goodPos,elecBadTrials);
                    else
                        goodPosElec = setdiff(goodPos,elecBadTrials1{c});
                    end
                    disp(['No. of good trials for elec' num2str(combineElec(c)) ': ' num2str(length(goodPosElec))]);
                    analogData = analogData(goodPosElec,:);
                    if notchData
                        analogDataNotched = analogDataNotched(goodPosElec,:);
                        analogData = analogDataNotched;
                    end
                    totalNumStim = totalNumStim + length(goodPosElec);
                    
                    if plotStyle==3 % line spectrum
                        [~,SElecBL(:,c),~,~,f2] = getSpectrum(analogData(:,BLPos),mtmParams,movingWin,takeLogTrial,BLMin,BLMax,timeVals,plotStyle);
                        [~,SElecST(:,c),~,~,f2] = getSpectrum(analogData(:,STPos),mtmParams,movingWin,takeLogTrial,BLMin,BLMax,timeVals,plotStyle);
                    else
                        if spectrumType==1 % raw spectrum
                            [~,SElec(:,:,c),~,t2,f2] = getSpectrum(analogData,mtmParams,movingWin,takeLogTrial,BLMin,BLMax,timeVals,plotStyle);
                        else
                            [~,~,dSElec(:,:,c),t2,f2] = getSpectrum(analogData,mtmParams,movingWin,takeLogTrial,BLMin,BLMax,timeVals,plotStyle);
                        end
                    end
                    
                end
                
                
                if plotStyle==3 % line spectrum
                    meanSpecBL = mean(SElecBL,2);
                    meanSpecST = mean(SElecST,2);
                    meanDiffSpec = meanSpecST- meanSpecBL;
                else
                    if spectrumType==1 % raw spectrum
                        meanSpec = mean(SElec,3);
                    else
                        meanDiffSpec = mean(dSElec,3);
                    end
                end

                disp(['No. of trials for combined electrodes: ' num2str(totalNumStim)]);
            end

            if refDiffType ~= 1
                numDiffElec = length(diffElec);
                analogDataD = [];
                totalNumStimD = 0;
                for c = 1:numDiffElec
                    clear analogData analogDataNotched
                    diffElecString = ['elec' num2str(diffElec(c)) '.mat'];
                    switch refType
                    case 1
                        load(fullfile(folderLFP,diffElecString));                    
                    case 2
                        diffElecStr = num2str(diffElec(c));
                        if length(diffElecStr) == 2
                            diffElecString1 = ['elec' diffElecStr(1) '.mat'];
                            diffElecString2 = ['elec' diffElecStr(2) '.mat'];
                        elseif length(diffElecStr) == 3
                            diffElecString1 = ['elec' diffElecStr(1) '.mat'];
                            diffElecString2 = ['elec' diffElecStr(2:3) '.mat'];
                        elseif length(diffElecStr) == 4
                            diffElecString1 = ['elec' diffElecStr(1:2) '.mat'];
                            diffElecString2 = ['elec' diffElecStr(3:4) '.mat'];
                        end
                        clear analogData analogDataNotched
                        load(fullfile(folderLFP,diffElecString1));
                        aD1 = analogData;
                        if notchData
                            aDN1 = analogDataNotched;
                        end
                        clear analogData analogDataNotched
                        load(fullfile(folderLFP,diffElecString2));
                        analogData = aD1 - analogData;
                        if notchData
                            analogDataNotched = aDN1 - analogDataNotched;
                        end
                    case 3
                        load(fullfile(folderAverage,diffElecString));
                        aD1 = analogData;
                        if notchData
                            aDN1 = analogDataNotched;
                        end
                        clear analogData analogDataNotched
                        load(fullfile(folderAverage,'avgRef.mat'));
                        analogData = aD1 - analogData;
                        if notchData
                            analogDataNotched = aDN1 - analogDataNotched;
                        end
                    case 4
                        load(fullfile(folderCSD,diffElecString));
                    otherwise
                        load(fullfile(folderLFP,diffElecString));
                    end
                    
                    % Get the relevant trials here
                    if intersectTrials
                        goodPosElecD = setdiff(goodPos,elecBadTrialsD);
                    else
                        goodPosElecD = setdiff(goodPos,elecBadTrials1D{c});
                    end
                    disp(['No. of good trials for ' num2str(diffElec(c)) ': ' num2str(length(goodPosElecD))]);
                    analogData = analogData(goodPosElecD,:);
                    
                    if notchData
                        analogDataNotched = analogDataNotched(goodPosElecD,:);
                        analogDataD = analogDataNotched;
                    else
                        analogDataD = analogData;
                    end
                    
                    totalNumStimD = totalNumStimD + length(goodPosElecD);
                    
                    if plotStyle==3 % line spectrum
                        [~,SElecBLD(:,c),~,~,f2] = getSpectrum(analogDataD(:,BLPos),mtmParams,movingWin,takeLogTrial,BLMin,BLMax,timeVals,plotStyle);
                        [~,SElecSTD(:,c),~,~,f2] = getSpectrum(analogDataD(:,STPos),mtmParams,movingWin,takeLogTrial,BLMin,BLMax,timeVals,plotStyle);
                    else
                        if spectrumType==1 % raw spectrum
                            [~,SElecD(:,:,c),~,t2,f2] = getSpectrum(analogDataD,mtmParams,movingWin,takeLogTrial,BLMin,BLMax,timeVals,plotStyle);
                        else
                            [~,~,dSElecD(:,:,c),t2,f2] = getSpectrum(analogDataD,mtmParams,movingWin,takeLogTrial,BLMin,BLMax,timeVals,plotStyle);
                        end
                    end
                end
                
                if plotStyle==3 % line spectrum
                    meanSpecBLD = mean(SElecBLD,2);
                    meanSpecSTD = mean(SElecSTD,2);
                    meanDiffSpecD = meanSpecSTD- meanSpecBLD;
                else
                    if spectrumType==1 % raw spectrum
                        meanSpecD = mean(SElecD,3);
                    else
                        meanDiffSpecD = mean(dSElecD,3);
                    end
                end
                
                disp(['No. of trials for combined electrodes (difference): ' num2str(totalNumStimD)]);
            end
                
            
            
            if (refCombineType == 1 && refDiffType == 1) % no combination or difference
            
                %--------------------------------------------------------------
                % Implement the method set by tfMethod variable

                if (tfMethod == 1) % MTM
                    %----------------------------------------------------------
                    if (plotStyle == 3) % line plot (power vs f)
                        if (spectrumType == 1) % raw
                            [SBL,f]=mtspectrumc((analogData(:,BLPos))',mtmParams); % baseline period
                            [SST,f]=mtspectrumc((analogData(:,STPos))',mtmParams); % stimulus period
                            
                            if refDiffType ~= 1
                                [SBLd,f]=mtspectrumc((analogDataD(:,BLPos))',mtmParams); % baseline period
                                [SSTd,f]=mtspectrumc((analogDataD(:,STPos))',mtmParams); % stimulus period
                                
                                SBL = abs(SBL./SBLd);
                                SST = abs(SST./SSTd);
                            end
                            
                            plot(tfplotHandles(k,j),f,conv2Log(SST),'color',plotColor,'Linewidth',plotLineWidth);
                            set(tfplotHandles(k,j),'Nextplot','add');
                            plot(tfplotHandles(k,j),f,conv2Log(SBL),'color','g','Linewidth',plotLineWidth);
                            if holdOnState
                                set(tfplotHandles(k,j),'Nextplot','add');
                            else
                                set(tfplotHandles(k,j),'Nextplot','replace');
                            end

                        elseif (spectrumType == 2) % difference from baseline
                            if xsST == xsBL 
                                [SBL,f]=mtspectrumc((analogData(:,BLPos))',mtmParams);
                                [SST,f]=mtspectrumc((analogData(:,STPos))',mtmParams);
                                
                                if refDiffType ~= 1
                                    [SBLd,f]=mtspectrumc((analogDataD(:,BLPos))',mtmParams); % baseline period
                                    [SSTd,f]=mtspectrumc((analogDataD(:,STPos))',mtmParams); % stimulus period

                                    SBL = abs(SBL./SBLd);
                                    SST = abs(SST./SSTd);
                                end
                                
                                plot(tfplotHandles(k,j),f,10*(conv2Log(SST)-conv2Log(SBL)),'color',plotColor,'Linewidth',plotLineWidth);
                            else
                                disp('Choose same baseline and stimulus periods..');
                            end
                            if holdOnState
                                set(tfplotHandles(k,j),'Nextplot','add');
                            else
                                set(tfplotHandles(k,j),'Nextplot','replace');
                            end
                        end

                    %----------------------------------------------------------    
                    elseif (plotStyle == 1) % pcolor plot

                        if (spectrumType == 1) % raw
                            [S1,t,f]=mtspecgramc((analogData)',movingWin,mtmParams);
                            
                            if refDiffType ~= 1
                                [S1d,t,f]=mtspecgramc((analogDataD)',movingWin,mtmParams);
                                S1 = abs(S1./S1d);
                            end
                            
                            t = t + timeVals(1); % shift the t values to the actual time
                            pcolor(tfplotHandles(k,j),t,f,conv2Log(S1')); shading(tfplotHandles(k,j),'interp');
                            caxis([cmin cmax]);

                        elseif (spectrumType == 2) % difference from baseline
                            [S1,t,f]=mtspecgramc((analogData)',movingWin,mtmParams);
                            
                            if refDiffType ~= 1
                                [S1d,t,f]=mtspecgramc((analogDataD)',movingWin,mtmParams);
                                S1 = abs(S1./S1d);
                            end
                            
                            t = t + timeVals(1); % shift the t values to the actual time

                            % baseline period calculation
                            tBL = (t>=BLMin) & (t<=BLMax); % baseline time indices
                            S1BL = S1(tBL,:); % part of spectrum corresponding to the baseline period
                            mlogS1BL = mean(conv2Log(S1BL),1); % mean log power 
                            % across these time points at every frequency

                            % difference spectrum calculation
                            dS1 = 10*(conv2Log(S1) - repmat(mlogS1BL,size(S1,1),1)); % in dB 
                            % subtract the mean baseline power at each
                            % frequency from the raw power at that frequency at
                            % every time point

                            % plot the difference spectrum
                            pcolor(tfplotHandles(k,j),t,f,dS1'); shading(tfplotHandles(k,j),'interp');
                            caxis([cmin cmax]);
                        end

                    %----------------------------------------------------------    
                    elseif (plotStyle == 2) % imagesc plot

                        if (spectrumType == 1) % raw
                            [S1,t,f]=mtspecgramc((analogData)',movingWin,mtmParams);
                            
                            if refDiffType ~= 1
                                [S1d,t,f]=mtspecgramc((analogDataD)',movingWin,mtmParams);
                                S1 = abs(S1./S1d);
                            end
                            
                            t = t + timeVals(1); % shift the t values to the actual time
                            imagesc(t,f,conv2Log(S1'),'Parent',tfplotHandles(k,j)); 
                            axis(tfplotHandles(k,j),'xy'); % for imagesc every 
                            % parameter has to be set using the plot handles 
                            % explicitly 
                            set(tfplotHandles(k,j),'cLim',[cmin cmax]);

                        elseif (spectrumType == 2) % difference from baseline
                            [S1,t,f]=mtspecgramc((analogData)',movingWin,mtmParams);
                            
                            if refDiffType ~= 1
                                [S1d,t,f]=mtspecgramc((analogDataD)',movingWin,mtmParams);
                                S1 = abs(S1./S1d);
                            end
                            
                            t = t + timeVals(1); % shift the t values to the actual time

                            % baseline period calculation
                            tBL = (t>=BLMin) & (t<=BLMax); % baseline time indices
                            S1BL = S1(tBL,:); % part of spectrum corresponding to the baseline period
                            mlogS1BL = mean(conv2Log(S1BL),1); % mean log power 
                            % across these time points at every frequency

                            % difference spectrum calculation
                            dS1 = 10*(conv2Log(S1) - repmat(mlogS1BL,size(S1,1),1)); % in dB 
                            % subtract the mean baseline power at each
                            % frequency from the raw power at that frequency at
                            % every time point

                            % plot the difference spectrum                       
                            imagesc(t,f,dS1','Parent',tfplotHandles(k,j)); 
                            axis(tfplotHandles(k,j),'xy'); % for imagesc every 
                            % parameter has to be set using the plot handles 
                            % explicitly 
                            set(tfplotHandles(k,j),'cLim',[cmin cmax]);
                        end
                    end
                    
                    %======================================================

                elseif (tfMethod == 2) % MP Algorithm
                    
                    disp(['number of atoms for reconstruction:' num2str(numAtomsMP)]);
                    
                    % Vinay - create a directory (if not already present)
                    % to store the saved mp energy data
                    refString = {'bipolar','average','csd'};
                    if isunix
                        if refType == 1
                            folderSaveMP = [folderName 'mpSavedSpectrum/'];
                        else
                            folderSaveMP = [folderName 'mpSavedSpectrum/' refString{refType-1} '/'];
                        end
                        makeDirectory(folderSaveMP);                        
                    else
                        if refType == 1
                            folderSaveMP = [folderName 'mpSavedSpectrum\'];
                        else
                            folderSaveMP = [folderName 'mpSavedSpectrum\' refString{refType-1} '\'];
                        end
                        makeDirectory(folderSaveMP);
                    end
                    
%                     tagPos = (['Ca' num2str(a(3)) 'e' num2str(e(3)) 's' num2str(s(3)) 'f' num2str(f(3)) 'o' num2str(o(3)) 'c' num2str(c(3)) 't' num2str(t(3)) 'p' num2str(p(3)) 'r' num2str(r(3)) ...
%                                    '_Ra' num2str(a(2)) 'e' num2str(e(2)) 's' num2str(s(2)) 'f' num2str(f(2)) 'o' num2str(o(2)) 'c' num2str(c(2)) 't' num2str(t(2)) 'p' num2str(p(2)) 'r' num2str(r(2)) ...
%                                    '_Sa' num2str(a(1)) 'e' num2str(e(1)) 's' num2str(s(1)) 'f' num2str(f(1)) 'o' num2str(o(1)) 'c' num2str(c(1)) 't' num2str(t(1)) 'p' num2str(p(1)) 'r' num2str(r(1))]);
                    
                    if strncmp(protocolName,'GRF',3)
                        tagPos = ['L' tagPos1];
                    elseif strncmp(protocolName,'CRS',3)
                        tagPos = ['C' tagPos3 'R' tagPos2 'S' tagPos1];
                    end

                                                                          
                    if useDiff
                        
                        if isunix
                            folderSaveMP = [folderSaveMP 'diff/'];
                        else
                            folderSaveMP = [folderSaveMP 'diff\'];
                        end
                        makeDirectory(folderSaveMP);
                        
                        
                        if strncmp(protocolName,'GRF',3)
                            mpSavedSpectrumName = (['mpSpectrum' num2str(numAtomsMP) 'atoms' '_notch' num2str(notchData)...
                                'pos' tagPos '_protoGRF'...
                                '_refType' num2str(refType)...
                                '_diff' num2str(channelString) num2str(analogChannelString2)...
                                '_vary' num2str(paramNum1) num2str(paramNum2)]);
                            
                        elseif strncmp(protocolName,'CRS',3)    
                            mpSavedSpectrumName = (['mpSpectrum' num2str(numAtomsMP) 'atoms' '_notch' num2str(notchData)...
                                'pos' tagPos '_proto' num2str(loadProtocolNumber)...
                                '_refType' num2str(refType)...
                                '_diff' num2str(channelString) num2str(analogChannelString2)...
                                '_vary' num2str(gaborNum1) num2str(paramNum1) num2str(gaborNum2) num2str(paramNum2)]);
                        end
                                                    
                        % Vinay - save or load saved MP data only if the
                        % saveMPFlag is set, else just do all the
                        % computations
                        if (exist([folderSaveMP mpSavedSpectrumName '.mat'],'file') && saveMPFlag)
                            disp('loading saved spectrum data...');
                            load([folderSaveMP mpSavedSpectrumName '.mat']);
                            L = size(analogData,2);
                        else     
                            clear tag gaborInfo X
                            tag = 'temp/';

                            disp('Running MP for the bipolar case...');

                            % Import the data
                            X = analogData';
                            L = size(analogData,2);
                            signalRange = [1 L]; % full range
                            importData(X,mpFolderTemp,tag,signalRange,Fs);

                            % perform Gabor decomposition
                            Numb_points = L; % length of the signal
                            Max_iterations = 100; % number of iterations

                            disp(['MP decomposition with max iterations =' num2str(Max_iterations)]);

                            runGabor(mpFolderTemp,tag,Numb_points, Max_iterations);
                            disp('MP decomposition done');

                            gaborInfo = getGaborData(mpFolderTemp,tag,1);

                            gaborInfoGoodPos = gaborInfo; % changed from gaborInfo(goodPos) because all are goodPos trials already

                            wrap = [];
                            atomList = (1:numAtomsMP);
                            mpSpectrum = [];
                            disp(['Reconstructing Energy from:' num2str(numAtomsMP) 'atoms, and'  num2str(length(goodPos)) 'trials']);
                            for m=1:length(goodPos) % length(goodPos) is the same as the number of trials, so left unchanged here 
                                disp(['trial number:' num2str(m) '(actual trial -)' num2str(goodPos(m))]);
                                rEnergy = reconstructEnergyFromAtomsMPP(gaborInfoGoodPos{m}.gaborData,L,wrap,atomList);
                                if m == 1
                                    mpSpectrum = rEnergy;
                                else
                                    mpSpectrum = mpSpectrum + rEnergy;
                                end
                            end
                            mpSpectrum = mpSpectrum/length(goodPos);
                        
                            if saveMPFlag
                                disp('saving spectrum data...');
                                save([folderSaveMP mpSavedSpectrumName], 'mpSpectrum');
                            end
                        end
                        
                    else
                        
                        if isunix
                            folderSaveMP = [folderSaveMP num2str(channelString) '/'];
                        else
                            folderSaveMP = [folderSaveMP num2str(channelString) '\'];
                        end
                        makeDirectory(folderSaveMP);
                        
                        
                        if strncmp(protocolName,'GRF',3)
                            mpSavedSpectrumName = (['mpSpectrum' num2str(numAtomsMP) 'atoms' '_notch' num2str(notchData)...
                                'pos_' tagPos '_protoGRF' '_refType' num2str(refType)...
                                '_refType' num2str(refType)...
                                '_vary' num2str(paramNum1) num2str(paramNum2)]);
                        elseif strncmp(protocolName,'CRS',3)
                            mpSavedSpectrumName = (['mpSpectrum' num2str(numAtomsMP) 'atoms' '_notch' num2str(notchData)...
                                'pos_' tagPos '_proto' num2str(loadProtocolNumber) '_refType' num2str(refType)...
                                '_refType' num2str(refType)...
                                '_vary' num2str(gaborNum1) num2str(paramNum1) num2str(gaborNum2) num2str(paramNum2)]);
                        end
                        
                        
                        if (exist([folderSaveMP mpSavedSpectrumName '.mat'],'file') && saveMPFlag)
                            disp('loading saved spectrum data...');
                            load([folderSaveMP mpSavedSpectrumName '.mat']);
                            L = size(analogData,2);
                        else
                            clear tag gaborInfo
                            tag = channelString;
                            mpDataExists = 0;
                            if ispc
                                if exist([mpFolder tag '\gaborInfo.mat'],'file')
                                    load([mpFolder tag '\gaborInfo.mat']);
                                    mpDataExists = 1;
                                end
                            else
                                if exist([mpFolder tag '/gaborInfo.mat'],'file')
                                    load([mpFolder tag '/gaborInfo.mat']);
                                    mpDataExists = 1;
                                end
                            end
                            
                            if ~mpDataExists
                                clear tag gaborInfo X
                                tag = 'temp/';

                                disp(['Running MP...Reference Type:' refString{refType-1}]);

                                % Import the data
                                X = analogData';
                                L = size(analogData,2);
                                signalRange = [1 L]; % full range
                                importData(X,mpFolderTemp,tag,signalRange,Fs);

                                % perform Gabor decomposition
                                Numb_points = L; % length of the signal
                                Max_iterations = 100; % number of iterations

                                disp(['MP decomposition with max iterations =' num2str(Max_iterations)]);

                                runGabor(mpFolderTemp,tag,Numb_points, Max_iterations);
                                disp('MP decomposition done');
                                
                                gaborInfo = getGaborData(mpFolderTemp,tag,1);
                            end

                            gaborInfoGoodPos = gaborInfo; % changed from gaborInfo(goodPos) because all are goodPos trials already

                            L = size(analogData,2);
                            wrap = [];
                            atomList = (1:numAtomsMP);
                            mpSpectrum = [];
                            disp(['Reconstructing Energy from:' num2str(numAtomsMP) 'atoms, and'  num2str(length(goodPos)) 'trials']);
                            for m=1:length(goodPos)% length(goodPos) is the same as the number of trials, so left unchanged here
                                disp(['trial number:' num2str(m) '(actual trial -)' num2str(goodPos(m))]);
                                rEnergy = reconstructEnergyFromAtomsMPP(gaborInfoGoodPos{m}.gaborData,L,wrap,atomList);
                                if m == 1
                                    mpSpectrum = rEnergy;
                                else
                                    mpSpectrum = mpSpectrum + rEnergy;
                                end
                            end
                            mpSpectrum = mpSpectrum/length(goodPos);
                            
                            if saveMPFlag
                                disp('saving spectrum data...');
                                save([folderSaveMP mpSavedSpectrumName], 'mpSpectrum');
                            end
                        
                        end
                    end
                    
                    mpSpectrum = mpSpectrum'; % transpose so that first index corresponds to 'time' and second to 'freq'
                    
                    t = timeVals;
                    f = 0:Fs/L:Fs/2; 
                    
                    % plot the MP Spectrum
                    if plotStyle == 3 % line plot
                        disp('No line plot for MP method');
                        
                    elseif plotStyle == 1 % pcolor plot
                        
                        if (spectrumType == 1) % raw
                            pcolor(tfplotHandles(k,j),t,f,conv2Log(mpSpectrum')); 
                            shading(tfplotHandles(k,j),'interp');
                            caxis([cmin cmax]);
                            
                        elseif(spectrumType == 2) % difference from baseline                            
                            % baseline period calculation
                            tBL = (t>=BLMin) & (t<=BLMax); % baseline time indices
                            mpSpectrumBL = mpSpectrum(tBL,:); % part of spectrum corresponding to the baseline period
                            mlogmpSpectrumBL = mean(conv2Log(mpSpectrumBL),1); % mean log power 
                            % across these time points at every frequency

                            % difference spectrum calculation
                            dmpSpectrum = 10*(conv2Log(mpSpectrum) - repmat(mlogmpSpectrumBL,size(mpSpectrum,1),1)); % in dB 
                            % subtract the mean baseline power at each
                            % frequency from the raw power at that frequency at
                            % every time point

                            % plot the difference spectrum
                            pcolor(tfplotHandles(k,j),t,f,dmpSpectrum'); 
                            shading(tfplotHandles(k,j),'interp');
                            caxis([cmin cmax]);
                        end
                        
                    elseif plotStyle == 2 % imagesc plot
                        
                        if (spectrumType == 1) % raw
                            imagesc(t,f,conv2Log(mpSpectrum'),'Parent',tfplotHandles(k,j)); 
                            axis(tfplotHandles(k,j),'xy'); % for imagesc every 
                            % parameter has to be set using the plot handles 
                            % explicitly 
                            set(tfplotHandles(k,j),'cLim',[cmin cmax]);
                            
                        elseif(spectrumType == 2) % difference from baseline                            
                            % baseline period calculation
                            tBL = (t>=BLMin) & (t<=BLMax); % baseline time indices
                            mpSpectrumBL = mpSpectrum(tBL,:); % part of spectrum corresponding to the baseline period
                            mlogmpSpectrumBL = mean(conv2Log(mpSpectrumBL),1); % mean log power 
                            % across these time points at every frequency

                            % difference spectrum calculation
                            dmpSpectrum = 10*(conv2Log(mpSpectrum) - repmat(mlogmpSpectrumBL,size(mpSpectrum,1),1)); 
                            % subtract the mean baseline power at each
                            % frequency from the raw power at that frequency at
                            % every time point

                            % plot the difference spectrum
                            imagesc(t,f,dmpSpectrum','Parent',tfplotHandles(k,j)); 
                            axis(tfplotHandles(k,j),'xy'); % for imagesc every 
                            % parameter has to be set using the plot handles 
                            % explicitly 
                            set(tfplotHandles(k,j),'cLim',[cmin cmax]);
                        end
                        
                    end
                    
                
%------------------------------------------------------------------------                  
                    elseif (tfMethod == 3) % HHT
                    
                    disp(['HHT calculation...Noise ratio: ' num2str(Nstd) ', Ensemble number: ' num2str(NE)]);
                    
                    % Vinay - create a directory (if not already present)
                    % to store the saved hht energy data
                    refString = {'bipolar','average','csd'};
                    if isunix
                        if refType == 1
                            folderSaveHHT = [folderName 'hhtSavedSpectrum/'];
                        else
                            folderSaveHHT = [folderName 'hhtSavedSpectrum/' refString{refType-1} '/'];
                        end
                        makeDirectory(folderSaveHHT);                        
                    else
                        if refType == 1
                            folderSaveHHT = [folderName 'hhtSavedSpectrum\'];
                        else
                            folderSaveHHT = [folderName 'hhtSavedSpectrum\' refString{refType-1} '\'];
                        end
                        makeDirectory(folderSaveHHT);
                    end
                    
                    % Required parameters for implementation
                    t0=timeVals(1);         % starting time value of the signal
                    t1=timeVals(end);       % last time value of the signal
                    fw0=0;                  % minimum frequency
                    fw1=250;                % maximum frequency
                    fres=251;               % frequency resolution i.e. no. of frequencies between the maximum and minimum frequencies
                    tres=length(timeVals);  % time resolution i.e. no. of time points between the maximum and minimum time values
                    tw0=t0;                 % minimum time for zooming
                    tw1=t1;                 % maximum time for zooming
                    lscale=[];              % => linear axis
                    nfilter=0;              % number of points to be median filtered
                    numIMFs = 9;            % number of IMFs
                    ifmethod='hilbert';     % instantaneous frequency method - default
                    normmethod='hilbert';   % normalization method - default
                    fastHHT = 1;            % Use FEEMD if set, EEMD if 0
                    
%                     tagPos = (['Ca' num2str(a(3)) 'e' num2str(e(3)) 's' num2str(s(3)) 'f' num2str(f(3)) 'o' num2str(o(3)) 'c' num2str(c(3)) 't' num2str(t(3)) 'p' num2str(p(3)) 'r' num2str(r(3)) ...
%                                    '_Ra' num2str(a(2)) 'e' num2str(e(2)) 's' num2str(s(2)) 'f' num2str(f(2)) 'o' num2str(o(2)) 'c' num2str(c(2)) 't' num2str(t(2)) 'p' num2str(p(2)) 'r' num2str(r(2)) ...
%                                    '_Sa' num2str(a(1)) 'e' num2str(e(1)) 's' num2str(s(1)) 'f' num2str(f(1)) 'o' num2str(o(1)) 'c' num2str(c(1)) 't' num2str(t(1)) 'p' num2str(p(1)) 'r' num2str(r(1))]);
                    
                    if strncmp(protocolName,'GRF',3)
                        tagPos = ['L' tagPos1];
                    elseif strncmp(protocolName,'CRS',3)
                        tagPos = ['C' tagPos3 'R' tagPos2 'S' tagPos1];
                    end
                    
                    clear E
                                                                          
                    if useBipolar
                        
                        if isunix
                            folderSaveHHT = [folderSaveHHT 'bipolar/'];
                        else
                            folderSaveHHT = [folderSaveHHT 'bipolar\'];
                        end
                        makeDirectory(folderSaveHHT);
                        
                        
                        if strncmp(protocolName,'GRF',3)
                            hhtSavedSpectrumName = (['hhtSpectrum_Nstd' num2str(Nstd) '_NE' num2str(NE) '_numIMFs' num2str(numIMFs) '_notch' num2str(notchData)...
                                'pos' tagPos '_protoGRF'...
                                '_refType' num2str(refType)...
                                '_diff' num2str(channelString) num2str(analogChannelString2)...
                                '_vary' num2str(paramNum1) num2str(paramNum2)]);
                        elseif strncmp(protocolName,'CRS',3)
                            hhtSavedSpectrumName = (['hhtSpectrum_Nstd' num2str(Nstd) '_NE' num2str(NE) '_numIMFs' num2str(numIMFs) '_notch' num2str(notchData)...
                                'pos' tagPos '_proto' num2str(protocolNumber)...
                                '_refType' num2str(refType)...
                                '_diff' num2str(channelString) num2str(analogChannelString2)...
                                '_vary' num2str(gaborNum1) num2str(paramNum1) num2str(gaborNum2) num2str(paramNum2)]);
                        end
                        
                        
                        % Vinay - save or load saved HHT data only if the
                        % saveHHTFlag is set, else just do all the
                        % computations
                        if (exist([folderSaveHHT hhtSavedSpectrumName '.mat'],'file') && saveHHTFlag)
                            disp('loading saved spectrum data...');
                            load([folderSaveHHT hhtSavedSpectrumName '.mat']);
                        
                        else
                            disp(['Computing HHT spectrum...no. of trials: ' num2str(numTrials)]);
                            % Initialize the energy distribution matrices
                            E = zeros(fres,tres);
                            numTrials = length(goodPos);
                            % Since goodPos trials have already been
                            % extracted in this case, indexing 
                            % analogData(goodPos(hi),:) changed to 
                            % analogData(hi,:) below
                            for hi=1:numTrials
                                disp(['trial number:' num2str(hi) '(actual trial -)' num2str(goodPos(hi))]);
                                if ~fastHHT
                                    imfList = eemd(analogData(hi,:),Nstd,NE);
                                else
                                    imfList = rcada_eemd(analogData(hi,:),Nstd,NE,numIMFs);
                                end
                                [P,t,f] = nnspe(imfList(:,2:(numIMFs+1)),t0,t1,fres,tres,fw0,fw1,tw0,tw1,ifmethod,normmethod,nfilter,lscale);
                                E=E+P;
                            end
                            
                            E = E/numTrials; % averaged across trials

                            if saveHHTFlag
                                save([folderSaveHHT hhtSavedSpectrumName '.mat'],'E','t','f'); % save the HHT spectrum matrices
                            end
                            
                        end
                        
                    else
                        
                        if isunix
                            folderSaveHHT = [folderSaveHHT num2str(channelString) '/'];
                        else
                            folderSaveHHT = [folderSaveHHT num2str(channelString) '\'];
                        end
                        makeDirectory(folderSaveHHT);
                        
                        
                        if strncmp(protocolName,'GRF',3)
                            hhtSavedSpectrumName = (['hhtSpectrum_Nstd' num2str(Nstd) '_NE' num2str(NE) '_numIMFs' num2str(numIMFs) '_notch' num2str(notchData)...
                                'pos_' tagPos '_protoGRF'...
                                '_refType' num2str(refType)...
                                '_vary' num2str(paramNum1) num2str(paramNum2)]);
                        elseif strncmp(protocolName,'CRS',3)
                            hhtSavedSpectrumName = (['hhtSpectrum_Nstd' num2str(Nstd) '_NE' num2str(NE) '_numIMFs' num2str(numIMFs) '_notch' num2str(notchData)...
                                'pos_' tagPos '_proto' num2str(protocolNumber)...
                                '_refType' num2str(refType)...
                                '_vary' num2str(gaborNum1) num2str(paramNum1) num2str(gaborNum2) num2str(paramNum2)]);
                        end
                        
                        
                        if (exist([folderSaveHHT hhtSavedSpectrumName '.mat'],'file') && saveHHTFlag)
                            disp('loading saved spectrum data...');
                            load([folderSaveHHT hhtSavedSpectrumName '.mat']);
                        else
                            
                            % Initialize the energy distribution matrices
                            E = zeros(fres,tres);
                            numTrials = length(goodPos);
                            disp(['Computing HHT spectrum...no. of trials: ' num2str(numTrials)]);
                            numIMFs = 7;
                            % Since goodPos trials have already been
                            % extracted in this case, indexing 
                            % analogData(goodPos(hi),:) changed to 
                            % analogData(hi,:) below
                            for hi=1:numTrials
                                disp(['trial number:' num2str(hi) '(actual trial -)' num2str(goodPos(hi))]);
                                if ~fastHHT
                                    imfList = eemd(analogData(hi,:),Nstd,NE);
                                    [P,t,f] = nnspe(imfList(:,2:(numIMFs+1)),t0,t1,fres,tres,fw0,fw1,tw0,tw1,ifmethod,normmethod,nfilter,lscale);
                                else
                                    imfList = rcada_eemd(analogData(hi,:),Nstd,NE,numIMFs);
                                    imfList = imfList';
                                    [P,t,f] = nnspe(imfList(:,1:numIMFs),t0,t1,fres,tres,fw0,fw1,tw0,tw1,ifmethod,normmethod,nfilter,lscale);
                                end
                                
                                E=E+P;
                            end
                            
                            E = E/numTrials; % averaged across trials

                            if saveHHTFlag
                                save([folderSaveHHT hhtSavedSpectrumName '.mat'],'E','t','f'); % save the HHT spectrum matrices
                            end
                        
                        end
                    end
                    
                    
                    if gaussFtr~=0
                        %%%%%%%%%%%%%%%%%____Smoothing of the HHT spectra____%%%%%%%%%%%%%%%%%%%%%%

                        gaussFilter = fspecial('gaussian', gaussFtr, 1); % a Gaussian low-pass filter

                        % filter the energy distribution E
                        E = filter2(gaussFilter, E);
                    end
                    
                    E = E'; % transpose for plotting

                    
                    % plot the HHT Spectrum
                    if plotStyle == 3 % line plot

%                         analogData = mean(analogData(goodPos,:),1);
                        % Since goodPos trials have already been extracted
                        analogData = mean(analogData,1);
                        
                        disp('Computing IMFs and displaying them...');
                        if ~fastHHT
                            imfList = eemd(analogData,Nstd,NE);
                        else
                            imfList = rcada_eemd(analogData,Nstd,NE,numIMFs);
                        end
                        imfList = imfList';
                        
%                         figure(3);
                        nimf = 4;
                        set(tfplotHandles(k,j),'Nextplot','add');
                        IMFcolors = {'r','c','b','k'};
                        for m=1:nimf
%                           subplot(nimf,1,m); 
                          plot(tfplotHandles(k,j),t,imfList(m,:),'color',IMFcolors{m},'LineWidth',1.5);
                        end
%                         cFigC = get(cFig,'all');
%                         tfC = get(tfplotHandles(k,j),'Children');
%                         copyobj(cFig,tfplotHandles(k,j));
%                         subplot(nimf,1,1); title('IMFs');

                            if holdOnState
                                set(tfplotHandles(k,j),'Nextplot','add');
                            else
                                set(tfplotHandles(k,j),'Nextplot','replace');
                            end
                    
                        %--------------------------------------------------
                    
                    elseif plotStyle == 1 % pcolor plot
                        
                        if (spectrumType == 1) % raw
                            pcolor(tfplotHandles(k,j),t,f,conv2Log(E)); 
                            shading(tfplotHandles(k,j),'interp');
                            caxis([cmin cmax]);
                            
                        elseif(spectrumType == 2) % difference from baseline                            
                            % baseline period calculation
                            tBL = (t>=BLMin) & (t<=BLMax); % baseline time indices
                            EBL = E(tBL,:); % part of spectrum corresponding to the baseline period
                            mlogEBL = mean(conv2Log(E),1); % mean log power 
                            % across these time points at every frequency

                            % difference spectrum calculation
                            dhhtSpectrum = 10*(conv2Log(E) - repmat(mlogEBL,size(E,1),1)); % in dB 
                            % subtract the mean baseline power at each
                            % frequency from the raw power at that frequency at
                            % every time point

                            % plot the difference spectrum
                            pcolor(tfplotHandles(k,j),t,f,dhhtSpectrum); 
                            shading(tfplotHandles(k,j),'interp');
                            caxis([cmin cmax]);
                        end
                        
                        %--------------------------------------------------
                        
                    elseif plotStyle == 2 % imagesc plot
                        
                        if (spectrumType == 1) % raw
                            imagesc(t,f,conv2Log(E'),'Parent',tfplotHandles(k,j)); 
                            axis(tfplotHandles(k,j),'xy'); % for imagesc every 
                            % parameter has to be set using the plot handles 
                            % explicitly 
                            set(tfplotHandles(k,j),'cLim',[cmin cmax]);
                            
                        elseif(spectrumType == 2) % difference from baseline                            
                            % baseline period calculation
                            tBL = (t>=BLMin) & (t<=BLMax); % baseline time indices
                            EBL = E(tBL,:); % part of spectrum corresponding to the baseline period
                            mlogEBL = mean(conv2Log(EBL),1); % mean log power 
                            % across these time points at every frequency

                            % difference spectrum calculation
                            dhhtSpectrum = 10*(conv2Log(E) - repmat(mlogEBL,size(E,1),1)); 
                            % subtract the mean baseline power at each
                            % frequency from the raw power at that frequency at
                            % every time point

                            % plot the difference spectrum
                            imagesc(t,f,dhhtSpectrum','Parent',tfplotHandles(k,j)); 
                            axis(tfplotHandles(k,j),'xy'); % for imagesc every 
                            % parameter has to be set using the plot handles 
                            % explicitly 
                            set(tfplotHandles(k,j),'cLim',[cmin cmax]);
                        end
                        
                    end
                   
                end
                
                
            else
                
                if refDiffType ~= 1
                    
                    if spectrumType==1 % raw
                        if plotStyle==3 % line spectrum
                            meanSpecST = meanSpecST - meanSpecSTD;
                            meanSpecBL = meanSpecBL - meanSpecBLD;
                        else
                            meanSpec = 10*(meanSpec - meanSpecD); % in dB
                        end
                    else % difference
                        if plotStyle==3 % line spectrum
                            meanSpecST = meanSpecST - meanSpecSTD;
                            meanSpecBL = meanSpecBL - meanSpecBLD;
                        else
                            meanDiffSpec = meanDiffSpec - meanDiffSpecD; % meanDiffSpec and meanDiffSpecD are already in dB
                        end
                    end
                    
                    
                    if plotStyle==3 % line spectrum plot
                        
                        if spectrumType==1 % raw
                            plot(tfplotHandles(k,j),f2,meanSpecST,'color',plotColor,'Linewidth',plotLineWidth);
                            set(tfplotHandles(k,j),'Nextplot','add');
                            plot(tfplotHandles(k,j),f2,meanSpecBL,'color','g','Linewidth',plotLineWidth);
                            if holdOnState
                                set(tfplotHandles(k,j),'Nextplot','add');
                            else
                                set(tfplotHandles(k,j),'Nextplot','replace');
                            end
                        else
                            plot(tfplotHandles(k,j),f2,10*(meanSpecST-meanSpecBL),'color',plotColor,'Linewidth',plotLineWidth); % in dB
                            if holdOnState
                                set(tfplotHandles(k,j),'Nextplot','add');
                            else
                                set(tfplotHandles(k,j),'Nextplot','replace');
                            end
                        end
                        
                    elseif plotStyle==1 % pcolor
                        
                        if spectrumType==1 % raw
                            colormap('default');
                            pcolor(tfplotHandles(k,j),t2,f2,meanSpec'); 
                            shading(tfplotHandles(k,j),'interp');
                            caxis([cmin cmax]);
                        else
                            colormap('default');
                            pcolor(tfplotHandles(k,j),t2,f2,meanDiffSpec'); 
                            shading(tfplotHandles(k,j),'interp');
                            caxis([cmin cmax]);
                        end
                        
                    elseif plotStyle==2 %imagesc
                        
                        if spectrumType==1 %raw
                            colormap('default');
                            imagesc(t2,f2,meanSpec','Parent',tfplotHandles(k,j)); 
                            axis(tfplotHandles(k,j),'xy'); % for imagesc every 
                            % parameter has to be set using the plot handles 
                            % explicitly 
                            set(tfplotHandles(k,j),'cLim',[cmin cmax]);
                        else
                            colormap('default');
                            imagesc(t2,f2,meanDiffSpec','Parent',tfplotHandles(k,j)); 
                            axis(tfplotHandles(k,j),'xy'); % for imagesc every 
                            % parameter has to be set using the plot handles 
                            % explicitly 
                            set(tfplotHandles(k,j),'cLim',[cmin cmax]);
                        end
                    end
                    
                else
                    
                    if plotStyle==3 % line spectrum plot
                        
                        if spectrumType==1 % raw
                            plot(tfplotHandles(k,j),f2,meanSpecST,'color',plotColor,'Linewidth',plotLineWidth);
                            set(tfplotHandles(k,j),'Nextplot','add');
                            plot(tfplotHandles(k,j),f2,meanSpecBL,'color','g','Linewidth',plotLineWidth);
                            if holdOnState
                                set(tfplotHandles(k,j),'Nextplot','add');
                            else
                                set(tfplotHandles(k,j),'Nextplot','replace');
                            end
                        else
                            plot(tfplotHandles(k,j),f2,10*(meanSpecST-meanSpecBL),'color',plotColor,'Linewidth',plotLineWidth); % in dB
                            if holdOnState
                                set(tfplotHandles(k,j),'Nextplot','add');
                            else
                                set(tfplotHandles(k,j),'Nextplot','replace');
                            end
                        end
                        
                    elseif plotStyle==1 % pcolor
                        
                        if spectrumType==1 % raw
                            colormap('default');
                            pcolor(tfplotHandles(k,j),t2,f2,meanSpec'); 
                            shading(tfplotHandles(k,j),'interp');
                            caxis([cmin cmax]);
                        else
                            colormap('default');
                            pcolor(tfplotHandles(k,j),t2,f2,meanDiffSpec'); 
                            shading(tfplotHandles(k,j),'interp');
                            caxis([cmin cmax]);
                        end
                        
                    elseif plotStyle==2 %imagesc
                        
                        if spectrumType==1 %raw
                            colormap('default');
                            imagesc(t2,f2,meanSpec','Parent',tfplotHandles(k,j)); 
                            axis(tfplotHandles(k,j),'xy'); % for imagesc every 
                            % parameter has to be set using the plot handles 
                            % explicitly 
                            set(tfplotHandles(k,j),'cLim',[cmin cmax]);
                        else
                            colormap('default');
                            imagesc(t2,f2,meanDiffSpec','Parent',tfplotHandles(k,j)); 
                            axis(tfplotHandles(k,j),'xy'); % for imagesc every 
                            % parameter has to be set using the plot handles 
                            % explicitly 
                            set(tfplotHandles(k,j),'cLim',[cmin cmax]);
                        end
                    end
                end
                
            end
  
                
                if strncmp(protocolName,'CRS',3)
                    titleGabor1 = titleGabor(gaborNum1);
                    titleGabor2 = titleGabor(gaborNum2);
                elseif strncmp(protocolName,'GRF',3)
                    titleGabor1 = 'gabor';
                    titleGabor2 = 'gabor';
                end

                % Display title
                if (j==1 && k==1)
                    title(tfplotHandles(k,j),[titleGabor1 '-' titleParam1 'vs' titleGabor2 '-' titleParam2 ':' num2str(titleList1(k)) ',' num2str(titleList2(j))],'FontSize',titleFontSize);
                else
                    title(tfplotHandles(k,j),[num2str(titleList1(k)) ',' num2str(titleList2(j))],'FontSize',titleFontSize);
                end

                % show number of trials for this condition
                if refCombineType==1
                    text(0.1,0.9,['n = ' num2str(length(goodPos))],'unit','normalized','fontsize',7,'Parent',tfplotHandles(k,j));
                elseif refDiffType==1
                    text(0.1,0.9,['n = ' num2str(totalNumStim)],'unit','normalized','fontsize',7,'Parent',tfplotHandles(k,j));
                else
                    text(0.1,0.9,['Combine n = ' num2str(totalNumStim) '; Diff n = ' num2str(totalNumStimD)],'unit','normalized','fontsize',7,'Parent',tfplotHandles(k,j));
                end
                    
                
            end  
            end
        end
        
end


%--------------------------------------------------------------------------
%====================Spectrum calculation==================================
function [S1,mlogS1,dS1,t2,f2] = getSpectrum(data,mtmParams,movingWin,takeLogTrial,BLMin,BLMax,timeVals,specType)

if ~exist ('takeLogTrial','var')
    takeLogTrial = 0;
end

if ~exist ('BLMin','var')
    BLMin = -0.5;
end

if ~exist ('BLMax','var')
    BLMax = 0;
end

if ~exist('mtmParams','var')
    mtmParams.Fs = 2000;
    mtmParams.tapers=[2 3]; % [1 1] case is simply stft with dpss window
    mtmParams.trialave=0;
    mtmParams.err=0;
    mtmParams.pad=-1;
end

if ~exist('movingWin','var')
    movingWin = [0.5 0.01];
end

if takeLogTrial
    
    if specType==3 % line spectrum i.e. normal stft
        
        mtmParams.trialave = 0; % don't average spectrum across trials
        [S1,f2]=mtspectrumc(data',mtmParams);
        
        logS1 = conv2Log(S1); % log power for each trial
        mlogS1 = mean(logS1,2); % mean of log power across trials
        
        dS1 = [];
        t2 = [];
        
    else
        
        mtmParams.trialave = 0; % don't average spectrum across trials
        [S1,t2,f2]=mtspecgramc(data',movingWin,mtmParams);

        logS1 = conv2Log(S1); % log power for each trial
        trialmlogS1 = mean(logS1,3); % mean of log power across trials

        t2 = t2 + timeVals(1); % shift the t values to the actual time
        tBL = (t2>=BLMin) & (t2<=BLMax); % baseline time indices
        trialmlogS1BL = trialmlogS1(tBL,:); % baseline trial mean log power

        mlogS1BL = mean(trialmlogS1BL,1); % mean baseline log power

        % difference spectrum calculation
        dS1 = 10*(trialmlogS1 - repmat(mlogS1BL,size(trialmlogS1,1),1));
        
        
        % mean log power across trials
        mlogS1 = trialmlogS1;
    end
    
else
    
    if specType==3 % line spectrum i.e. normal stft
        
        mtmParams.trialave = 1; % average spectrum across trials
        [S1,f2]=mtspectrumc(data',mtmParams);
        
        mlogS1 = conv2Log(S1); % log power for each trial
        
        dS1 = [];
        t2= [];
        
    else
        
        mtmParams.trialave = 1; % average spectrum across trials
        [S1,t2,f2]=mtspecgramc(data',movingWin,mtmParams);
        t2 = t2 + timeVals(1); % shift the t values to the actual time

        tBL = (t2>=BLMin) & (t2<=BLMax); % baseline time indices
        S1BL = S1(tBL,:); % part of spectrum corresponding to the baseline period
        mlogS1BL = mean(conv2Log(S1BL),1); % mean log power across these time points at every frequency

        % difference spectrum calculation
        dS1 = 10*(conv2Log(S1) - repmat(mlogS1BL,size(S1,1),1)); % in dB
        
        % mean log power across trials
        mlogS1 = conv2Log(S1);
    end
    
end

end

