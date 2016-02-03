% Make a list of all stimulus conditions
% Vinay Shirhatti, 10 Dec 2015
% *************************************************************************

runForIndex = []; % run for all indices if empty

usePipes = 1;

subjectName = 'alpa';
gridType = 'Microelectrode';

[expDates,protocolNames,stimType] = eval(['allProtocols' upper(subjectName(1)) subjectName(2:end) gridType]);

labUbuntu = 1;
if ispc
    folderSourceString = 'K:\'; % Vinay - changed source directory
else
    if labUbuntu
        folderSourceString = '/media/vinay/SRLHD02M/';
    else
        folderSourceString = '/media/store/';
    end
end

if isempty(runForIndex)
    runForIndex = 1:length(protocolNames);
end

% ignoreProtocolsList = [1:10];
% 
% runForIndex = diff(runForIndex,ignoreProtocolsList);
writeFile = 'protocolsConditionsDetailsSheet.csv';
headings = {'index','expDate','protocol','gabor','azimuth','elevation','sigma','spatialFreq','orientation','contrast','temporalFreq','spatialPhase','radius'};
% headings = {'index';'expDate';'protocol';'azimuth';'elevation';'sigma';'spatialFreq';'orientation';'contrast';'temporalFreq';'spatialPhase';'radius'};
% xlswrite(xlFile,headings,1,'A1');
% tb = table(headings);
% writetable(tb,xlFile);
fid = fopen(writeFile,'w');

fprintf(fid,'%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s\n',headings{1,:});

for i=1:length(runForIndex)
    clear azimuthString elevationString sigmaString spatialFreqString orientationString contrastString temporalFreqString spatialPhaseString radiusString
    j = runForIndex(i);
    disp(['working on index: ' num2str(j)]);
    expDate = expDates{j};
    protocolName = protocolNames{j};
    
    folderName = fullfile(folderSourceString,'data',subjectName,gridType,expDate,protocolName);
    folderExtract = fullfile(folderName,'extractedData');
    parameterCombinationsFile = fullfile(folderExtract,'parameterCombinations.mat');
    
    if exist(parameterCombinationsFile,'file')
        [~,aValsUnique,eValsUnique,sValsUnique,fValsUnique,oValsUnique,cValsUnique,tValsUnique,rValsUnique,pValsUnique] = loadParameterCombinations(folderExtract);

        if isempty(rValsUnique)
            rValsUnique = sValsUnique.*3;
        end

        if isempty(pValsUnique)
            pValsUnique = 0;
        end

        if strncmp(protocolName,'GRF',3)
            azimuthString = getStringFromValues(aValsUnique,1);
            elevationString = getStringFromValues(eValsUnique,1);
            sigmaString = getStringFromValues(sValsUnique,1);
            spatialFreqString = getStringFromValues(fValsUnique,1);
            orientationString = getStringFromValues(oValsUnique,1);
            contrastString = getStringFromValues(cValsUnique,1);
            temporalFreqString = getStringFromValues(tValsUnique,1);
            spatialPhaseString = getStringFromValues(pValsUnique,1);
            radiusString = getStringFromValues(rValsUnique,1);
        elseif strncmp(protocolName,'CRS',3)
            for gaborNum=1:3
                azimuthString{gaborNum} = getStringFromValues(aValsUnique,1, gaborNum);
                elevationString{gaborNum} = getStringFromValues(eValsUnique,1, gaborNum);
                sigmaString{gaborNum} = getStringFromValues(sValsUnique,1, gaborNum);
                spatialFreqString{gaborNum} = getStringFromValues(fValsUnique,1, gaborNum);
                orientationString{gaborNum} = getStringFromValues(oValsUnique,1, gaborNum);
                contrastString{gaborNum} = getStringFromValues(cValsUnique,1, gaborNum);
                temporalFreqString{gaborNum} = getStringFromValues(tValsUnique,1, gaborNum);
                spatialPhaseString{gaborNum} = getStringFromValues(pValsUnique,1, gaborNum);
                radiusString{gaborNum} = getStringFromValues(rValsUnique,1, gaborNum);
            end
        end

        if strncmp(protocolName,'GRF',3)    
            entries = {num2str(j),expDate,protocolName,'leftGabor',azimuthString,elevationString,sigmaString,spatialFreqString,orientationString,...
                contrastString,temporalFreqString,spatialPhaseString,radiusString};
    %             tb = table(entries);
    %             writetable(tb,xlFile,'-append');
    %             dlmwrite(xlFile,entries,'delimiter',',','-append');
            fprintf(fid,'%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s\n',entries{1,:});
            
        elseif strncmp(protocolName,'CRS',3)
            for gaborNum=1:3
                switch gaborNum
                    case 1
                        entries = {num2str(j),expDate,protocolName,'Surround',azimuthString{gaborNum},elevationString{gaborNum},sigmaString{gaborNum},...
                            spatialFreqString{gaborNum},orientationString{gaborNum},...
                            contrastString{gaborNum},temporalFreqString{gaborNum},spatialPhaseString{gaborNum},radiusString{gaborNum}};
%                         tb = table(entries);
%                         writetable(tb,xlFile,'-append');
%                         dlmwrite(xlFile,entries,'delimiter',',','-append');
                        fprintf(fid,'%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s\n',entries{1,:});
                    case 2
                        entries = {'','','','Ring',azimuthString{gaborNum},elevationString{gaborNum},sigmaString{gaborNum},...
                            spatialFreqString{gaborNum},orientationString{gaborNum},...
                            contrastString{gaborNum},temporalFreqString{gaborNum},spatialPhaseString{gaborNum},radiusString{gaborNum}};
%                         tb = table(entries);
%                         writetable(tb,xlFile,'-append');
%                         dlmwrite(xlFile,entries,'delimiter',',','-append');
                        fprintf(fid,'%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s\n',entries{1,:});
                    case 3
                        entries = {'','','','Centre',azimuthString{gaborNum},elevationString{gaborNum},sigmaString{gaborNum},...
                            spatialFreqString{gaborNum},orientationString{gaborNum},...
                            contrastString{gaborNum},temporalFreqString{gaborNum},spatialPhaseString{gaborNum},radiusString{gaborNum}};
%                         tb = table(entries);
%                         writetable(tb,xlFile,'-append');
%                         dlmwrite(xlFile,entries,'delimiter',',','-append');
                        fprintf(fid,'%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s\n',entries{1,:});
                end
                
            end
        end
    fprintf(fid,'\n');
    end
    
end

fclose(fid);
