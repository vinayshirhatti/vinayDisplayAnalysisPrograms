% Draw the Stimuli shown in a protocol
%
%
% Vinay Shirhatti, 04 March 2015
%--------------------------------------------------------------------------

function drawProtocolStimuli(subjectName,expDate,protocolName,folderSourceString,gridType,loadProtocolNumber)

% Load defaults
if ~exist('folderSourceString','var')  folderSourceString='/media/store/';        end
if ~exist('gridType','var')            gridType='EEG';                  end
if ~exist('loadProtocolNumber','var')  loadProtocolNumber = 11;         end % Vinay - anything greater than 10 will read the protocolNumber from the saved data

% path to saved data
folderName = fullfile(folderSourceString,'data',subjectName,gridType,expDate,protocolName);
folderExtract = fullfile(folderName,'extractedData');

% Get Combinations
[parameterCombinations,aValsUnique,eValsUnique,sValsUnique,fValsUnique,oValsUnique,cValsUnique,...
    tValsUnique,rValsUnique,pValsUnique] = loadParameterCombinations(folderExtract);


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

rFactor = 3;
fontSizeTiny = 8;

if strncmp(protocolName,'GRF',3)
    
    aLen = length(aValsUnique);
    eLen = length(eValsUnique);
    sLen = length(sValsUnique);
    fLen = length(fValsUnique);
    oLen = length(oValsUnique);
    cLen = length(cValsUnique);
    tLen = length(tValsUnique);
    
%     numStim = aLen*eLen*sLen*fLen*oLen*cLen*tLen;
    
    aB = []; eB = []; sB = []; fB = []; oB = []; cB = []; tB = [];
    index = [];
%     index1 = [];
    count = 1;
    for a = 1:aLen
        for e = 1:eLen
            for s = 1:sLen
                for f = 1:fLen
                    for o = 1:oLen
                        for c = 1:cLen
                            for t = 1:tLen
                                
                                aB = cat(1,aB,aValsUnique(a));
                                eB = cat(1,eB,eValsUnique(e));
                                sB = cat(1,sB,sValsUnique(s));
                                fB = cat(1,fB,fValsUnique(f));
                                oB = cat(1,oB,oValsUnique(o));
                                cB = cat(1,cB,cValsUnique(c));
                                tB = cat(1,tB,tValsUnique(t));
                                
%                                 index1 = cat(1,index1,[a,e,s,f,o,c,t]);
                                
                                index{count} = [a,e,s,f,o,c,t];
                                
                                count = count+1;
                            end
                        end
                    end
                end
            end
        end
    end

    aB = squeeze(reshape(aB,tLen,cLen,oLen,fLen,sLen,eLen,aLen)); % this reshape works fine if the lengths are passed in this reverse order
    eB = squeeze(reshape(eB,tLen,cLen,oLen,fLen,sLen,eLen,aLen));
    sB = squeeze(reshape(sB,tLen,cLen,oLen,fLen,sLen,eLen,aLen));
    fB = squeeze(reshape(fB,tLen,cLen,oLen,fLen,sLen,eLen,aLen));
    oB = squeeze(reshape(oB,tLen,cLen,oLen,fLen,sLen,eLen,aLen));
    cB = squeeze(reshape(cB,tLen,cLen,oLen,fLen,sLen,eLen,aLen));
    tB = squeeze(reshape(tB,tLen,cLen,oLen,fLen,sLen,eLen,aLen));
    
    index = squeeze(reshape(index,tLen,cLen,oLen,fLen,sLen,eLen,aLen));
    
    
    figure;
%     hStim = subplot(1,1,1);
    
    numDim = sum(size(index)~=1);
    if numDim == 1
        for i = 1:size(index,1)
            gaborBackground.azimuthDeg = aB(i);
            gaborBackground.elevationDeg = eB(i);
            gaborBackground.sigmaDeg = sB(i);
            gaborBackground.spatialFreqCPD = fB(i);
            gaborBackground.orientationDeg = oB(i);
            gaborBackground.contrastPC = cB(i);
            gaborBackground.temporalFreqHz = tB(i);
            gaborBackground.spatialPhaseDeg = 0;
            gaborBackground.radiusDeg = 45/rFactor;
            
            gaborRing.azimuthDeg = 0;
            gaborRing.elevationDeg = 0;
            gaborRing.sigmaDeg = 0;
            gaborRing.spatialFreqCPD = 0;
            gaborRing.orientationDeg = 0;
            gaborRing.contrastPC = 0;
            gaborRing.temporalFreqHz = 0;
            gaborRing.spatialPhaseDeg = 0;
            gaborRing.radiusDeg = 0;
            
            hS{i} = subplot(1,size(index,1),i);
            titleString = getTitleGRF(gaborBackground);
            drawStimulus(hS{i},gaborBackground,gaborRing,titleString);
        end
        
    elseif numDim == 2
        for i = 1:size(index,1)
            for j = 1:size(index,2)
                gaborBackground.azimuthDeg = aB(i,j);
                gaborBackground.elevationDeg = eB(i,j);
                gaborBackground.sigmaDeg = sB(i,j);
                gaborBackground.spatialFreqCPD = fB(i,j);
                gaborBackground.orientationDeg = oB(i,j);
                gaborBackground.contrastPC = cB(i,j);
                gaborBackground.temporalFreqHz = tB(i,j);
                gaborBackground.spatialPhaseDeg = 0;
                gaborBackground.radiusDeg = 45/rFactor;
                
                gaborRing.azimuthDeg = 0;
                gaborRing.elevationDeg = 0;
                gaborRing.sigmaDeg = 0;
                gaborRing.spatialFreqCPD = 0;
                gaborRing.orientationDeg = 0;
                gaborRing.contrastPC = 0;
                gaborRing.temporalFreqHz = 0;
                gaborRing.spatialPhaseDeg = 0;
                gaborRing.radiusDeg = 0;
                
                numRows = size(index,1);
                numCols = size(index,2);
                hS{i,j} = subplot(numRows,numCols,(i-1)*numCols+j);
                [~,xString,yString] = getTitleGRF(gaborBackground);
                drawStimulus(hS{i,j},gaborBackground,gaborRing,[]);
                if i==1 % top row
                    title(xString,'Fontsize',fontSizeTiny,'Parent',hS{i,j});
                end
                if j==1 % first column
                    ylabel(yString,'Fontsize',fontSizeTiny,'Rotation',0,'Parent',hS{i,j});
                end
                    
            end
        end
    end
    
    
elseif strncmp(protocolName,'CRS',3)
    
    
    for i = 1:3
        aLen(i) = length(aValsUnique{i});
        eLen(i) = length(eValsUnique{i});
        sLen(i) = length(sValsUnique{i});
        fLen(i) = length(fValsUnique{i});
        oLen(i) = length(oValsUnique{i});
        cLen(i) = length(cValsUnique{i});
        tLen(i) = length(tValsUnique{i});
        pLen(i) = length(pValsUnique{i});
        rLen(i) = length(rValsUnique{i});
    end
    
    if protocolNumber == 10 
        % AFP - radius of Ring is tied to the radius of Centre i.e.
        % although there are many rValsUnique for the Ring Gabor, each
        % value is valid only for certain radii of the centre. Therefore
        % the Rradius has to be checked for validity for every Cradius
        validAnnulusWidth = [0,0.25,0.5,1,2,4]; % set of annulus widths
    end
    
    aB = []; eB = []; sB = []; fB = []; oB = []; cB = []; tB = []; pB = []; rB = [];
    aR = []; eR = []; sR = []; fR = []; oR = []; cR = []; tR = []; pR = []; rR = [];
    index = [];
%     index1 = [];
    count = 1;
    
    switch protocolNumber
        
        case {1, 2, 6, 7, 8}
            for a = 1:aLen(1)
                for e = 1:eLen(1)
                    for s = 1:sLen(1)
                        for f = 1:fLen(1)
                            for o = 1:oLen(1)
                                for c = 1:cLen(1)
                                    for t = 1:tLen(1)
                                        for p = 1:pLen(1)
                                            for r = 1:rLen(1)
                                                
                                                for a2 = 1:aLen(2)
                                                    for e2 = 1:eLen(2)
                                                        for s2 = 1:sLen(2)
                                                            for f2 = 1:fLen(2)
                                                                for o2 = 1:oLen(2)
                                                                    for c2 = 1:cLen(2)
                                                                        for t2 = 1:tLen(2)
                                                                            for p2 = 1:pLen(2)
                                                                                for r2 = 1:rLen(2)
                                                                                    
                                                                                    for r3 = 1:rLen(3)

                                                                                        aB = cat(1,aB,aValsUnique{1}(a));
                                                                                        eB = cat(1,eB,eValsUnique{1}(e));
                                                                                        sB = cat(1,sB,sValsUnique{1}(s));
                                                                                        fB = cat(1,fB,fValsUnique{1}(f));
                                                                                        oB = cat(1,oB,oValsUnique{1}(o));
                                                                                        cB = cat(1,cB,cValsUnique{1}(c));
                                                                                        tB = cat(1,tB,tValsUnique{1}(t));
                                                                                        pB = cat(1,pB,pValsUnique{1}(p));
                                                                                        rB{count} = [0 rValsUnique{1}(r)/rFactor];

                                                                                        aR = cat(1,aR,aValsUnique{2}(a2));
                                                                                        eR = cat(1,eR,eValsUnique{2}(e2));
                                                                                        sR = cat(1,sR,sValsUnique{2}(s2));
                                                                                        fR = cat(1,fR,fValsUnique{2}(f2));
                                                                                        oR = cat(1,oR,oValsUnique{2}(o2));
                                                                                        cR = cat(1,cR,cValsUnique{2}(c2));
                                                                                        tR = cat(1,tR,tValsUnique{2}(t2));
                                                                                        pR = cat(1,pR,pValsUnique{2}(p2));
                                                                                        rR{count} = [rValsUnique{3}(r3)/rFactor rValsUnique{2}(r2)/rFactor];


                                                                                        index{count} = [a,e,s,f,o,c,t,p,r,a2,e2,s2,f2,o2,c2,t2,p2,r2,r3];

                                                                                        count = count+1;
                                                                                            
                                                                                    end
                                                                                        
                                                                                end

                                                                            end
                                                                        end
                                                                    end
                                                                end
                                                            end
                                                        end
                                                    end
                                                end
                                            end

                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end

            aB = squeeze(reshape(aB,rLen(3),rLen(2),pLen(2),tLen(2),cLen(2),oLen(2),fLen(2),sLen(2),eLen(2),aLen(2),rLen(1),pLen(1),tLen(1),cLen(1),oLen(1),fLen(1),sLen(1),eLen(1),aLen(1))); % this reshape works fine if the lengths are passed in this reverse order
            eB = squeeze(reshape(eB,rLen(3),rLen(2),pLen(2),tLen(2),cLen(2),oLen(2),fLen(2),sLen(2),eLen(2),aLen(2),rLen(1),pLen(1),tLen(1),cLen(1),oLen(1),fLen(1),sLen(1),eLen(1),aLen(1)));
            sB = squeeze(reshape(sB,rLen(3),rLen(2),pLen(2),tLen(2),cLen(2),oLen(2),fLen(2),sLen(2),eLen(2),aLen(2),rLen(1),pLen(1),tLen(1),cLen(1),oLen(1),fLen(1),sLen(1),eLen(1),aLen(1)));
            fB = squeeze(reshape(fB,rLen(3),rLen(2),pLen(2),tLen(2),cLen(2),oLen(2),fLen(2),sLen(2),eLen(2),aLen(2),rLen(1),pLen(1),tLen(1),cLen(1),oLen(1),fLen(1),sLen(1),eLen(1),aLen(1)));
            oB = squeeze(reshape(oB,rLen(3),rLen(2),pLen(2),tLen(2),cLen(2),oLen(2),fLen(2),sLen(2),eLen(2),aLen(2),rLen(1),pLen(1),tLen(1),cLen(1),oLen(1),fLen(1),sLen(1),eLen(1),aLen(1)));
            cB = squeeze(reshape(cB,rLen(3),rLen(2),pLen(2),tLen(2),cLen(2),oLen(2),fLen(2),sLen(2),eLen(2),aLen(2),rLen(1),pLen(1),tLen(1),cLen(1),oLen(1),fLen(1),sLen(1),eLen(1),aLen(1)));
            tB = squeeze(reshape(tB,rLen(3),rLen(2),pLen(2),tLen(2),cLen(2),oLen(2),fLen(2),sLen(2),eLen(2),aLen(2),rLen(1),pLen(1),tLen(1),cLen(1),oLen(1),fLen(1),sLen(1),eLen(1),aLen(1)));
            pB = squeeze(reshape(pB,rLen(3),rLen(2),pLen(2),tLen(2),cLen(2),oLen(2),fLen(2),sLen(2),eLen(2),aLen(2),rLen(1),pLen(1),tLen(1),cLen(1),oLen(1),fLen(1),sLen(1),eLen(1),aLen(1)));
            rB = squeeze(reshape(rB,rLen(3),rLen(2),pLen(2),tLen(2),cLen(2),oLen(2),fLen(2),sLen(2),eLen(2),aLen(2),rLen(1),pLen(1),tLen(1),cLen(1),oLen(1),fLen(1),sLen(1),eLen(1),aLen(1)));

            aR = squeeze(reshape(aR,rLen(3),rLen(2),pLen(2),tLen(2),cLen(2),oLen(2),fLen(2),sLen(2),eLen(2),aLen(2),rLen(1),pLen(1),tLen(1),cLen(1),oLen(1),fLen(1),sLen(1),eLen(1),aLen(1))); % this reshape works fine if the lengths are passed in this reverse order
            eR = squeeze(reshape(eR,rLen(3),rLen(2),pLen(2),tLen(2),cLen(2),oLen(2),fLen(2),sLen(2),eLen(2),aLen(2),rLen(1),pLen(1),tLen(1),cLen(1),oLen(1),fLen(1),sLen(1),eLen(1),aLen(1)));
            sR = squeeze(reshape(sR,rLen(3),rLen(2),pLen(2),tLen(2),cLen(2),oLen(2),fLen(2),sLen(2),eLen(2),aLen(2),rLen(1),pLen(1),tLen(1),cLen(1),oLen(1),fLen(1),sLen(1),eLen(1),aLen(1)));
            fR = squeeze(reshape(fR,rLen(3),rLen(2),pLen(2),tLen(2),cLen(2),oLen(2),fLen(2),sLen(2),eLen(2),aLen(2),rLen(1),pLen(1),tLen(1),cLen(1),oLen(1),fLen(1),sLen(1),eLen(1),aLen(1)));
            oR = squeeze(reshape(oR,rLen(3),rLen(2),pLen(2),tLen(2),cLen(2),oLen(2),fLen(2),sLen(2),eLen(2),aLen(2),rLen(1),pLen(1),tLen(1),cLen(1),oLen(1),fLen(1),sLen(1),eLen(1),aLen(1)));
            cR = squeeze(reshape(cR,rLen(3),rLen(2),pLen(2),tLen(2),cLen(2),oLen(2),fLen(2),sLen(2),eLen(2),aLen(2),rLen(1),pLen(1),tLen(1),cLen(1),oLen(1),fLen(1),sLen(1),eLen(1),aLen(1)));
            tR = squeeze(reshape(tR,rLen(3),rLen(2),pLen(2),tLen(2),cLen(2),oLen(2),fLen(2),sLen(2),eLen(2),aLen(2),rLen(1),pLen(1),tLen(1),cLen(1),oLen(1),fLen(1),sLen(1),eLen(1),aLen(1)));
            pR = squeeze(reshape(pR,rLen(3),rLen(2),pLen(2),tLen(2),cLen(2),oLen(2),fLen(2),sLen(2),eLen(2),aLen(2),rLen(1),pLen(1),tLen(1),cLen(1),oLen(1),fLen(1),sLen(1),eLen(1),aLen(1)));
            rR = squeeze(reshape(rR,rLen(3),rLen(2),pLen(2),tLen(2),cLen(2),oLen(2),fLen(2),sLen(2),eLen(2),aLen(2),rLen(1),pLen(1),tLen(1),cLen(1),oLen(1),fLen(1),sLen(1),eLen(1),aLen(1)));

            
            index = squeeze(reshape(index,rLen(3),rLen(2),pLen(2),tLen(2),cLen(2),oLen(2),fLen(2),sLen(2),eLen(2),aLen(2),rLen(1),pLen(1),tLen(1),cLen(1),oLen(1),fLen(1),sLen(1),eLen(1),aLen(1)));
            
            
        case {3, 4, 5, 9}
            
            for a = 1:aLen(2)
                for e = 1:eLen(2)
                    for s = 1:sLen(2)
                        for f = 1:fLen(2)
                            for o = 1:oLen(2)
                                for c = 1:cLen(2)
                                    for t = 1:tLen(2)
                                        for p = 1:pLen(2)
                                            for r = 1:rLen(2)
                                                
                                                for a2 = 1:aLen(3)
                                                    for e2 = 1:eLen(3)
                                                        for s2 = 1:sLen(3)
                                                            for f2 = 1:fLen(3)
                                                                for o2 = 1:oLen(3)
                                                                    for c2 = 1:cLen(3)
                                                                        for t2 = 1:tLen(3)
                                                                            for p2 = 1:pLen(3)
                                                                                for r2 = 1:rLen(3)

                                                                                        aB = cat(1,aB,aValsUnique{2}(a));
                                                                                        eB = cat(1,eB,eValsUnique{2}(e));
                                                                                        sB = cat(1,sB,sValsUnique{2}(s));
                                                                                        fB = cat(1,fB,fValsUnique{2}(f));
                                                                                        oB = cat(1,oB,oValsUnique{2}(o));
                                                                                        cB = cat(1,cB,cValsUnique{2}(c));
                                                                                        tB = cat(1,tB,tValsUnique{2}(t));
                                                                                        pB = cat(1,pB,pValsUnique{2}(p));
                                                                                        rB{count} = [0 rValsUnique{2}(r)/rFactor];
                                                                                        
                                                                                        aR = cat(1,aR,aValsUnique{3}(a2));
                                                                                        eR = cat(1,eR,eValsUnique{3}(e2));
                                                                                        sR = cat(1,sR,sValsUnique{3}(s2));
                                                                                        fR = cat(1,fR,fValsUnique{3}(f2));
                                                                                        oR = cat(1,oR,oValsUnique{3}(o2));
                                                                                        cR = cat(1,cR,cValsUnique{3}(c2));
                                                                                        tR = cat(1,tR,tValsUnique{3}(t2));
                                                                                        pR = cat(1,pR,pValsUnique{3}(p2));
                                                                                        rR{count} = [0 rValsUnique{3}(r2)/rFactor];
                                                                                        

                                                                                        index{count} = [a,e,s,f,o,c,t,p,r,a2,e2,s2,f2,o2,c2,t2,p2,r2];

                                                                                        count = count+1;
                                                                                    
                                                                                end
                                                                            end
                                                                        end
                                                                    end
                                                                end
                                                            end
                                                        end
                                                    end
                                                end
                                                
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end

            aB = squeeze(reshape(aB,rLen(3),pLen(3),tLen(3),cLen(3),oLen(3),fLen(3),sLen(3),eLen(3),aLen(3),rLen(2),pLen(2),tLen(2),cLen(2),oLen(2),fLen(2),sLen(2),eLen(2),aLen(2))); % this reshape works fine if the lengths are passed in this reverse order
            eB = squeeze(reshape(eB,rLen(3),pLen(3),tLen(3),cLen(3),oLen(3),fLen(3),sLen(3),eLen(3),aLen(3),rLen(2),pLen(2),tLen(2),cLen(2),oLen(2),fLen(2),sLen(2),eLen(2),aLen(2)));
            sB = squeeze(reshape(sB,rLen(3),pLen(3),tLen(3),cLen(3),oLen(3),fLen(3),sLen(3),eLen(3),aLen(3),rLen(2),pLen(2),tLen(2),cLen(2),oLen(2),fLen(2),sLen(2),eLen(2),aLen(2)));
            fB = squeeze(reshape(fB,rLen(3),pLen(3),tLen(3),cLen(3),oLen(3),fLen(3),sLen(3),eLen(3),aLen(3),rLen(2),pLen(2),tLen(2),cLen(2),oLen(2),fLen(2),sLen(2),eLen(2),aLen(2)));
            oB = squeeze(reshape(oB,rLen(3),pLen(3),tLen(3),cLen(3),oLen(3),fLen(3),sLen(3),eLen(3),aLen(3),rLen(2),pLen(2),tLen(2),cLen(2),oLen(2),fLen(2),sLen(2),eLen(2),aLen(2)));
            cB = squeeze(reshape(cB,rLen(3),pLen(3),tLen(3),cLen(3),oLen(3),fLen(3),sLen(3),eLen(3),aLen(3),rLen(2),pLen(2),tLen(2),cLen(2),oLen(2),fLen(2),sLen(2),eLen(2),aLen(2)));
            tB = squeeze(reshape(tB,rLen(3),pLen(3),tLen(3),cLen(3),oLen(3),fLen(3),sLen(3),eLen(3),aLen(3),rLen(2),pLen(2),tLen(2),cLen(2),oLen(2),fLen(2),sLen(2),eLen(2),aLen(2)));
            pB = squeeze(reshape(pB,rLen(3),pLen(3),tLen(3),cLen(3),oLen(3),fLen(3),sLen(3),eLen(3),aLen(3),rLen(2),pLen(2),tLen(2),cLen(2),oLen(2),fLen(2),sLen(2),eLen(2),aLen(2)));
            rB = squeeze(reshape(rB,rLen(3),pLen(3),tLen(3),cLen(3),oLen(3),fLen(3),sLen(3),eLen(3),aLen(3),rLen(2),pLen(2),tLen(2),cLen(2),oLen(2),fLen(2),sLen(2),eLen(2),aLen(2)));

            aR = squeeze(reshape(aR,rLen(3),pLen(3),tLen(3),cLen(3),oLen(3),fLen(3),sLen(3),eLen(3),aLen(3),rLen(2),pLen(2),tLen(2),cLen(2),oLen(2),fLen(2),sLen(2),eLen(2),aLen(2))); % this reshape works fine if the lengths are passed in this reverse order
            eR = squeeze(reshape(eR,rLen(3),pLen(3),tLen(3),cLen(3),oLen(3),fLen(3),sLen(3),eLen(3),aLen(3),rLen(2),pLen(2),tLen(2),cLen(2),oLen(2),fLen(2),sLen(2),eLen(2),aLen(2)));
            sR = squeeze(reshape(sR,rLen(3),pLen(3),tLen(3),cLen(3),oLen(3),fLen(3),sLen(3),eLen(3),aLen(3),rLen(2),pLen(2),tLen(2),cLen(2),oLen(2),fLen(2),sLen(2),eLen(2),aLen(2)));
            fR = squeeze(reshape(fR,rLen(3),pLen(3),tLen(3),cLen(3),oLen(3),fLen(3),sLen(3),eLen(3),aLen(3),rLen(2),pLen(2),tLen(2),cLen(2),oLen(2),fLen(2),sLen(2),eLen(2),aLen(2)));
            oR = squeeze(reshape(oR,rLen(3),pLen(3),tLen(3),cLen(3),oLen(3),fLen(3),sLen(3),eLen(3),aLen(3),rLen(2),pLen(2),tLen(2),cLen(2),oLen(2),fLen(2),sLen(2),eLen(2),aLen(2)));
            cR = squeeze(reshape(cR,rLen(3),pLen(3),tLen(3),cLen(3),oLen(3),fLen(3),sLen(3),eLen(3),aLen(3),rLen(2),pLen(2),tLen(2),cLen(2),oLen(2),fLen(2),sLen(2),eLen(2),aLen(2)));
            tR = squeeze(reshape(tR,rLen(3),pLen(3),tLen(3),cLen(3),oLen(3),fLen(3),sLen(3),eLen(3),aLen(3),rLen(2),pLen(2),tLen(2),cLen(2),oLen(2),fLen(2),sLen(2),eLen(2),aLen(2)));
            pR = squeeze(reshape(pR,rLen(3),pLen(3),tLen(3),cLen(3),oLen(3),fLen(3),sLen(3),eLen(3),aLen(3),rLen(2),pLen(2),tLen(2),cLen(2),oLen(2),fLen(2),sLen(2),eLen(2),aLen(2)));
            rR = squeeze(reshape(rR,rLen(3),pLen(3),tLen(3),cLen(3),oLen(3),fLen(3),sLen(3),eLen(3),aLen(3),rLen(2),pLen(2),tLen(2),cLen(2),oLen(2),fLen(2),sLen(2),eLen(2),aLen(2)));

            
            index = squeeze(reshape(index,rLen(3),pLen(3),tLen(3),cLen(3),oLen(3),fLen(3),sLen(3),eLen(3),aLen(3),rLen(2),pLen(2),tLen(2),cLen(2),oLen(2),fLen(2),sLen(2),eLen(2),aLen(2)));
            
            
        case 10

            for a = 1:aLen(1)
                for e = 1:eLen(1)
                    for s = 1:sLen(1)
                        for f = 1:fLen(1)
                            for o = 1:oLen(1)
                                for c = 1:cLen(1)
                                    for t = 1:tLen(1)
                                        for p = 1:pLen(1)
                                            for r = 1:rLen(1)
                                                
                                                for a2 = 1:aLen(2)
                                                    for e2 = 1:eLen(2)
                                                        for s2 = 1:sLen(2)
                                                            for f2 = 1:fLen(2)
                                                                for o2 = 1:oLen(2)
                                                                    for c2 = 1:cLen(2)
                                                                        for t2 = 1:tLen(2)
                                                                            for p2 = 1:pLen(2)
                                                                                
                                                                                for r3 = 1:rLen(3)
                                                                                    
                                                                                    for r2 = 1:rLen(2)

                                                                                            annWidth = rValsUnique{2}(r2) - rValsUnique{3}(r3);
                                                                                            
                                                                                            if ~isempty(find(validAnnulusWidth == annWidth, 1))
                                                                                            
                                                                                                aB = cat(1,aB,aValsUnique{1}(a));
                                                                                                eB = cat(1,eB,eValsUnique{1}(e));
                                                                                                sB = cat(1,sB,sValsUnique{1}(s));
                                                                                                fB = cat(1,fB,fValsUnique{1}(f));
                                                                                                oB = cat(1,oB,oValsUnique{1}(o));
                                                                                                cB = cat(1,cB,cValsUnique{1}(c));
                                                                                                tB = cat(1,tB,tValsUnique{1}(t));
                                                                                                pB = cat(1,pB,pValsUnique{1}(p));
                                                                                                rB{count} = [0 rValsUnique{1}(r)/rFactor];

                                                                                                aR = cat(1,aR,aValsUnique{2}(a2));
                                                                                                eR = cat(1,eR,eValsUnique{2}(e2));
                                                                                                sR = cat(1,sR,sValsUnique{2}(s2));
                                                                                                fR = cat(1,fR,fValsUnique{2}(f2));
                                                                                                oR = cat(1,oR,oValsUnique{2}(o2));
                                                                                                cR = cat(1,cR,cValsUnique{2}(c2));
                                                                                                tR = cat(1,tR,tValsUnique{2}(t2));
                                                                                                pR = cat(1,pR,pValsUnique{2}(p2));
                                                                                                rR{count} = [rValsUnique{3}(r3)/rFactor rValsUnique{2}(r2)/rFactor];


                                                                                                index{count} = [a,e,s,f,o,c,t,p,r,a2,e2,s2,f2,o2,c2,t2,p2,r2,r3];

                                                                                                count = count+1;
                                                                                            end

                                                                                            
                                                                                    end

                                                                                end

                                                                            end
                                                                        end
                                                                    end
                                                                end
                                                            end
                                                        end
                                                    end
                                                end
                                            end

                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
            
            rLen(2) = length(validAnnulusWidth);


            aB = (squeeze(reshape(aB,rLen(2),rLen(3),pLen(2),tLen(2),cLen(2),oLen(2),fLen(2),sLen(2),eLen(2),aLen(2),rLen(1),pLen(1),tLen(1),cLen(1),oLen(1),fLen(1),sLen(1),eLen(1),aLen(1))))'; % this reshape works fine if the lengths are passed in this reverse order
            eB = (squeeze(reshape(eB,rLen(2),rLen(3),pLen(2),tLen(2),cLen(2),oLen(2),fLen(2),sLen(2),eLen(2),aLen(2),rLen(1),pLen(1),tLen(1),cLen(1),oLen(1),fLen(1),sLen(1),eLen(1),aLen(1))))';
            sB = (squeeze(reshape(sB,rLen(2),rLen(3),pLen(2),tLen(2),cLen(2),oLen(2),fLen(2),sLen(2),eLen(2),aLen(2),rLen(1),pLen(1),tLen(1),cLen(1),oLen(1),fLen(1),sLen(1),eLen(1),aLen(1))))';
            fB = (squeeze(reshape(fB,rLen(2),rLen(3),pLen(2),tLen(2),cLen(2),oLen(2),fLen(2),sLen(2),eLen(2),aLen(2),rLen(1),pLen(1),tLen(1),cLen(1),oLen(1),fLen(1),sLen(1),eLen(1),aLen(1))))';
            oB = (squeeze(reshape(oB,rLen(2),rLen(3),pLen(2),tLen(2),cLen(2),oLen(2),fLen(2),sLen(2),eLen(2),aLen(2),rLen(1),pLen(1),tLen(1),cLen(1),oLen(1),fLen(1),sLen(1),eLen(1),aLen(1))))';
            cB = (squeeze(reshape(cB,rLen(2),rLen(3),pLen(2),tLen(2),cLen(2),oLen(2),fLen(2),sLen(2),eLen(2),aLen(2),rLen(1),pLen(1),tLen(1),cLen(1),oLen(1),fLen(1),sLen(1),eLen(1),aLen(1))))';
            tB = (squeeze(reshape(tB,rLen(2),rLen(3),pLen(2),tLen(2),cLen(2),oLen(2),fLen(2),sLen(2),eLen(2),aLen(2),rLen(1),pLen(1),tLen(1),cLen(1),oLen(1),fLen(1),sLen(1),eLen(1),aLen(1))))';
            pB = (squeeze(reshape(pB,rLen(2),rLen(3),pLen(2),tLen(2),cLen(2),oLen(2),fLen(2),sLen(2),eLen(2),aLen(2),rLen(1),pLen(1),tLen(1),cLen(1),oLen(1),fLen(1),sLen(1),eLen(1),aLen(1))))';
            rB = (squeeze(reshape(rB,rLen(2),rLen(3),pLen(2),tLen(2),cLen(2),oLen(2),fLen(2),sLen(2),eLen(2),aLen(2),rLen(1),pLen(1),tLen(1),cLen(1),oLen(1),fLen(1),sLen(1),eLen(1),aLen(1))))';

            aR = (squeeze(reshape(aR,rLen(2),rLen(3),pLen(2),tLen(2),cLen(2),oLen(2),fLen(2),sLen(2),eLen(2),aLen(2),rLen(1),pLen(1),tLen(1),cLen(1),oLen(1),fLen(1),sLen(1),eLen(1),aLen(1))))'; % this reshape works fine if the lengths are passed in this reverse order
            eR = (squeeze(reshape(eR,rLen(2),rLen(3),pLen(2),tLen(2),cLen(2),oLen(2),fLen(2),sLen(2),eLen(2),aLen(2),rLen(1),pLen(1),tLen(1),cLen(1),oLen(1),fLen(1),sLen(1),eLen(1),aLen(1))))';
            sR = (squeeze(reshape(sR,rLen(2),rLen(3),pLen(2),tLen(2),cLen(2),oLen(2),fLen(2),sLen(2),eLen(2),aLen(2),rLen(1),pLen(1),tLen(1),cLen(1),oLen(1),fLen(1),sLen(1),eLen(1),aLen(1))))';
            fR = (squeeze(reshape(fR,rLen(2),rLen(3),pLen(2),tLen(2),cLen(2),oLen(2),fLen(2),sLen(2),eLen(2),aLen(2),rLen(1),pLen(1),tLen(1),cLen(1),oLen(1),fLen(1),sLen(1),eLen(1),aLen(1))))';
            oR = (squeeze(reshape(oR,rLen(2),rLen(3),pLen(2),tLen(2),cLen(2),oLen(2),fLen(2),sLen(2),eLen(2),aLen(2),rLen(1),pLen(1),tLen(1),cLen(1),oLen(1),fLen(1),sLen(1),eLen(1),aLen(1))))';
            cR = (squeeze(reshape(cR,rLen(2),rLen(3),pLen(2),tLen(2),cLen(2),oLen(2),fLen(2),sLen(2),eLen(2),aLen(2),rLen(1),pLen(1),tLen(1),cLen(1),oLen(1),fLen(1),sLen(1),eLen(1),aLen(1))))';
            tR = (squeeze(reshape(tR,rLen(2),rLen(3),pLen(2),tLen(2),cLen(2),oLen(2),fLen(2),sLen(2),eLen(2),aLen(2),rLen(1),pLen(1),tLen(1),cLen(1),oLen(1),fLen(1),sLen(1),eLen(1),aLen(1))))';
            pR = (squeeze(reshape(pR,rLen(2),rLen(3),pLen(2),tLen(2),cLen(2),oLen(2),fLen(2),sLen(2),eLen(2),aLen(2),rLen(1),pLen(1),tLen(1),cLen(1),oLen(1),fLen(1),sLen(1),eLen(1),aLen(1))))';
            rR = (squeeze(reshape(rR,rLen(2),rLen(3),pLen(2),tLen(2),cLen(2),oLen(2),fLen(2),sLen(2),eLen(2),aLen(2),rLen(1),pLen(1),tLen(1),cLen(1),oLen(1),fLen(1),sLen(1),eLen(1),aLen(1))))';

            
            index = (squeeze(reshape(index,rLen(2),rLen(3),pLen(2),tLen(2),cLen(2),oLen(2),fLen(2),sLen(2),eLen(2),aLen(2),rLen(1),pLen(1),tLen(1),cLen(1),oLen(1),fLen(1),sLen(1),eLen(1),aLen(1))))';
            
            
    end
      
    figure;
    
    numDim = sum(size(index)~=1);
    if numDim == 1
        for i = 1:size(index,1)
            gaborBackground.azimuthDeg = aB(i);
            gaborBackground.elevationDeg = eB(i);
            gaborBackground.sigmaDeg = sB(i);
            gaborBackground.spatialFreqCPD = fB(i);
            gaborBackground.orientationDeg = oB(i);
            gaborBackground.contrastPC = cB(i);
            gaborBackground.temporalFreqHz = tB(i);
            gaborBackground.spatialPhaseDeg = pB(i);
            gaborBackground.radiusDeg = rB{i};
            
            gaborRing.azimuthDeg = aR(i);
            gaborRing.elevationDeg = eR(i);
            gaborRing.sigmaDeg = sR(i);
            gaborRing.spatialFreqCPD = fR(i);
            gaborRing.orientationDeg = oR(i);
            gaborRing.contrastPC = cR(i);
            gaborRing.temporalFreqHz = tR(i);
            gaborRing.spatialPhaseDeg = pR(i);
            gaborRing.radiusDeg = rR{i};
            
            hS{i} = subplot(1,size(index,1),i);
            titleString = getTitleCRS(gaborBackground,gaborRing,protocolNumber,rFactor);
            drawStimulus(hS{i},gaborBackground,gaborRing,titleString,protocolNumber);
        end
        
    elseif numDim == 2
        for i = 1:size(index,1)
            for j = 1:size(index,2)
                gaborBackground.azimuthDeg = aB(i,j);
                gaborBackground.elevationDeg = eB(i,j);
                gaborBackground.sigmaDeg = sB(i,j);
                gaborBackground.spatialFreqCPD = fB(i,j);
                gaborBackground.orientationDeg = oB(i,j);
                gaborBackground.contrastPC = cB(i,j);
                gaborBackground.temporalFreqHz = tB(i,j);
                gaborBackground.spatialPhaseDeg = pB(i,j);
                gaborBackground.radiusDeg = rB{i,j};
                
                gaborRing.azimuthDeg = aR(i,j);
                gaborRing.elevationDeg = eR(i,j);
                gaborRing.sigmaDeg = sR(i,j);
                gaborRing.spatialFreqCPD = fR(i,j);
                gaborRing.orientationDeg = oR(i,j);
                gaborRing.contrastPC = cR(i,j);
                gaborRing.temporalFreqHz = tR(i,j);
                gaborRing.spatialPhaseDeg = pR(i,j);
                gaborRing.radiusDeg = rR{i,j};
                
                numRows = size(index,1);
                numCols = size(index,2);
                hS{i,j} = subplot(numRows,numCols,(i-1)*numCols+j);
                [~,xString,yString] = getTitleCRS(gaborBackground,gaborRing,protocolNumber,rFactor);
                drawStimulus(hS{i,j},gaborBackground,gaborRing,[],protocolNumber);
                
                if i==1 % top row
                    title(xString,'Fontsize',fontSizeTiny,'Parent',hS{i,j});
                end
                if j==1 % first column
                    ylabel(yString,'Fontsize',fontSizeTiny,'Rotation',0,'Parent',hS{i,j});
                end
            end
        end
        
    elseif numDim == 3
        for i = 1:size(index,1)
            for j = 1:size(index,2)
                for k = 1:size(index,3)
                    gaborBackground.azimuthDeg = aB(i,j,k);
                    gaborBackground.elevationDeg = eB(i,j,k);
                    gaborBackground.sigmaDeg = sB(i,j,k);
                    gaborBackground.spatialFreqCPD = fB(i,j,k);
                    gaborBackground.orientationDeg = oB(i,j,k);
                    gaborBackground.contrastPC = cB(i,j,k);
                    gaborBackground.temporalFreqHz = tB(i,j,k);
                    gaborBackground.spatialPhaseDeg = pB(i,j,k);
                    gaborBackground.radiusDeg = rB{i,j,k};

                    gaborRing.azimuthDeg = aR(i,j,k);
                    gaborRing.elevationDeg = eR(i,j,k);
                    gaborRing.sigmaDeg = sR(i,j,k);
                    gaborRing.spatialFreqCPD = fR(i,j,k);
                    gaborRing.orientationDeg = oR(i,j,k);
                    gaborRing.contrastPC = cR(i,j,k);
                    gaborRing.temporalFreqHz = tR(i,j,k);
                    gaborRing.spatialPhaseDeg = pR(i,j,k);
                    gaborRing.radiusDeg = rR{i,j,k};

                    numRows = size(index,1);
                    numCols = size(index,2);
                    numRepeats = size(index,3);
                    hS{i,j,k} = subplot(numRows*numRepeats,numCols,(k-1)*numRows*numCols+(i-1)*numCols+j);
                    titleString = getTitleCRS(gaborBackground,gaborRing,protocolNumber,rFactor);
                    drawStimulus(hS{i,j,k},gaborBackground,gaborRing,titleString,protocolNumber);
                    axis off;
                    axis tight;
                end
            end
        end
    end

end

end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Supporting Functions
%--------------------------------------------------------------------------


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

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function protocolNumber = getProtocolNumber(folderExtract)
load (fullfile(folderExtract,'stimResults'));
protocolNumber = stimResults.protocolNumber;
end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function gaborsDisplayed = getGaborsDisplayed(folderExtract)
load (fullfile(folderExtract,'stimResults'));
gaborsDisplayed = stimResults.side;
end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function drawStimulus(h,gaborBackground,gaborRing,titleString,protocolNumber)

if ~exist('gaborRing','var')
    gaborRing.azimuthDeg = 0;
    gaborRing.elevationDeg = 0;
    gaborRing.sigmaDeg = 0;
    gaborRing.spatialFreqCPD = 0;
    gaborRing.orientationDeg = 0;
    gaborRing.contrastPC = 0;
    gaborRing.temporalFreqHz = 0;
    gaborRing.spatialPhaseDeg = 0;
    gaborRing.radiusDeg = 0;
end

if ~exist('protocolNumber','var')
    protocolNumber = 11;
end

% gridLims = [-3 -0.5 -2.5 0];
gridLims = [-6 -1.5 -4 0.5];
% gridLims = [-21.5 -1.5 -19.5 0.5];
gridLimsNormalized(1) = -(gridLims(2)-gridLims(1))/2;
gridLimsNormalized(2) = (gridLims(2)-gridLims(1))/2;
gridLimsNormalized(3) = -(gridLims(4)-gridLims(3))/2;
gridLimsNormalized(4) = (gridLims(4)-gridLims(3))/2;

aPoints=gridLimsNormalized(1):1/30:gridLimsNormalized(2);
ePoints=gridLimsNormalized(3):1/30:gridLimsNormalized(4);

aVals = aPoints;
eVals = ePoints;

if protocolNumber == 3 || protocolNumber == 4 || protocolNumber == 5 % dual protocols
    innerphase = gaborBackground.spatialPhaseDeg; % for DPP
else
    innerphase = gaborRing.spatialPhaseDeg; % for PRP
end

gaborPatch = makeGRGStimulusWithPhase(gaborRing,gaborBackground,aVals,eVals,innerphase);

% Changed to show absolute contrasts and not normalized values - 26 Jan'13 
gabP = gaborPatch./100;
hGab = imshow(gabP,'Parent',h); %colorbar;
fontSizeTiny = 7;
title(titleString,'Fontsize',fontSizeTiny);

end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function [titleString,xString,yString] = getTitleGRF(gaborBackground)

sf = gaborBackground.spatialFreqCPD;
o = gaborBackground.orientationDeg;
titleString = ['SF = ' num2str(sf) ' cyc/deg; Ori = ' num2str(o) ' deg'];
xString = ['SF = ' num2str(sf) ' cyc/deg'];
yString = ['Ori = ' num2str(o) ' deg'];

end

function [titleString,xString,yString] = getTitleCRS(gaborBackground,gaborRing,protocolNumber,rFactor)

switch protocolNumber
    case 1 % RP
        titleString = ['R rad = ' num2str(gaborRing.radiusDeg(2)) 'deg; C/S cont = ' num2str(gaborBackground.contrastPC) ' %'];
        xString = ['R rad = ' num2str(gaborRing.radiusDeg) ' deg'];
        yString = ['C/S cont = ' num2str(gaborBackground.contrastPC) ' %'];
    case 2 % CRP
        titleString = ['C/S cont = ' num2str(gaborBackground.contrastPC) ' %; R cont = ' num2str(gaborRing.contrastPC) ' %'];
        yString = ['R cont = ' num2str(gaborRing.contrastPC) ' %'];
        xString = ['C/S cont = ' num2str(gaborBackground.contrastPC) ' %'];
    case 3 % DCP
%         titleString = ['C cont = ' num2str(gaborBackground.contrastPC) ' %; R cont = ' num2str(gaborRing.contrastPC) ' %'];
        titleString = ['C rad = ' num2str(gaborRing.radiusDeg(2)*rFactor) ' deg'];
        xString = ['R cont = ' num2str(gaborBackground.contrastPC) ' %'];
        yString = ['C cont = ' num2str(gaborRing.contrastPC) ' %'];
    case 4 % DOP
        titleString = ['C ori = ' num2str(gaborRing.orientationDeg) ' deg; R ori = ' num2str(gaborBackground.orientationDeg) ' deg'];
        xString = ['R ori = ' num2str(gaborBackground.orientationDeg) ' deg'];
        yString = ['C ori = ' num2str(gaborRing.orientationDeg) ' deg'];
    case 5 % DPP
        titleString = ['C phase = ' num2str(gaborRing.spatialPhaseDeg) ' deg; R phase = ' num2str(gaborBackground.spatialPhaseDeg) ' deg'];
        xString = ['R phase = ' num2str(gaborBackground.spatialPhaseDeg) ' deg'];
        yString = ['C phase = ' num2str(gaborRing.spatialPhaseDeg) ' deg'];
    case 6 % ORP
        titleString = ['C/S ori = ' num2str(gaborBackground.orientationDeg) ' deg; R ori = ' num2str(gaborRing.orientationDeg) ' deg'];
        xString = ['R ori = ' num2str(gaborBackground.orientationDeg) ' deg'];
        yString = ['C/S ori = ' num2str(gaborRing.orientationDeg) ' deg'];
    case 7 % PRP
        titleString = ['C/S phase = ' num2str(gaborBackground.spatialPhaseDeg) ' deg; R phase = ' num2str(gaborRing.spatialPhaseDeg) ' deg'];
        xString = ['R phase = ' num2str(gaborBackground.spatialPhaseDeg) ' deg'];
        yString = ['C/S phase = ' num2str(gaborRing.spatialPhaseDeg) ' deg'];
    case 8 % DRP
        titleString = ['C/S TF = ' num2str(gaborBackground.temporalFreqHz) ' deg/s; R TF = ' num2str(gaborRing.temporalFreqHz) ' deg/s'];
        xString = ['R TF = ' num2str(gaborBackground.temporalFreqHz) ' deg/s'];
        yString = ['C/S TF = ' num2str(gaborRing.temporalFreqHz) ' deg/s'];
    case 9 % COS
        titleString = ['C ori = ' num2str(gaborBackground.orientationDeg) ' deg; R ori = ' num2str(gaborRing.orientationDeg) ' deg; C TF = ' num2str(gaborBackground.temporalFreqHz) ' deg/s; R TF = ' num2str(gaborRing.temporalFreqHz) ' deg/s'];
    case 10 % AFP
        titleString = ['C rad = ' num2str(gaborRing.radiusDeg(1)*rFactor) 'deg; Ann width = ' num2str((gaborRing.radiusDeg(2)-gaborRing.radiusDeg(1))*rFactor) ' deg'];
        xString = ['Ann width = ' num2str((gaborRing.radiusDeg(2)-gaborRing.radiusDeg(1))*rFactor) ' deg'];
        yString = ['C rad = ' num2str(gaborRing.radiusDeg(1)*rFactor) 'deg'];
        
end
        
end
