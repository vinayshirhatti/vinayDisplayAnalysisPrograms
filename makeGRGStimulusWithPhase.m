% makeGRGStimulus(gabor0,gabor1,aVals,eVals)
% GRG stands for Grating-RING, composed of a gabor/grating in the center
% (gabor0), and the annulus gabor1. 

function gaborPatch = makeGRGStimulusWithPhase(gabor0,gabor1,aVals,eVals,innerphase)

[gaborPatch0,aperature0] = makeGaborStimulusWithPhase(gabor0,aVals,eVals,0,innerphase);
gaborPatch1 = makeGaborStimulus(gabor1,aVals,eVals);

gaborPatch = gaborPatch1;
changeTheseValues=find(aperature0==1);
gaborPatch(changeTheseValues) = gaborPatch0(changeTheseValues);
end