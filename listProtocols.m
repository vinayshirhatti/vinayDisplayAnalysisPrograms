% List of all protocols
% Vinay Shirhatti, 24 March 2015
% =========================================================================

function [listProtocolName,subjectNames,expDates,protocolNames,stimTypes,gridLayout] = listProtocols

% StimTypes
% timeStartFromBaseLineList(1) = -0.55; deltaTList(1) = 1.024; % in seconds
% timeStartFromBaseLineList(2) = -1.148; deltaTList(2) = 2.048;
% timeStartFromBaseLineList(3) = -1.596; deltaTList(3) = 4.096; % [Vinay] - for a longer stim time; Stim ON and OFF 1500ms each
% timeStartFromBaseLineList(4) = -0.848; deltaTList(4) = 2.048; % [Vinay] - 800ms stim, 700ms interstim
% timeStartFromBaseLineList(5) = -1.096; deltaTList(5) = 4.096; % [Vinay] - for a longer stim time

% 230315
listProtocolName{1} = 'GRF_001#photodiode230315'; % ainp3/4/5 => photodiode1/2/3; 6 contrasts, 6 TFs
expDates{1} = '230315'; protocolNames{1} = 'GRF_001'; stimTypes{1} = 4; % 800ms ON, 700ms OFF
subjectNames{1} = 'photodiode';
gridLayout{1} = 10; % Layout 10 for only ainp

listProtocolName{2} = 'CRS_001#AFP#photodiode230315'; % ainp3/4/5 => photodiode1/2/3
% AFP, S,C - 6 contrasts, 6 TFs
expDates{2} = '230315'; protocolNames{2} = 'CRS_001'; stimTypes{2} = 4; % 800ms ON, 700ms OFF
subjectNames{2} = 'photodiode';
gridLayout{2} = 10; % Layout 10 for only ainp

listProtocolName{3} = 'CRS_002#DCP#photodiode230315'; % ainp3/4/5 => photodiode1/2/3
% DCP, C,R - 5 contrasts (independent), 5 TFs (matched)
expDates{3} = '230315'; protocolNames{3} = 'CRS_002'; stimTypes{3} = 4; % 800ms ON, 700ms OFF
subjectNames{3} = 'photodiode';
gridLayout{3} = 10; % Layout 10 for only ainp

% 190315
listProtocolName{4} = 'GRF_001#ST190315'; % ST; 6 SF, 4 ori
expDates{4} = '190315'; protocolNames{4} = 'GRF_001'; stimTypes{4} = 4; % 800ms ON, 700ms OFF
subjectNames{4} = 'ST';
gridLayout{4} = 1; % 64 electrodes


% 240315
listProtocolName{5} = 'GRF_001#SFOri#GM240315'; % ST; 6 SF, 4 ori
expDates{5} = '240315'; protocolNames{5} = 'GRF_001'; stimTypes{5} = 4; % 800ms ON, 700ms OFF
subjectNames{5} = 'GM';
% 64 channels; ainp 1/2 => eyeX/Y, ainp 3/4/5 => photodiode1/2/3
% SF, Ori Tuning; 6 SFs (0.5 to 16), 4 orientations (0,45,90,135)
gridLayout{5} = 1; % 64 electrodes

listProtocolName{6} = 'CRS_001#Size#GM240315';
expDates{6} = '240315'; protocolNames{6} = 'CRS_001'; stimTypes{6} = 4; % 800ms ON, 700ms OFF
subjectNames{6} = 'GM';
% 64 channels; ainp 1/2 => eyeX/Y, ainp 3/4/5 => photodiode1/2/3
% Size tuning; using DCP; centre radius - 0,1,2,4,8,16,32 (full screen) deg
% Ring at 0% contrast
gridLayout{6} = 1; % 64 electrodes

listProtocolName{7} = 'CRS_002#AFP#GM240315';
expDates{7} = '240315'; protocolNames{7} = 'CRS_002'; stimTypes{7} = 4; % 800ms ON, 700ms OFF
subjectNames{7} = 'GM';
% 64 channels; ainp 1/2 => eyeX/Y, ainp 3/4/5 => photodiode1/2/3
% Annulus Fixed Protocol; centre radius: 1,2,4,8;
% Ring radius i.e. annulus width: 0.25, 0.5, 1, 2, 4
gridLayout{7} = 1; % 64 electrodes


% 250315
listProtocolName{8} = 'CRS_001#Size#GM250315';
expDates{8} = '250315'; protocolNames{8} = 'CRS_001'; stimTypes{8} = 4; % 800ms ON, 700ms OFF
subjectNames{8} = 'GM';
% 64 channels; ainp 1/2 => eyeX/Y, ainp 3/4/5 => photodiode1/2/3
% Size tuning; using DCP; centre radius - 0,1,2,4,8,16,32 (full screen) deg
% Ring at 0% contrast
gridLayout{8} = 1; % 64 electrodes

listProtocolName{9} = 'CRS_002#AFP#GM250315';
expDates{9} = '250315'; protocolNames{9} = 'CRS_002'; stimTypes{9} = 4; % 800ms ON, 700ms OFF
subjectNames{9} = 'GM';
% 64 channels; ainp 1/2 => eyeX/Y, ainp 3/4/5 => photodiode1/2/3
% Annulus Fixed Protocol; centre radius: 1,2,4,8;
% Ring radius i.e. annulus width: 0.25, 0.5, 1, 2, 4
gridLayout{9} = 1; % 64 electrodes

listProtocolName{10} = 'CRS_003#DCP#GM250315';
expDates{10} = '250315'; protocolNames{10} = 'CRS_003'; stimTypes{10} = 4; % 800ms ON, 700ms OFF
subjectNames{10} = 'GM';
% 64 channels; ainp 1/2 => eyeX/Y, ainp 3/4/5 => photodiode1/2/3
% Dual Contrast Protocol; centre radius: 2 deg;
% Centre & Ring contrast: 4x4 contrast. Contrast mapping was linear
% (0,33.4,66.6,100)
gridLayout{10} = 1; % 64 electrodes


% 250315
listProtocolName{11} = 'GRF_001#SFOri#AB270315'; % ST; 6 SF, 4 ori
expDates{11} = '270315'; protocolNames{11} = 'GRF_001'; stimTypes{11} = 4; % 800ms ON, 700ms OFF
subjectNames{11} = 'AB';
% 64 channels; ainp 1/2 => eyeX/Y, ainp 3/4/5 => photodiode1/2/3
% SF, Ori Tuning; 6 SFs (0.5 to 16), 4 orientations (0,45,90,135)
gridLayout{11} = 1; % 64 electrodes

listProtocolName{12} = 'GRF_002#SFOri#AB270315'; % ST; 6 SF, 4 ori
expDates{12} = '270315'; protocolNames{12} = 'GRF_002'; stimTypes{12} = 4; % 800ms ON, 700ms OFF
subjectNames{12} = 'AB';
% 64 channels; ainp 1/2 => eyeX/Y, ainp 3/4/5 => photodiode1/2/3
% SF, Ori Tuning; 6 SFs (0.5 to 16), 4 orientations (0,45,90,135)
gridLayout{12} = 1; % 64 electrodes

listProtocolName{13} = 'CRS_001#Size#AB270315';
expDates{13} = '270315'; protocolNames{13} = 'CRS_001'; stimTypes{13} = 4; % 800ms ON, 700ms OFF
subjectNames{13} = 'AB';
% 64 channels; ainp 1/2 => eyeX/Y, ainp 3/4/5 => photodiode1/2/3
% Size tuning; using DCP; centre radius - 0,1,2,4,8,16,32 (full screen) deg
% Ring at 0% contrast
gridLayout{13} = 1; % 64 electrodes

listProtocolName{14} = 'CRS_002#AFP#AB270315';
expDates{14} = '270315'; protocolNames{14} = 'CRS_002'; stimTypes{14} = 4; % 800ms ON, 700ms OFF
subjectNames{14} = 'AB';
% 64 channels; ainp 1/2 => eyeX/Y, ainp 3/4/5 => photodiode1/2/3
% Annulus Fixed Protocol; centre radius: 1,2,4,8;
% Ring radius i.e. annulus width: 0.25, 0.5, 1, 2, 4
gridLayout{14} = 1; % 64 electrodes


% % 261014
% listProtocolName{15} = 'CRS_001#murty261014_sizeTuning';
% expDates{15} = '261014'; protocolNames{15} = 'CRS_001'; stimTypes{15} = 4; % Murty; EEG Electrodes: 10-20 system 21 electrodes
% % CRS (DCP) protocol for size tuning;  7 radii - 0,1,2,4,8,16,32 deg
% subjectNames{15} = 'murty';
% gridLayout{15} = 3; % 21 electrodes


% 310315
subjectNames{15} = 'AD'; expDates{15} = '310315'; protocolNames{15} = 'GRF_001'; stimTypes{15} = 4;
% 64 channels; ainp 1/2 => eyeX/Y, ainp 3/4/5 => photodiode1/2/3
% SF, Ori Tuning; 6 SFs (0.5 to 16), 4 orientations (0,45,90,135)
gridLayout{15} = 1; % 64 electrodes

subjectNames{16} = 'AD'; expDates{16} = '310315'; protocolNames{16} = 'GRF_002'; stimTypes{16} = 4;
% 64 channels; ainp 1/2 => eyeX/Y, ainp 3/4/5 => photodiode1/2/3
% SF, Ori Tuning; 6 SFs (0.5 to 16), 4 orientations (0,45,90,135)
% The fixation spot was 0.1 deg in this case, as against 0.05 in the
% previous
gridLayout{16} = 1; % 64 electrodes

subjectNames{17} = 'AD'; expDates{17} = '310315'; protocolNames{17} = 'CRS_001'; stimTypes{17} = 4;
% 64 channels; ainp 1/2 => eyeX/Y, ainp 3/4/5 => photodiode1/2/3
% Size tuning; using DCP; centre radius - 0,1,2,4,8,16,32 (full screen) deg
% Ring at 0% contrast. Fixation: 0.05 deg; Ori = 45 deg
gridLayout{17} = 1; % 64 electrodes

subjectNames{18} = 'AD'; expDates{18} = '310315'; protocolNames{18} = 'CRS_002'; stimTypes{18} = 4;
% 64 channels; ainp 1/2 => eyeX/Y, ainp 3/4/5 => photodiode1/2/3
% Size tuning; using DCP; centre radius - 0,1,2,4,8,16,32 (full screen) deg
% Ring at 0% contrast. Fixation: 0.05 deg; Ori = 135 deg
gridLayout{18} = 1; % 64 electrodes

subjectNames{19} = 'AD'; expDates{19} = '310315'; protocolNames{19} = 'CRS_003'; stimTypes{19} = 4;
% 64 channels; ainp 1/2 => eyeX/Y, ainp 3/4/5 => photodiode1/2/3
% Annulus Fixed Protocol; centre radius: 1,2,4,8;
% Ring radius i.e. annulus width: 0,0.25, 0.5, 1, 2, 4; Ori = 135 deg
gridLayout{19} = 1; % 64 electrodes

% 030415
subjectNames{20} = 'AB'; expDates{20} = '030415'; protocolNames{20} = 'CRS_001'; stimType{20} = 4;
% 64 channels; ainp 1/2 => eyeX/Y, ainp 3/4/5 => photodiode1/2/3
% Dual Contrast Protocol; 4 x 4 contrasts (C and R)
% Centre radius = 2 deg. Fixation spot = 0.05 deg for all protocols on this
% day. Optimal SF = 2 cyc/deg (as per bipolar responses in SF tuning);
% optimal ori = 45 deg (as per both single and bipolar responses in Ori
% tuning
gridLayout{20} = 1; % 64 electrodes

subjectNames{21} = 'AB'; expDates{21} = '030415'; protocolNames{21} = 'CRS_002'; stimType{21} = 4;
% 64 channels; ainp 1/2 => eyeX/Y, ainp 3/4/5 => photodiode1/2/3
% Dual Contrast Protocol; 4 x 4 contrasts (C and R)
% Centre radius = 4 deg.
gridLayout{21} = 1; % 64 electrodes

subjectNames{22} = 'AB'; expDates{22} = '030415'; protocolNames{22} = 'CRS_003'; stimType{22} = 4;
% 64 channels; ainp 1/2 => eyeX/Y, ainp 3/4/5 => photodiode1/2/3
% Contrast Ring Protocol; 4 x 4 contrasts (C and R)
% Centre radius = 2 deg, Ring radius = 4 deg => annulus of 2 deg
gridLayout{22} = 1; % 64 electrodes

subjectNames{23} = 'AB'; expDates{23} = '030415'; protocolNames{23} = 'CRS_004'; stimType{23} = 4;
% 64 channels; ainp 1/2 => eyeX/Y, ainp 3/4/5 => photodiode1/2/3
gridLayout{23} = 1; % 64 electrodes

subjectNames{24} = 'AB'; expDates{24} = '030415'; protocolNames{24} = 'CRS_005'; stimType{24} = 4;
% 64 channels; ainp 1/2 => eyeX/Y, ainp 3/4/5 => photodiode1/2/3
gridLayout{24} = 1; % 64 electrodes

subjectNames{25} = 'AB'; expDates{25} = '030415'; protocolNames{25} = 'CRS_006'; stimType{25} = 4;
% 64 channels; ainp 1/2 => eyeX/Y, ainp 3/4/5 => photodiode1/2/3
gridLayout{25} = 1; % 64 electrodes

subjectNames{26} = 'AB'; expDates{26} = '030415'; protocolNames{26} = 'CRS_007'; stimType{26} = 4;
% 64 channels; ainp 1/2 => eyeX/Y, ainp 3/4/5 => photodiode1/2/3
gridLayout{26} = 1; % 64 electrodes

subjectNames{27} = 'AB'; expDates{27} = '030415'; protocolNames{27} = 'CRS_008'; stimType{27} = 4;
% 64 channels; ainp 1/2 => eyeX/Y, ainp 3/4/5 => photodiode1/2/3
gridLayout{27} = 1; % 64 electrodes

subjectNames{28} = 'AB'; expDates{28} = '030415'; protocolNames{28} = 'CRS_009'; stimType{28} = 4;
% 64 channels; ainp 1/2 => eyeX/Y, ainp 3/4/5 => photodiode1/2/3
gridLayout{28} = 1; % 64 electrodes



% 040415
subjectNames{29} = 'SB'; expDates{29} = '040415'; protocolNames{29} = 'GRF_001'; stimType{29} = 4;
% 64 channels; ainp 1/2 => eyeX/Y, ainp 3/4/5 => photodiode1/2/3
gridLayout{29} = 1; % 64 electrodes

subjectNames{30} = 'SB'; expDates{30} = '040415'; protocolNames{30} = 'CRS_001'; stimType{30} = 4;
% 64 channels; ainp 1/2 => eyeX/Y, ainp 3/4/5 => photodiode1/2/3
gridLayout{30} = 1; % 64 electrodes

subjectNames{31} = 'SB'; expDates{31} = '040415'; protocolNames{31} = 'CRS_002'; stimType{31} = 4;
% 64 channels; ainp 1/2 => eyeX/Y, ainp 3/4/5 => photodiode1/2/3
gridLayout{31} = 1; % 64 electrodes

subjectNames{32} = 'SB'; expDates{32} = '040415'; protocolNames{32} = 'CRS_003'; stimType{32} = 4;
% 64 channels; ainp 1/2 => eyeX/Y, ainp 3/4/5 => photodiode1/2/3
gridLayout{32} = 1; % 64 electrodes

subjectNames{33} = 'SB'; expDates{33} = '040415'; protocolNames{33} = 'CRS_004'; stimType{33} = 4;
% 64 channels; ainp 1/2 => eyeX/Y, ainp 3/4/5 => photodiode1/2/3
gridLayout{33} = 1; % 64 electrodes


% 080415
subjectNames{34} = 'AD'; expDates{34} = '080415'; protocolNames{34} = 'CRS_001'; stimType{34} = 4;
% 64 channels; ainp 1/2 => eyeX/Y, ainp 3/4/5 => photodiode1/2/3
gridLayout{34} = 1; % 64 electrodes

subjectNames{35} = 'AD'; expDates{35} = '080415'; protocolNames{35} = 'CRS_002'; stimType{35} = 4;
% 64 channels; ainp 1/2 => eyeX/Y, ainp 3/4/5 => photodiode1/2/3
gridLayout{35} = 1; % 64 electrodes

subjectNames{36} = 'AD'; expDates{36} = '080415'; protocolNames{36} = 'CRS_003'; stimType{36} = 4;
% 64 channels; ainp 1/2 => eyeX/Y, ainp 3/4/5 => photodiode1/2/3
gridLayout{36} = 1; % 64 electrodes

subjectNames{37} = 'AD'; expDates{37} = '080415'; protocolNames{37} = 'CRS_004'; stimType{37} = 4;
% 64 channels; ainp 1/2 => eyeX/Y, ainp 3/4/5 => photodiode1/2/3
gridLayout{37} = 1; % 64 electrodes

subjectNames{38} = 'AD'; expDates{38} = '080415'; protocolNames{38} = 'CRS_005'; stimType{38} = 4;
% 64 channels; ainp 1/2 => eyeX/Y, ainp 3/4/5 => photodiode1/2/3
gridLayout{38} = 1; % 64 electrodes

subjectNames{39} = 'AD'; expDates{39} = '080415'; protocolNames{39} = 'CRS_006'; stimType{39} = 4;
% 64 channels; ainp 1/2 => eyeX/Y, ainp 3/4/5 => photodiode1/2/3
gridLayout{39} = 1; % 64 electrodes

subjectNames{40} = 'AD'; expDates{40} = '080415'; protocolNames{40} = 'CRS_007'; stimType{40} = 4;
% 64 channels; ainp 1/2 => eyeX/Y, ainp 3/4/5 => photodiode1/2/3
gridLayout{40} = 1; % 64 electrodes

subjectNames{41} = 'AD'; expDates{41} = '080415'; protocolNames{41} = 'CRS_008'; stimType{41} = 4;
% 64 channels; ainp 1/2 => eyeX/Y, ainp 3/4/5 => photodiode1/2/3
gridLayout{41} = 1; % 64 electrodes

subjectNames{42} = 'AD'; expDates{42} = '080415'; protocolNames{42} = 'CRS_009'; stimType{42} = 4;
% 64 channels; ainp 1/2 => eyeX/Y, ainp 3/4/5 => photodiode1/2/3
gridLayout{42} = 1; % 64 electrodes


% 090415
subjectNames{43} = 'AP'; expDates{43} = '090415'; protocolNames{43} = 'GRF_002'; stimType{43} = 4;
% 64 channels; ainp 1/2 => eyeX/Y, ainp 3/4/5 => photodiode1/2/3
% SF, Ori Tuning; 6 SFs (0.5 to 16), 4 orientations (0,45,90,135)
gridLayout{43} = 1; % 64 electrodes


% 130415
subjectNames{44} = 'GM'; expDates{44} = '130415'; protocolNames{44} = 'GRF_001'; stimType{44} = 4;
% 64 channels; ainp 1/2 => eyeX/Y, ainp 3/4/5 => photodiode1/2/3
% SF, Ori Tuning; 6 SFs (0.5 to 16), 4 orientations (0,45,90,135)
gridLayout{44} = 1; % 64 electrodes

subjectNames{45} = 'GM'; expDates{45} = '130415'; protocolNames{45} = 'CRS_002'; stimType{44} = 4;
% 64 channels; ainp 1/2 => eyeX/Y, ainp 3/4/5 => photodiode1/2/3
% Dual Contrast Protocol, 4x4 C, R contrasts - 0,25,50,100%; C rad = 2deg
% Optimal SF = 2 cyc/deg; Ori = 135deg
gridLayout{45} = 1; % 64 electrodes

subjectNames{46} = 'GM'; expDates{46} = '130415'; protocolNames{46} = 'CRS_003'; stimType{46} = 4;
% 64 channels; ainp 1/2 => eyeX/Y, ainp 3/4/5 => photodiode1/2/3
% Dual Contrast Protocol, 4x4 C, R contrasts - 0,25,50,100%; C rad = 4deg
gridLayout{46} = 1; % 64 electrodes

subjectNames{47} = 'GM'; expDates{47} = '130415'; protocolNames{47} = 'CRS_004'; stimType{47} = 4;
% 64 channels; ainp 1/2 => eyeX/Y, ainp 3/4/5 => photodiode1/2/3
% Contrast Ring Protocol, 4x4 C (S), R contrasts - 0,25,50,100%; 
% C rad = 2deg, R rad = 4deg => Annulus width = 2 deg
gridLayout{47} = 1; % 64 electrodes

subjectNames{48} = 'GM'; expDates{48} = '130415'; protocolNames{48} = 'CRS_005'; stimType{48} = 4;
% 64 channels; ainp 1/2 => eyeX/Y, ainp 3/4/5 => photodiode1/2/3
% Contrast Ring Protocol, 3x3 C (S), R contrasts - 0,50,100%; 
% C rad = 2deg, R rad = 2.5,3deg => Annulus width = 0.5, 1 deg
gridLayout{48} = 1; % 64 electrodes

subjectNames{49} = 'GM'; expDates{49} = '130415'; protocolNames{49} = 'CRS_006'; stimType{49} = 4;
% 64 channels; ainp 1/2 => eyeX/Y, ainp 3/4/5 => photodiode1/2/3
% Dual Orientation Protocol, C Ori = 135deg; R(S) Ori = 105,120,135,150,165
% deg; C rad = 2deg, R rad = 45deg (fullscreen)
gridLayout{49} = 1; % 64 electrodes

subjectNames{50} = 'GM'; expDates{50} = '130415'; protocolNames{50} = 'CRS_007'; stimType{50} = 4;
% 64 channels; ainp 1/2 => eyeX/Y, ainp 3/4/5 => photodiode1/2/3
% Dual Orientation Protocol, C Ori = 0,135deg; R(S) Ori = 0,45,90,135 deg
% 135 deg => preferred, 0 deg => non-preferred
% C rad = 2deg, R rad = 45deg (fullscreen)
gridLayout{50} = 1; % 64 electrodes

subjectNames{51} = 'GM'; expDates{51} = '130415'; protocolNames{51} = 'CRS_008'; stimType{51} = 4;
% 64 channels; ainp 1/2 => eyeX/Y, ainp 3/4/5 => photodiode1/2/3
% Dual Phase Protocol, C phase = 0deg; R(S) phase = 0,90,180,270 deg
% C rad = 2deg, R rad = 45deg (fullscreen)
gridLayout{51} = 1; % 64 electrodes

subjectNames{52} = 'GM'; expDates{52} = '130415'; protocolNames{52} = 'CRS_009'; stimType{52} = 4;
% 64 channels; ainp 1/2 => eyeX/Y, ainp 3/4/5 => photodiode1/2/3
% Phase Ring Protocol, C (S) phase = 0deg; R phase = 0,90,180,270 deg
% contrast = 75%,100%
% C rad = 2deg, R rad = 4deg => Annulus width = 2 deg
gridLayout{52} = 1; % 64 electrodes

subjectNames{53} = 'GM'; expDates{53} = '130415'; protocolNames{53} = 'CRS_010'; stimType{53} = 4;
% 64 channels; ainp 1/2 => eyeX/Y, ainp 3/4/5 => photodiode1/2/3
% Orientation Ring Protocol, C ori = 135deg, R ori = 105,120,135,150,165deg
% C rad = 2deg, R rad = 4deg => Annulus width = 2 deg
gridLayout{53} = 1; % 64 electrodes



% 140415
subjectNames{54} = 'GM'; expDates{54} = '140415'; protocolNames{54} = 'CRS_001'; stimType{54} = 4;
% 64 channels; ainp 1/2 => eyeX/Y, ainp 3/4/5 => photodiode1/2/3
gridLayout{54} = 1; % 64 electrodes

subjectNames{55} = 'GM'; expDates{55} = '140415'; protocolNames{55} = 'CRS_002'; stimType{55} = 4;
% 64 channels; ainp 1/2 => eyeX/Y, ainp 3/4/5 => photodiode1/2/3
gridLayout{55} = 1; % 64 electrodes

subjectNames{56} = 'GM'; expDates{56} = '140415'; protocolNames{56} = 'CRS_003'; stimType{56} = 4;
% 64 channels; ainp 1/2 => eyeX/Y, ainp 3/4/5 => photodiode1/2/3
gridLayout{56} = 1; % 64 electrodes

subjectNames{57} = 'GM'; expDates{57} = '140415'; protocolNames{57} = 'CRS_004'; stimType{57} = 4;
% 64 channels; ainp 1/2 => eyeX/Y, ainp 3/4/5 => photodiode1/2/3
gridLayout{57} = 1; % 64 electrodes

subjectNames{58} = 'GM'; expDates{58} = '140415'; protocolNames{58} = 'CRS_005'; stimType{58} = 5;
% 64 channels; ainp 1/2 => eyeX/Y, ainp 3/4/5 => photodiode1/2/3
gridLayout{58} = 1; % 64 electrodes

subjectNames{59} = 'GM'; expDates{59} = '140415'; protocolNames{59} = 'CRS_006'; stimType{59} = 4;
% 64 channels; ainp 1/2 => eyeX/Y, ainp 3/4/5 => photodiode1/2/3
gridLayout{59} = 1; % 64 electrodes



% 150415
subjectNames{60} = 'SB'; expDates{60} = '150415'; protocolNames{60} = 'CRS_001'; stimType{60} = 4;
% 64 channels; ainp 1/2 => eyeX/Y, ainp 3/4/5 => photodiode1/2/3
gridLayout{60} = 1; % 64 electrodes

subjectNames{61} = 'SB'; expDates{61} = '150415'; protocolNames{61} = 'CRS_002'; stimType{61} = 4;
% 64 channels; ainp 1/2 => eyeX/Y, ainp 3/4/5 => photodiode1/2/3
gridLayout{61} = 1; % 64 electrodes

subjectNames{62} = 'SB'; expDates{62} = '150415'; protocolNames{62} = 'CRS_003'; stimType{62} = 4;
% 64 channels; ainp 1/2 => eyeX/Y, ainp 3/4/5 => photodiode1/2/3
gridLayout{62} = 1; % 64 electrodes

subjectNames{63} = 'SB'; expDates{63} = '150415'; protocolNames{63} = 'CRS_004'; stimType{63} = 4;
% 64 channels; ainp 1/2 => eyeX/Y, ainp 3/4/5 => photodiode1/2/3
gridLayout{63} = 1; % 64 electrodes

subjectNames{64} = 'SB'; expDates{64} = '150415'; protocolNames{64} = 'CRS_005'; stimType{64} = 4;
% 64 channels; ainp 1/2 => eyeX/Y, ainp 3/4/5 => photodiode1/2/3
gridLayout{64} = 1; % 64 electrodes

subjectNames{65} = 'SB'; expDates{65} = '150415'; protocolNames{65} = 'CRS_006'; stimType{65} = 4;
% 64 channels; ainp 1/2 => eyeX/Y, ainp 3/4/5 => photodiode1/2/3
gridLayout{65} = 1; % 64 electrodes

subjectNames{66} = 'SB'; expDates{66} = '150415'; protocolNames{66} = 'CRS_007'; stimType{66} = 4;
% 64 channels; ainp 1/2 => eyeX/Y, ainp 3/4/5 => photodiode1/2/3
gridLayout{66} = 1; % 64 electrodes

subjectNames{67} = 'SB'; expDates{67} = '150415'; protocolNames{67} = 'CRS_008'; stimType{67} = 5;
% 64 channels; ainp 1/2 => eyeX/Y, ainp 3/4/5 => photodiode1/2/3
gridLayout{67} = 1; % 64 electrodes


% 070515
subjectNames{68} = 'alpa'; expDates{68} = '070515'; protocolNames{68} = 'GRF_001'; stimType{68} = 4;
% 27 channels; ainp 1/2/6 => eyeX/Y/PD, ainp 3/4/5 => photodiode1/2/3
% SF, Ori Tuning; 6 SFs (0.5 to 16), 4 orientations (0,45,90,135)
gridLayout{68} = 1; % should be type 2

subjectNames{69} = 'alpa'; expDates{69} = '070515'; protocolNames{69} = 'GRF_002'; stimType{69} = 4;
% 27 channels; ainp 1/2/6 => eyeX/Y/PD, ainp 3/4/5 => photodiode1/2/3
% SF, Ori Tuning; 6 SFs (0.5 to 16), 4 orientations (0,45,90,135)
gridLayout{69} = 1; % should be type 2

subjectNames{70} = 'alpa'; expDates{70} = '070515'; protocolNames{70} = 'GRF_003'; stimType{70} = 4;
% 27 channels; ainp 1/2/6 => eyeX/Y/PD, ainp 3/4/5 => photodiode1/2/3
% SF, Ori Tuning; 6 SFs (0.5 to 16), 4 orientations (0,45,90,135)
gridLayout{70} = 1; % should be type 2


for n = 1:length(listProtocolName)
    disp([num2str(n) ' : ' listProtocolName{n}]);
end


