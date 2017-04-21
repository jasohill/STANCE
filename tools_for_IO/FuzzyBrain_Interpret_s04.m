% Fuzzy Model of Anatomical Normal Brain: subject 4
% Fuzzy Models:
% Background   'subject04_bck_v.mnc'
% CSF          'subject04_csf_v.mnc'
% Grey Matter  'subject04_gm_v.mnc'
% White Matter 'subject04_wm_v.mnc'
% Fat          'subject04_fat_v.mnc'
% Muscle       'subject04_muscles_v.mnc'
% Muscle/Skin  'subject04_muscles_skin_v.mnc'
% Skull        'subject04_skull_v.mnc'
% Connective   'subject04_fat2_v.mnc'
% Dura         'subject04_dura_v.mnc'
% Marrow       'subject04_marrow_v.mnc'

clear, clc, close all;

% 4D phantom construction
FuzzyBrain = zeros(217, 181, 181, 12);
maxValues  = zeros(1,12);
1
[background,bck_scaninfo]        = loadminc('subject04_bck_v.mnc');
maxValues(1) = max(background(:));
FuzzyBrain(:,:,:,1) = (double(background(1:2:433,1:2:361,1:2:361)) + double(background(1:2:433,1:2:361,2:2:362)));
FuzzyBrain(:,:,:,1) = FuzzyBrain(:,:,:,1) + double(background(1:2:433,2:2:362,1:2:361));
FuzzyBrain(:,:,:,1) = FuzzyBrain(:,:,:,1) + double(background(1:2:433,2:2:362,2:2:362));
FuzzyBrain(:,:,:,1) = FuzzyBrain(:,:,:,1) + double(background(2:2:434,1:2:361,1:2:361));
FuzzyBrain(:,:,:,1) = FuzzyBrain(:,:,:,1) + double(background(2:2:434,1:2:361,2:2:362));
FuzzyBrain(:,:,:,1) = FuzzyBrain(:,:,:,1) + double(background(2:2:434,2:2:362,1:2:361));
FuzzyBrain(:,:,:,1) = FuzzyBrain(:,:,:,1) + double(background(2:2:434,2:2:362,2:2:362));
clear background;
2
[csf,csf_scaninfo]               = loadminc('subject04_csf_v.mnc');
maxValues(2) = max(csf(:));
FuzzyBrain(:,:,:,2) = (double(csf(1:2:433,1:2:361,1:2:361)) + double(csf(1:2:433,1:2:361,2:2:362)));
FuzzyBrain(:,:,:,2) = FuzzyBrain(:,:,:,2) + double(csf(1:2:433,2:2:362,1:2:361));
FuzzyBrain(:,:,:,2) = FuzzyBrain(:,:,:,2) + double(csf(1:2:433,2:2:362,2:2:362));
FuzzyBrain(:,:,:,2) = FuzzyBrain(:,:,:,2) + double(csf(2:2:434,1:2:361,1:2:361));
FuzzyBrain(:,:,:,2) = FuzzyBrain(:,:,:,2) + double(csf(2:2:434,1:2:361,2:2:362));
FuzzyBrain(:,:,:,2) = FuzzyBrain(:,:,:,2) + double(csf(2:2:434,2:2:362,1:2:361));
FuzzyBrain(:,:,:,2) = FuzzyBrain(:,:,:,2) + double(csf(2:2:434,2:2:362,2:2:362));
clear csf;
3
[gm,gm_scaninfo]                 = loadminc('subject04_gm_v.mnc');
maxValues(3) = max(gm(:));
FuzzyBrain(:,:,:,3) = (double(gm(1:2:433,1:2:361,1:2:361)) + double(gm(1:2:433,1:2:361,2:2:362)));
FuzzyBrain(:,:,:,3) = FuzzyBrain(:,:,:,3) + double(gm(1:2:433,2:2:362,1:2:361));
FuzzyBrain(:,:,:,3) = FuzzyBrain(:,:,:,3) + double(gm(1:2:433,2:2:362,2:2:362));
FuzzyBrain(:,:,:,3) = FuzzyBrain(:,:,:,3) + double(gm(2:2:434,1:2:361,1:2:361));
FuzzyBrain(:,:,:,3) = FuzzyBrain(:,:,:,3) + double(gm(2:2:434,1:2:361,2:2:362));
FuzzyBrain(:,:,:,3) = FuzzyBrain(:,:,:,3) + double(gm(2:2:434,2:2:362,1:2:361));
FuzzyBrain(:,:,:,3) = FuzzyBrain(:,:,:,3) + double(gm(2:2:434,2:2:362,2:2:362));
clear gm;
4
[wm,wm_scaninfo]                 = loadminc('subject04_wm_v.mnc');
maxValues(4) = max(wm(:));
FuzzyBrain(:,:,:,4) = (double(wm(1:2:433,1:2:361,1:2:361)) + double(wm(1:2:433,1:2:361,2:2:362)));
FuzzyBrain(:,:,:,4) = FuzzyBrain(:,:,:,4) + double(wm(1:2:433,2:2:362,1:2:361));
FuzzyBrain(:,:,:,4) = FuzzyBrain(:,:,:,4) + double(wm(1:2:433,2:2:362,2:2:362));
FuzzyBrain(:,:,:,4) = FuzzyBrain(:,:,:,4) + double(wm(2:2:434,1:2:361,1:2:361));
FuzzyBrain(:,:,:,4) = FuzzyBrain(:,:,:,4) + double(wm(2:2:434,1:2:361,2:2:362));
FuzzyBrain(:,:,:,4) = FuzzyBrain(:,:,:,4) + double(wm(2:2:434,2:2:362,1:2:361));
FuzzyBrain(:,:,:,4) = FuzzyBrain(:,:,:,4) + double(wm(2:2:434,2:2:362,2:2:362));
clear wm;
5
[fat,fat_scaninfo]               = loadminc('subject04_fat_v.mnc');
maxValues(5) = max(fat(:));
FuzzyBrain(:,:,:,5) = (double(fat(1:2:433,1:2:361,1:2:361)) + double(fat(1:2:433,1:2:361,2:2:362)));
FuzzyBrain(:,:,:,5) = FuzzyBrain(:,:,:,5) + double(fat(1:2:433,2:2:362,1:2:361));
FuzzyBrain(:,:,:,5) = FuzzyBrain(:,:,:,5) + double(fat(1:2:433,2:2:362,2:2:362));
FuzzyBrain(:,:,:,5) = FuzzyBrain(:,:,:,5) + double(fat(2:2:434,1:2:361,1:2:361));
FuzzyBrain(:,:,:,5) = FuzzyBrain(:,:,:,5) + double(fat(2:2:434,1:2:361,2:2:362));
FuzzyBrain(:,:,:,5) = FuzzyBrain(:,:,:,5) + double(fat(2:2:434,2:2:362,1:2:361));
FuzzyBrain(:,:,:,5) = FuzzyBrain(:,:,:,5) + double(fat(2:2:434,2:2:362,2:2:362));
clear fat;
6
[muscles,muscles_scaninfo]       = loadminc('subject04_muscles_v.mnc');
maxValues(6) = max(muscles(:));
FuzzyBrain(:,:,:,6) = (double(muscles(1:2:433,1:2:361,1:2:361)) + double(muscles(1:2:433,1:2:361,2:2:362)));
FuzzyBrain(:,:,:,6) = FuzzyBrain(:,:,:,6) + double(muscles(1:2:433,2:2:362,1:2:361));
FuzzyBrain(:,:,:,6) = FuzzyBrain(:,:,:,6) + double(muscles(1:2:433,2:2:362,2:2:362));
FuzzyBrain(:,:,:,6) = FuzzyBrain(:,:,:,6) + double(muscles(2:2:434,1:2:361,1:2:361));
FuzzyBrain(:,:,:,6) = FuzzyBrain(:,:,:,6) + double(muscles(2:2:434,1:2:361,2:2:362));
FuzzyBrain(:,:,:,6) = FuzzyBrain(:,:,:,6) + double(muscles(2:2:434,2:2:362,1:2:361));
FuzzyBrain(:,:,:,6) = FuzzyBrain(:,:,:,6) + double(muscles(2:2:434,2:2:362,2:2:362));
clear muscles;
7
[skin,skin_scaninfo]             = loadminc('subject04_muscles_skin_v.mnc');
maxValues(7) = max(skin(:));
FuzzyBrain(:,:,:,7) = (double(skin(1:2:433,1:2:361,1:2:361)) + double(skin(1:2:433,1:2:361,2:2:362)));
FuzzyBrain(:,:,:,7) = FuzzyBrain(:,:,:,7) + double(skin(1:2:433,2:2:362,1:2:361));
FuzzyBrain(:,:,:,7) = FuzzyBrain(:,:,:,7) + double(skin(1:2:433,2:2:362,2:2:362));
FuzzyBrain(:,:,:,7) = FuzzyBrain(:,:,:,7) + double(skin(2:2:434,1:2:361,1:2:361));
FuzzyBrain(:,:,:,7) = FuzzyBrain(:,:,:,7) + double(skin(2:2:434,1:2:361,2:2:362));
FuzzyBrain(:,:,:,7) = FuzzyBrain(:,:,:,7) + double(skin(2:2:434,2:2:362,1:2:361));
FuzzyBrain(:,:,:,7) = FuzzyBrain(:,:,:,7) + double(skin(2:2:434,2:2:362,2:2:362));
clear skin;
8
[skull,skull_scaninfo]           = loadminc('subject04_skull_v.mnc');
maxValues(8) = max(skull(:));
FuzzyBrain(:,:,:,8) = (double(skull(1:2:433,1:2:361,1:2:361)) + double(skull(1:2:433,1:2:361,2:2:362)));
FuzzyBrain(:,:,:,8) = FuzzyBrain(:,:,:,8) + double(skull(1:2:433,2:2:362,1:2:361));
FuzzyBrain(:,:,:,8) = FuzzyBrain(:,:,:,8) + double(skull(1:2:433,2:2:362,2:2:362));
FuzzyBrain(:,:,:,8) = FuzzyBrain(:,:,:,8) + double(skull(2:2:434,1:2:361,1:2:361));
FuzzyBrain(:,:,:,8) = FuzzyBrain(:,:,:,8) + double(skull(2:2:434,1:2:361,2:2:362));
FuzzyBrain(:,:,:,8) = FuzzyBrain(:,:,:,8) + double(skull(2:2:434,2:2:362,1:2:361));
FuzzyBrain(:,:,:,8) = FuzzyBrain(:,:,:,8) + double(skull(2:2:434,2:2:362,2:2:362));
clear skull;
9
[vessels,vessels_scaninfo]       = loadminc('subject04_vessels.mnc');
maxValues(9) = max(vessels(:));
FuzzyBrain(:,:,:,9) = (double(vessels(1:2:433,1:2:361,1:2:361)) + double(vessels(1:2:433,1:2:361,2:2:362)));
FuzzyBrain(:,:,:,9) = FuzzyBrain(:,:,:,9) + double(vessels(1:2:433,2:2:362,1:2:361));
FuzzyBrain(:,:,:,9) = FuzzyBrain(:,:,:,9) + double(vessels(1:2:433,2:2:362,2:2:362));
FuzzyBrain(:,:,:,9) = FuzzyBrain(:,:,:,9) + double(vessels(2:2:434,1:2:361,1:2:361));
FuzzyBrain(:,:,:,9) = FuzzyBrain(:,:,:,9) + double(vessels(2:2:434,1:2:361,2:2:362));
FuzzyBrain(:,:,:,9) = FuzzyBrain(:,:,:,9) + double(vessels(2:2:434,2:2:362,1:2:361));
FuzzyBrain(:,:,:,9) = FuzzyBrain(:,:,:,9) + double(vessels(2:2:434,2:2:362,2:2:362));
clear vessels;
10
[connective,connective_scaninfo] = loadminc('subject04_fat2_v.mnc');
maxValues(10) = max(connective(:));
FuzzyBrain(:,:,:,10) = (double(connective(1:2:433,1:2:361,1:2:361)) + double(connective(1:2:433,1:2:361,2:2:362)));
FuzzyBrain(:,:,:,10) = FuzzyBrain(:,:,:,10) + double(connective(1:2:433,2:2:362,1:2:361));
FuzzyBrain(:,:,:,10) = FuzzyBrain(:,:,:,10) + double(connective(1:2:433,2:2:362,2:2:362));
FuzzyBrain(:,:,:,10) = FuzzyBrain(:,:,:,10) + double(connective(2:2:434,1:2:361,1:2:361));
FuzzyBrain(:,:,:,10) = FuzzyBrain(:,:,:,10) + double(connective(2:2:434,1:2:361,2:2:362));
FuzzyBrain(:,:,:,10) = FuzzyBrain(:,:,:,10) + double(connective(2:2:434,2:2:362,1:2:361));
FuzzyBrain(:,:,:,10) = FuzzyBrain(:,:,:,10) + double(connective(2:2:434,2:2:362,2:2:362));
clear connective;
11
[dura,dura_scaninfo]             = loadminc('subject04_dura_v.mnc');
maxValues(11) = max(dura(:));
FuzzyBrain(:,:,:,11) = (double(dura(1:2:433,1:2:361,1:2:361)) + double(dura(1:2:433,1:2:361,2:2:362)));
FuzzyBrain(:,:,:,11) = FuzzyBrain(:,:,:,11) + double(dura(1:2:433,2:2:362,1:2:361));
FuzzyBrain(:,:,:,11) = FuzzyBrain(:,:,:,11) + double(dura(1:2:433,2:2:362,2:2:362));
FuzzyBrain(:,:,:,11) = FuzzyBrain(:,:,:,11) + double(dura(2:2:434,1:2:361,1:2:361));
FuzzyBrain(:,:,:,11) = FuzzyBrain(:,:,:,11) + double(dura(2:2:434,1:2:361,2:2:362));
FuzzyBrain(:,:,:,11) = FuzzyBrain(:,:,:,11) + double(dura(2:2:434,2:2:362,1:2:361));
FuzzyBrain(:,:,:,11) = FuzzyBrain(:,:,:,11) + double(dura(2:2:434,2:2:362,2:2:362));
clear dura;
12
[marrow,marrow_scaninfo]         = loadminc('subject04_marrow_v.mnc');
maxValues(12) = max(marrow(:));
FuzzyBrain(:,:,:,12) = (double(marrow(1:2:433,1:2:361,1:2:361)) + double(marrow(1:2:433,1:2:361,2:2:362)));
FuzzyBrain(:,:,:,12) = FuzzyBrain(:,:,:,12) + double(marrow(1:2:433,2:2:362,1:2:361));
FuzzyBrain(:,:,:,12) = FuzzyBrain(:,:,:,12) + double(marrow(1:2:433,2:2:362,2:2:362));
FuzzyBrain(:,:,:,12) = FuzzyBrain(:,:,:,12) + double(marrow(2:2:434,1:2:361,1:2:361));
FuzzyBrain(:,:,:,12) = FuzzyBrain(:,:,:,12) + double(marrow(2:2:434,1:2:361,2:2:362));
FuzzyBrain(:,:,:,12) = FuzzyBrain(:,:,:,12) + double(marrow(2:2:434,2:2:362,1:2:361));
FuzzyBrain(:,:,:,12) = FuzzyBrain(:,:,:,12) + double(marrow(2:2:434,2:2:362,2:2:362));
clear marrow;

maxValues

FuzzyBrain = FuzzyBrain/8.0;
% correction
FuzzyBrain(isnan(FuzzyBrain)) = 0.0;
% validate
sum(sum(sum(sum(FuzzyBrain,4)-1.0)))

% phantom output
save('FuzzyBrain_s04.mat','FuzzyBrain');

%clear all;