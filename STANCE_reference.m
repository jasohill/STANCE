%Header
% PURPOSE: Provides optional global level reference structures or loading of 
% spreadsheet data for various resources of interest to the user, 
% ce.g. description of various atlas/parcellation ROIs.
% 
% Authors: Jason E. Hill & Ryan Hayden
% Updated: 26 SEPT 2016

% from aal_allvalues.txt in the STANCE/MNI folder
% RAVI MED STUDENT, Brain textbooks library

atlas.name = 'aal';
atlas.fullername = 'Automatic Anatomical Labeled atlas';
%% SECTION TITLE
label{1} = 'Precentral_L';      
label{2} = 'Precentral_R';      
label{3} = 'Frontal_Sup_L';     
label{4} = 'Frontal_Sup_R';
label{5} = 'Frontal_Sup_Orb_L';
label{6} = 'Frontal_Sup_Orb_R';
label{7} = 'Frontal_Mid_L';
label{8} = 'Frontal_Mid_R';
label{9} = 'Frontal_Mid_Orb_L';   
label{10} = 'Frontal_Mid_Orb_R';
label{11} = 'Frontal_Inf_Oper_L';
label{12} = 'Frontal_Inf_Oper_R';
label{13} = 'Frontal_Inf_Tri_L';  %%%   ?????
label{14} = 'Frontal_Inf_Tri_R';  %%%   ?????
label{15} = 'Frontal_Inf_Orb_L';  %%%   ???
label{16} = 'Frontal_Inf_Orb_R';  %%%   ???
label{17} = 'Rolandic_Oper_L';    %%%   ??
label{18} = 'Rolandic_Oper_R';    %%%   ??
label{19} = 'Supp_Motor_Area_L';  %%%   ???
label{20} = 'Supp_Motor_Area_R';  %%%   ???
label{21} = 'Olfactory_L';        %%%   ? sulculus?
label{22} = 'Olfactory_R';        %%%   ? sulculus?
label{23} = 'Frontal_Sup_Medial_L';
label{24} = 'Frontal_Sup_Medial_R'; %%%    ???? Orbital or orbitofrontal
label{25} = 'Frontal_Med_Orb_L';    %%%    ???? Same
label{26} = 'Frontal_Med_Orb_R';
label{27} = 'Rectus_L';
label{28} = 'Rectus_R';
label{29} = 'Insula_L';
label{30} = 'Insula_R';
label{31} = 'Cingulum_Ant_L';
label{32} = 'Cingulum_Ant_R';
label{33} = 'Cingulum_Mid_L';
label{34} = 'Cingulum_Mid_R';
label{35} = 'Cingulum_Post_L';
label{36} = 'Cingulum_Post_R';
label{37} = 'Hippocampus_L';
label{38} = 'Hippocampus_R';
label{39} = 'ParaHippocampal_L';
label{40} = 'ParaHippocampal_R';
label{41} = 'Amygdala_L';
label{42} = 'Amygdala_R';
label{43} = 'Calcarine_L';
label{44} = 'Calcarine_R';
label{45} = 'Cuneus_L';
label{46} = 'Cuneus_R';
label{47} = 'Lingual_L';
label{48} = 'Lingual_R';
label{49} = 'Occipital_Sup_L';
label{50} = 'Occipital_Sup_R';
label{51} = 'Occipital_Mid_L';
label{52} = 'Occipital_Mid_R';
label{53} = 'Occipital_Inf_L';
label{54} = 'Occipital_Inf_R';
label{55} = 'Fusiform_L';
label{56} = 'Fusiform_R';
label{57} = 'Postcentral_L';
label{58} = 'Postcentral_R';
label{59} = 'Parietal_Sup_L';
label{60} = 'Parietal_Sup_R';
label{61} = 'Parietal_Inf_L';
label{62} = 'Parietal_Inf_R';
label{63} = 'SupraMarginal_L'; %BA(40)
label{64} = 'SupraMarginal_R';
label{65} = 'Angular_L';
label{66} = 'Angular_L';
label{67} = 'Precuneus_L';
label{68} = 'Precuneus_R';
label{69} = 'Paracentral_Lobule_L';
label{70} = 'Paracentral_Lobule_R';
label{71} = 'Caudate_L';
label{72} = 'Caudate_R';
label{73} = 'Putamen_L';
label{74} = 'Putamen_R';
label{75} = 'Pallidum_L';
%label{76} = '';
%label{77} = '';
%label{78} = '';
%label{79} = '';
%label{80} = '';
%label{81} = '';
%label{82} = '';
%label{83} = '';
%label{84} = '';
%label{85} = '';
%label{} = '';
%label{} = '';
%label{} = '';
%label{} = '';
%label{} = '';

atlas.label = label;
%% ABBREVIATION
% DESCRIPTIVE TEXT %CHECK THE TEXT FILE
abbreviation{1} = 'FAG';
abbreviation{2} = 'FAD';
abbreviation{3} = 'F1G';
abbreviation{3} = 'F1D';
%abbreviation{4} = 'F
%abbreviation{5} = 'F
%abbreviation{6} = 'F
%abbreviation{7} = 'F
%abbreviation{8} = 'F
%abbreviation{9} = 'F
%abbreviation{10} = 'F
%abbreviation{11} = 'F
%abbreviation{12} = 'F
%abbreviation{13} = 'F
%abbreviation{14} = 'F
%abbreviation{15} = 'F
%abbreviation{16} = 'F
%abbreviation{17} = 'F
%abbreviation{18} = 'F
%abbreviation{19} = 'F
%abbreviation{20} = 'F
%abbreviation{21} = 'F
%abbreviation{22} = 'F
%abbreviation{23} = 'F
%abbreviation{24} = 'F
%abbreviation{25} = 'F
%abbreviation{26} = 'F
%abbreviation{27} = 'F
%abbreviation{28} = 'F
%abbreviation{29} = 'F
%abbreviation{30} = 'F
%abbreviation{31} = 'F
%abbreviation{32} = 'F
%abbreviation{33} = 'F
%abbreviation{34} = 'F
%abbreviation{35} = 'F
%abbreviation{36} = 'F
%abbreviation{37} = 'F
%abbreviation{38} = 'F
%abbreviation{39} = 'F
%abbreviation{40} = 'F
%abbreviation{41} = 'F
%abbreviation{41} = 'F
%abbreviation{42} = 'F
%abbreviation{43} = 'F
%abbreviation{44} = 'F
%abbreviation{45} = 'F
%abbreviation{46} = 'F
%abbreviation{47} = 'F
%abbreviation{48} = 'F
%abbreviation{49} = 'F
%abbreviation{50} = 'F
%abbreviation{51} = 'F
%abbreviation{52} = 'F
%abbreviation{53} = 'F
%abbreviation{54} = 'F
%abbreviation{55} = 'F
%abbreviation{56} = 'F
%abbreviation{57} = 'F
%abbreviation{58} = 'F
%abbreviation{59} = 'F
%abbreviation{60} = 'F
%abbreviation{61} = 'F
%abbreviation{62} = 'F
%abbreviation{63} = 'F
%abbreviation{64} = 'F
%abbreviation{65} = 'F
%abbreviation{66} = 'F
%abbreviation{67} = 'F
%abbreviation{68} = 'F
%abbreviation{69} = 'F
%abbreviation{70} = 'F
%abbreviation{71} = 'F
%abbreviation{72} = 'F
%abbreviation{73} = 'F
%abbreviation{74} = 'F
%abbreviation{75} = 'F
%...
atlas.abbreviation = abbreviation;
%% CODE
% DESCRIPTIVE TEXT
code(1) = 2001;
code(2) = 2002;
code(3) = 2101;
code(4) = 2102;
code(5) = 2111;
code(6) = 2112;
code(7) = 2201;
code(8) = 2202;
code(9) = 2211;
code(10) = 2212;
code(11) = 2301;
code(12) = 2302;
code(13) = 2311;
code(14) = 2312;
code(15) = 2321;
code(16) = 2322;
code(17) = 2331;
code(18) = 2332;
code(19) = 2401;
code(20) = 2402;
code(21) = 2501;
code(22) = 2502;
code(23) = 2601;
code(24) = 2602;
code(25) = 2611;
code(26) = 2612;
code(27) = 2701;
code(28) = 2702;
code(29) = 3001;
code(30) = 3002;
code(31) = 4001;
code(32) = 4002;
code(33) = 4012;
code(34) = 4021;
code(35) = 4022;
code(36) = 4101;
code(37) = 4101;
code(38) = 4102;
code(39) = 4111;
code(40) = 4112;
code(41) = 4201;
code(42) = 4202;
code(43) = 5001;
code(44) = 5002;
code(45) = 5011;
code(46) = 5012;
code(47) = 5021;
code(48) = 5022;
code(49) = 5101;
code(50) = 5102;
code(51) = 5202;
code(52) = 5202;
code(53) = 5301;
code(54) = 5302;
code(55) = 5401;
code(56) = 5402;
code(57) = 6001;
code(58) = 6002;
code(59) = 6101;
code(60) = 6102;
code(61) = 6201;
code(62) = 6202;
code(63) = 6211;
code(64) = 6212;
code(65) = 6221;
code(66) = 6222;
code(67) = 6301;
code(68) = 6302;
code(69) = 6401;
code(70) = 6402;
code(71) = 7001;
code(72) = 7002;
code(73) = 7011;
code(74) = 7012;
code(75) = 7021;

atlas.code = code;
%% FULL ANATOMICAL NAME
% full anatomical name
% Orbitofrontal = orbital division of frontal cortex
full_name{1} = 'left precentral gyrus';
full_name{2} = 'right precentral gyrus';
full_name{3} = 'left superior frontal gyrus on the lateral surface';
full_name{4} = 'right superior frontal gyrus on the lateral surface';
full_name{5} = 'left superior gyrus on the orbital surface of the frontal lobe';
full_name{6} = 'right superior gyrus on the orbital surface of the frontal lobe';
full_name{7} = 'left middle gyrus on the frontal lobe';
full_name{8} = 'right middle gyrus on the frontal lobe';
full_name{9} = 'left middle gyrus on the frontal lobe of the orbital division';
full_name{10} = 'right middle gyrus on the orbital surface of the frontal lobe';
full_name{11} = 'left frontal lobe on the inferior operculum';
full_name{12} = 'right frontal lobe on the inferior operculum';
full_name{13} = 'left triangular division of the frontal lobe of the inferior gyrus';
full_name{14} = 'right triangular division of the frontal lobe of the inferior gyrus';
full_name{15} = 'left orbital division of the inferior gyrus on the frontal lobe';
full_name{16} = 'right orbital division of the inferior gyrus on the frontal lobe';
full_name{17} = 'left rolandic operculum';
full_name{18} = 'right frontal inferior';
full_name{19} = 'left supplementary motor cortex';
full_name{20} = 'right supplemntary motor cortex';
full_name{21} = 'left olfactory sulculus';
full_name{22} = 'right olfactory sulculus';
full_name{23} = 'left superior frontal gyrus on the medial surface';
full_name{24} = 'right superior frontal gyrus on the medial surface';
full_name{25} = 'left medial orbitofrontal cortex';
full_name{26} = 'right frontal medial orbital sulculus';
full_name{27} = 'left gyrus rectus';
full_name{28} = 'right gyrus rectus';
full_name{29} = 'left insular cortex';
full_name{30} = 'right insular cortex';
full_name{31} = 'left anterior cingulate cortex'; %Found -ulate cortex?
full_name{32} = 'right anterior cingulate cortex';
full_name{33} = 'left middle cingulate cortex';
full_name{34} = 'right middle cingulate cortex';
full_name{35} = 'left posterior cingulate cortex';
full_name{36} = 'right posterior cingulate cortex';
full_name{37} = 'left hippocampus';
full_name{38} = 'right hippocampus';
full_name{39} = 'left parahippocampal gyrus';
full_name{40} = 'right parahippocampal gyrus';
full_name{41} = 'left amygdala';                  %(uncus)
full_name{42} = 'right amygdala';
full_name{43} = 'left calcarine fissure';
full_name{44} = 'right calcarine fissure';
full_name{45} = 'left cuneus';
full_name{46} = 'right cuneus';
full_name{47} = 'left lingual gyrus';
full_name{48} = 'right lingual gyrus';
full_name{49} = 'left superior occipital gyrus';
full_name{50} = 'right superior occipital gyrus';
full_name{51} = 'left middle occipital gyrus';
full_name{52} = 'right middle occipital gyrus';
full_name{53} = 'left inferior occipital gyrus';
full_name{54} = 'right inferior occipital gyrus';
full_name{55} = 'left fusiform gyrus';
full_name{56} = 'right fusiform gyrus';
full_name{57} = 'left postcentral gyrus';
full_name{58} = 'right postcentral gyrus';
full_name{59} = 'left superior parietal cortex';
full_name{60} = 'right superior parietal cortex';
full_name{61} = 'left inferior parietal cortex';
full_name{62} = 'right inferior parietal cortex';
full_name{63} = 'left supramarginal gyrus';
full_name{64} = 'right supramarginal gyrus';
full_name{65} = 'left angular gyrus';
full_name{66} = 'right angular gyrus';
full_name{67} = 'left precuneus';
full_name{68} = 'right precuneus';
full_name{69} = 'left paracentral lobule';
full_name{70} = 'right paracentral lobule';
full_name{71} = 'left caudate nucleus';
full_name{72} = 'right caudate nucleus';
full_name{73} = 'left putamen';
full_name{74} = 'right putamen';
full_name{75} = 'left pallidum';

atlas.fullname = full_name;
%% latin anatomical name
latin{1} = 'sinister gyrus precentralis';
latin{2} = 'dexter gyrus precentralis';
latin{3} = 'sinister gyrus frontalis superior';
latin{4} = 'dexter gyrus frontalis superior';
latin{5} = 'sinister gyrus frontalis superior';
latin{6} = 'dexter gyrus frontalis superior';
latin{7} = 'sinister gyrus frontalis medius';
latin{8} = 'dexter gyrus frontalis medius';
latin{9} = 'sinister medial frontal orbital gyrus';
latin{10} = 'dexter medial frontal orbital gyrus';
latin{11} = 'sinister frontal inferior operculum';
latin{12} = 'dexter frontal inferior operculum';
latin{13} = 'sinister frontal inferior';
latin{14} = 'dexter frontal inferior';
latin{15} = 'sinister frontal inferior orbital';
latin{16} = 'dexter frontal inferior orbital';
latin{17} = 'sinister rolandic operculum';
latin{18} = 'dexter rolandic operculum';
latin{19} = 'sinister frontal inferior';
latin{20} = 'dexter frontal inferior';
latin{21} = 'sinister olfactory sulculus';
latin{22} = 'dexter olfactory sulculus';
latin{23} = 'sinister ...'; 
latin{24} = 'dexter ...';
latin{25} = 'sinister frontal medial orbital sulculus';
latin{26} = 'dexter frontal medial orbital sulculus';
latin{27} = 'sinister gyrus rectus';
latin{28} = 'dexter gyrus rectus';
latin{29} = 'sinister insular cortex';
latin{30} = 'dexter insular cortex';
latin{31} = 'sinister anterior ...'; 
latin{32} = 'dexter anterior ...';
latin{33} = 'sinister middle ...';
latin{34} = 'dexter middle ...';
latin{35} = 'sinister posterior ...'; 
latin{36} = 'dexter posterior ...';
latin{37} = 'sinister hippocampus';
latin{38} = 'dexter hippocampus';
latin{39} = 'sinister ...';
latin{40} = 'dexter ...';
latin{41} = 'sinister corpus amygdaloideum';
latin{42} = 'dexter corpus amygdaloideum';
latin{43} = 'sinister calcarine sulcus';
latin{44} = 'dexter calcarine sulcus';
latin{45} = 'sinister cuneus';
latin{46} = 'dexter cuneus';
latin{47} = 'sinister ...';
latin{48} = 'dexter ...';
latin{49} = 'sinister ...';
latin{50} = 'dexter ...';
%latin{51} = 'sinister
%latin{52} = 'dexter
%latin{53} = 'sinister
%latin{54} = 'dexter
%latin{55} = 'sinister
%latin{56} = 'dexter
%latin{57} = 'sinister
%latin{58} = 'dexter
%latin{59} = 'sinister
%latin{60} = 'dexter
%latin{61} = 'sinister
%latin{62} = 'dexter
%latin{63} = 'sinister
%latin{64} = 'dexter
%latin{65} = 'sinister
%latin{66} = 'dexter
%latin{67} = 'sinister
%latin{68} = 'dexter
%latin{69} = 'sinister
%latin{70} = 'dexter
%latin{71} = 'sinister
%latin{72} = 'dexter
%latin{73} = 'sinister
%latin{74} = 'dexter
%latin{75} = 'sinister globus pallidus';
atlas.latin = latin;
%% also known as
alias{1} = 'motor strip';
alias{2} = 'motor strip';
%alias{3} = [];
%alias{4} = [];
%alias{5} = [];
%alias{6} = [];
%alias{7} = [];
%alias{8} = [];
%alias{9} = [];
%alias{10} = [];
%alias{11} = [];
%alias{12} = [];
%alias{13} = [];
%alias{14} = [];
%alias{15} = [];
%alias{16} = [];
%alias{17} = [];
%alias{18} = [];
%alias{19} = [];
%alias{20} = [];
%alias{21} = [];
%alias{22} = [];
%alias{23} = [];
%alias{24} = [];
%alias{25} = [];
%alias{26} = [];
%alias{27} = [];
%alias{28} = [];
%alias{29} = [];
%alias{55} = 'discontinuous occippitotemporal gyrus';
%alias{67} = 'quadrate lobule of Foville';
%alias{73} = 'corpus striatum';
atlas.alias = alias;
%% Brodmann area correspondance
BA{1} = 4;
BA{2} = 4;
BA{3} = [4 6 8];
BA{4} = [4 6 8];
% BA{5} = [10 11];
% BA{6} = [10];
% BA{7} = [45 48 44 9];
% BA{8} = [10];
% BA{9} = [46 45 10];
% BA{10} = [10 32 45 46 47];
% BA{11} = [38 46 47];
% BA{12} = [11 45 46 47];
% BA{13} = [45];
% BA{14} = [47];
% BA{15} = [11 45];
% BA{16} = [11 47];
% BA{17} = [6 22 37 42 48];
% BA{18} = [6 43 48];
% BA{19} = [6 8 9 32];
% BA{20} = [6 32];
% BA{21} = [25];
% BA{22} = [25];
% BA{23} = [6 9 10 32];
% BA{24} = [10 32];
% BA{25} = [11 12 13];
% BA{26} = [10 11 32];
% BA{27} = [10 11];
% BA{28} = [10 11];
% BA{29} = [48];
% BA{30} = [48];
% BA{31} = [10 24 25];
% BA{32} = [11 24 25 32];
% BA{33} = [8 24 32];
% BA{34} = [23 24 32];
% BA{35} = [23 26 29];
% BA{36} = [23 26 29 30];
% BA{37} = [27 28 30 34 36];
% BA{38} = [20 34 35 37 48];
% BA{39} = [30 37 48];
% BA{40} = [30 37 48];
% BA{41} = [34 28];
% BA{42} = [28 34 35 36];
% BA{43} = [17 19];
% BA{44} = [18 19];
% BA{45} = [17 18 23 26 30];
% BA{46} = [17 18 19];
% BA{47} = [18 19 37];
% BA{48} = [19];
% BA{49} = [17 18 19 27 30];
% BA{50} = [18 19];
% BA{51} = [];
% BA{52} = [];
% BA{53} = [];
% BA{54} = [];
% BA{55} = [];
% BA{56} = [];
% BA{57} = [];
% BA{58} = [];
% BA{59} = [];
% BA{60} = [];
% BA{61} = [];
% BA{62} = [];
% BA{63} = [];
% BA{64} = [];
% BA{65} = [];
% BA{66} = [];
% BA{67} = [];
% BA{68} = [];
% BA{69} = [];
% BA{70} = [];
% BA{71} = [];
% BA{72} = [];
% BA{73} = [];
% BA{74} = [];
% BA{75} = [];
% BA{76} = [];
% BA{77} = [];
% BA{78} = [];
% BA{79} = [];
% BA{} = [];
atlas.BA = BA;
%% Function Descriptions
% details about function (see e.g. http://www.wikipedia.org or http://medical-dictionary.thefreedictionary.com/
function_description{1} = 'Primary motor cortex. Lesions of the precentral gyrus result in paralysis of the contralateral side of the body.';
function_description{2} = 'Primary motor cortex. Lesions of the precentral gyrus result in paralysis of the contralateral side of the body.';
function_description{3} = 'Can induce laughter when stimulated with electric current. fMRI experiments have found evidence that the superior frontal gyrus is involved in self-awareness, in coordination with the action of the sensory system.';
function_description{4} = 'fMRI experiments have found evidence that the superior frontal gyrus is involved in self-awareness, in coordination with the action of the sensory system.';
%...
atlas.function_description = function_description;
%% References
%references{3,1} = % wikipedia
%references{3,2} = % fMRI: Goldberg et al.
%references{3,3} = % laughter: ...
%references{3,4} = https://radiopaedia.org/articles/ % Radiology
%references{3,5} = http://umich.edu/~cogneuro/jpg/Brodmann.html% broadmanA
atlas.references = references;
%% details about location
location_description{1} = 'The precentral gyrus lies in front of the postcentral gyrus - mostly on the lateral (convex) side of the cerebral hemispheres - from which it is separated by the central sulcus. Its anterior border is represented by the precentral sulcus, while inferiorly it borders to the lateral fissure (Sylvian fissure). Medially, it is contiguous with the paracentral lobule.';
location_description{2} = 'The precentral gyrus lies in front of the postcentral gyrus - mostly on the lateral (convex) side of the cerebral hemispheres - from which it is separated by the central sulcus. Its anterior border is represented by the precentral sulcus, while inferiorly it borders to the lateral fissure (Sylvian fissure). Medially, it is contiguous with the paracentral lobule.';
location_description{3} = 'The superior frontal gyrus (SFG) makes up about one third of the frontal lobe of the human brain. It is bounded laterally by the superior frontal sulcus. It is more of a region than a true gyrus.';
location_description{4} = 'The superior frontal gyrus (SFG) makes up about one third of the frontal lobe of the human brain. It is bounded laterally by the superior frontal sulcus. It is more of a region than a true gyrus.';
%location_description{5} = '';
%location_description{6} = '';
%location_description{7} = '';
%location_description{8} = '';
%location_description{9} = '';
%location_description{10} = '';
%location_description{11} = '';
%location_description{12} = '';
%location_description{13} = '';
%location_description{14} = '';
%location_description{15} = '';
%location_description{16} = '';
%location_description{17} = '';
%location_description{18} = '';
%location_description{19} = '';
%location_description{20} = '';
%location_description{21} = '';
%location_description{22} = '';
%location_description{23} = '';
%location_description{24} = '';
%location_description{25} = '';
%location_description{26} = '';
%location_description{27} = '';
%location_description{28} = '';
%location_description{29} = '';
%location_description{30} = '';
%location_description{31} = '';
%location_description{32} = '';
%location_description{33} = '';
%location_description{34} = '';
%location_description{35} = '';
%location_description{36} = '';
%location_description{37} = '';
%location_description{38} = '';
%location_description{39} = '';
%location_description{} = '';
%location_description{} = '';
%location_description{} = '';
%location_description{} = '';
%location_description{} = '';
%location_description{} = '';
%location_description{} = '';
%location_description{} = '';
%location_description{} = '';
%location_description{} = '';
%location_description{} = '';
%location_description{} = '';
%location_description{} = '';
%location_description{} = '';
%location_description{} = '';
%location_description{} = '';
%location_description{} = '';
%... 
atlas.location_description = location_description;
%% Color
color(1,:) = [203 203 203];
color(2,:) = [203 203 203];
color(3,:) = [202 148 148];
color(4,:) = [202 148 148];
%color(5,:) = [
%color(6,:) = [
%color(7,:) = [
%color(8,:) = [
%color(9,:) = [
%color(10,:) = [
%color(11,:) = [
%color(,:) = [

reference.atlas{1} = atlas;
%% Save 
save('reference.mat','reference');