function dydt = camODE_wFluor(t,y,iCalciumAmp,systemPrms,camPrms,fluorPrms)
% dydt = calmodulinModel_1(t,y,c,p)
% 
% t is time (which is irrelevant because this process is memoryless)
% y is the data 
%   y(1): [Ca]free -- concentration free calcium in cell
%   y(2:10): [CaM] in each group [00;10;20;01;11;21;02;12;22]
%   y(11): [CaFB] -- concentration calcium bound to fluorescent buffer
%
% systemPrms provides parameters required for computation
%   systemPrms.v is volume of compartment
%   systemPrms.rest is rest calcium concentration
%   systemPrms.beta is the rate (1/ms) of the yoked extrusion mechanisms
%
% camPrms provide cam parameters - [CaM] determined by iState
%   camPrms.konTC
%   camPrms.konRC
%   camPrms.koffTC
%   camPrms.koffRC
%   camPrms.konTN
%   camPrms.konRN
%   camPrms.koffTN
%   camPrms.koffRN
%
% fluorPrms provide fluorescent buffer parameters
%   fluorPrms.kon is association rate of buffer
%   fluorPrms.koff is dissociation rate of buffer
%   fluorPrms.bconc is the total concentration of the buffer
% 

%y(1) = cafree
% -------
%y(2) = N0_C0
%y(3) = N1_C0
%y(4) = N2_C0
%y(5) = N0_C1
%y(6) = N1_C1
%y(7) = N2_C1
%y(8) = N0_C2
%y(9) = N1_C2
%y(10) = N2_C2
% -------
%y(11) = fluorBound2Calcium

% Calcium current stuff
restCurrent = systemPrms.beta * systemPrms.rest;
curr = restCurrent + ica(t,iCalciumAmp)/(2*96485*systemPrms.v);
extrusion = systemPrms.beta * y(1);

% On/Off reactions with calmodulin
C0_on = 2*camPrms.konTC * y(1) * y(2:4); %[00,10,20]
C1_on = camPrms.konRC * y(1) * y(5:7); %[01,11,21]
C1_off = camPrms.koffTC * y(5:7); %[01,11,21]
C2_off = 2*camPrms.koffRC * y(8:10); %[02,12,22]

N0_on = 2*camPrms.konTN * y(1) * y([2 5 8]); %[00,01,02]
N1_on = camPrms.konRN * y(1) * y([3 6 9]); %[10,11,12]
N1_off = camPrms.koffTN * y([3 6 9]); %[10,11,12]
N2_off = 2*camPrms.koffRN * y([4 7 10]); %[20,21,22]

% Transition Matrix
cLobeTransitions = [-C0_on + C1_off, C0_on - C1_off - C1_on + C2_off, C1_on - C2_off];
nLobeTransitions = [-N0_on + N1_off, N0_on - N1_off - N1_on + N2_off, N1_on - N2_off]';
transitionMatrix = cLobeTransitions + nLobeTransitions;

% Fluorescent Buffer Reaction
fluorAssociation = y(1)*(fluorPrms.bt - y(11))*fluorPrms.kon;
fluorDissociation = y(11)*fluorPrms.koff;
fluorExchange = fluorAssociation - fluorDissociation;

% Calcium Movement
calciumToCaM = sum(C0_on + C1_on + N0_on + N1_on - C1_off - C2_off - N1_off - N2_off);
dCalcium = curr - extrusion - calciumToCaM - fluorExchange;

%Output
dydt = [dCalcium; transitionMatrix(:); fluorExchange];











