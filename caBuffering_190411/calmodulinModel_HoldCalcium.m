function dydt = calmodulinModel_HoldCalcium(t,y,p)
% dydt = calmodulinModel_1(t,y,p)
% 
% model where dydt for free calcium is set to 0 so [calcium] is clamped
%
% t is time (which is irrelevant because this process is memoryless)
% y is the data - nomenclature (N/C - lobe) (0/1/2 - # Ca ions bound)
%               - # Ca ions considered deterministic with T/R state of CaM
% 
%     y(1) = cafree
%     -------
%     y(2) = N0_C0
%     y(3) = N1_C0
%     y(4) = N2_C0
%     y(5) = N0_C1
%     y(6) = N1_C1
%     y(7) = N2_C1
%     y(8) = N0_C2
%     y(9) = N1_C2
%     y(10) = N2_C2
%
% p provides parameters required for computation
%   - The total concentration of CaM nested into initial state 
%   - Volume and other terms relating to calcium influx and efflux are
%     irrelevant because we're clamping calcium
%
%   p.konTC
%   p.konRC
%   p.koffTC
%   p.koffRC
%   p.konTN
%   p.konRN
%   p.koffTN
%   p.koffRN
%   

y = y(:);

% On/Off reactions with calmodulin
C0_on = 2*p.konTC * y(1) * y(2:4); %[00,10,20]
C1_on = p.konRC * y(1) * y(5:7); %[01,11,21]
C1_off = p.koffTC * y(5:7); %[01,11,21]
C2_off = 2*p.koffRC * y(8:10); %[02,12,22]

N0_on = 2*p.konTN * y(1) * y([2 5 8]); %[00,01,02]
N1_on = p.konRN * y(1) * y([3 6 9]); %[10,11,12]
N1_off = p.koffTN * y([3 6 9]); %[10,11,12]
N2_off = 2*p.koffRN * y([4 7 10]); %[20,21,22]

% Transition Matrix
cLobeTransitions = [-C0_on + C1_off, C0_on - C1_off - C1_on + C2_off, C1_on - C2_off];
nLobeTransitions = [-N0_on + N1_off, N0_on - N1_off - N1_on + N2_off, N1_on - N2_off]';
transitionMatrix = cLobeTransitions + nLobeTransitions;

% Ca Reactions
dydt = [0; transitionMatrix(:)];











