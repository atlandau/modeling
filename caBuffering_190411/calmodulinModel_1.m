function dydt = calmodulinModel_1(t,y,ica,p)
% dydt = calmodulinModel_1(t,y,c,p)
% 
% t is time (which is irrelevant because this process is memoryless)
% y is the data 
%   y(1): [Ca]free -- concentration free calcium in cell
%   y(2): [N
%
% ica an inline function giving calcium current at time t
%
% p provides parameters required for computation
%   p.v is volume of compartment
%   p.kon is association rate of buffer
%   p.koff is dissociation rate of buffer
%   p.bconc is the total concentration of the buffer
%   p.beta is the rate (1/ms) of the yoked extrusion mechanisms


% okay fuck
% ahhhhhhhhh
% shitballs


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

restCurrent = prms(

% On/Off reactions with calmodulin
C0_on = 2*p.konTC * cafree * y(2:4); %[00,10,20]
C1_on = p.konRC * cafree * y(5:7); %[01,11,21]
C1_off = p.koffTC * y(5:7); %[01,11,21]
C2_off = 2*p.koffRC * y(8:10); %[02,12,22]

N0_on = 2*p.konTN * cafree * y([2 5 8]); %[00,01,02]
N1_on = p.konRN * cafree * y([3 6 9]); %[10,11,12]
N1_off = p.koffTN * y([3 6 9]); %[10,11,12]
N2_off = 2*p.koffRN * y([4 7 10]); %[20,21,22]

% Transition Matrix
cLobeTransitions = [-C0_on + C1_off, C0_on - C1_off - C1_on + C2_off, C1_on - C2_off];
nLobeTransitions = [-N0_on + N1_off, N0_on - N1_off - N1_on + N2_off, N1_on - N2_off]';
transitionMatrix = cLobeTransitions + nLobeTransitions;

% Ca Reactions
caAssociation = sum(C0_on + C1_on + N0_on + N1_on);
caDissociation = sum(C1_off + C2_off + N1_off + N2_off);

dydt = [caDissociation - caAssociation; transitionMatrix(:)];











