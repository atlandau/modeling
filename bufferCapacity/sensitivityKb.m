
hpath = '/Users/LandauLand/Documents/Research/SabatiniLab/modeling/bufferCapacity';

%{
To check the sensitivity on a Kb estimate to Kd and resting calcium
relative to concentration of indicator
%}

getKB = @(B,kd,ca) B./(kd + ca);


KD = 800e-9; % Kd of Fluo-5f in solution
B = 300e-6; % Concentration of Fluo-5f at equilibrium


ca = 10e-9:10e-9:500e-6;
kb = getKB(B,KD,ca);

plot(ca,kb);
set(gca,'xscale','log')

maxKB = getKB(B,KD,0);
minKB = getKB(B,KD,500e-3);
line(xlim,[maxKB maxKB],'color','k');
line(xlim,[minKB minKB],'color','k');
line([KD KD],ylim,'color','k')


























