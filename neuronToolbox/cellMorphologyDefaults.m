function defaultCell = cellMorphologyDefaults()
% generate default values for a cortical L2/3 cell 
% based off of a few rough observations from 2P images
% Andrew Landau - September 2019

% Soma
defaultCell.soma.morphid = 0;
defaultCell.soma.type = 0;
defaultCell.soma.diameter = 25;

% Apical Main
defaultCell.aMain.morphid = 1;
defaultCell.aMain.type = 1;
defaultCell.aMain.diameter = 1.8;
defaultCell.aMain.length = 100;
defaultCell.aMain.link = 0;
defaultCell.aMain.location = 0;

% Apical Branches
defaultCell.aBranch1.morphid = 2;
defaultCell.aBranch1.type = 1;
defaultCell.aBranch1.diameter = 1;
defaultCell.aBranch1.length = 90;
defaultCell.aBranch1.link = 1;
defaultCell.aBranch1.location = 57;

defaultCell.aBranch2.morphid = 3;
defaultCell.aBranch2.type = 1;
defaultCell.aBranch2.diameter = 1;
defaultCell.aBranch2.length = 85;
defaultCell.aBranch2.link = 1;
defaultCell.aBranch2.location = 57;

defaultCell.aBranch3.morphid = 4;
defaultCell.aBranch3.type = 1;
defaultCell.aBranch3.diameter = 1;
defaultCell.aBranch3.length = 195;
defaultCell.aBranch3.link = 1;
defaultCell.aBranch3.location = 64;

defaultCell.aBranch4.morphid = 5;
defaultCell.aBranch4.type = 1;
defaultCell.aBranch4.diameter = 1.4;
defaultCell.aBranch4.length = 17;
defaultCell.aBranch4.link = 1;
defaultCell.aBranch4.location = 100;

defaultCell.aBranch5.morphid = 6;
defaultCell.aBranch5.type = 1;
defaultCell.aBranch5.diameter = 1;
defaultCell.aBranch5.length = 17;
defaultCell.aBranch5.link = 5;
defaultCell.aBranch5.location = 17;

defaultCell.aBranch6.morphid = 7;
defaultCell.aBranch6.type = 1;
defaultCell.aBranch6.diameter = 0.8;
defaultCell.aBranch6.length = 54;
defaultCell.aBranch6.link = 1;
defaultCell.aBranch6.location = 86;

defaultCell.aBranch7.morphid = 8;
defaultCell.aBranch7.type = 1;
defaultCell.aBranch7.diameter = 1;
defaultCell.aBranch7.length = 30;
defaultCell.aBranch7.link = 6;
defaultCell.aBranch7.location = 17;

defaultCell.aBranch8.morphid = 9;
defaultCell.aBranch8.type = 1;
defaultCell.aBranch8.diameter = 1;
defaultCell.aBranch8.length = 30;
defaultCell.aBranch8.link = 6;
defaultCell.aBranch8.location = 17;

% Spines
defaultCell.spine1.morphid = 10;
defaultCell.spine1.type = 2;
defaultCell.spine1.diameter = 0.8;
defaultCell.spine1.length = 250e6;
defaultCell.spine1.link = 1;
defaultCell.spine1.location = 50;

defaultCell.spine2.morphid = 11;
defaultCell.spine2.type = 2;
defaultCell.spine2.diameter = 0.8;
defaultCell.spine2.length = 250e6;
defaultCell.spine2.link = 9;
defaultCell.spine2.location = 15;