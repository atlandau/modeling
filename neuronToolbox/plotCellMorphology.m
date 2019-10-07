function [h,m] = plotCellMorphology(compStructure,cellStructure)

numFeatures = length(cell2mat(compStructure.gci')) + sum(cellStructure.type>2); % minor compartments + auxiliary features

h = struct(); % preallocate handle structure
m = struct(); % preallocate morphology structure
m(numFeatures).gci = 0;

stillBuilding = true;
linkIdx = -1; % Start with -1, which is usually soma but is always the "seed" compartment
cidx = 1;
while stillBuilding
   idx = find(cellStructure.link==linkIdx); % Find all components that link to linkIdx
   for nc = idx(:)'
       switch cellStructure.type(nc)
           case 0 % is soma
               m(cidx).diameter = cellStructure.diameter(nc);
               m(cidx).patchCoord = [];
           case 1
               
               
           case 2
       end
   end
end


% this sucks balls
% 
% the plan is to construct from seed up, using the link dependency to
% choose which compartments to make in each round of the while loop
%
% the initial structure will be in units of µm, therefore accurate and to
% scale. 
% 
% 




