function [odeCellModel,cellMorph,cellTable,cStructure] = generateCell(cellMorph,membranePhys,dx)
%
% -------------------------------------------------------------------------
% cellMorph provides instructions for constructing the cell morphology
%
%
% -------------------------------------------------------------------------
% membranePhys is basic physiological parameters
%   *TYPE*
%    - can either produce a passive membrane using *PASSIVE DEFAULTS* or
%      produce an active membrane using *ACTIVE DEFAULTS*
%    - user can overwrite values by adding them as fields to physPrms
%    - note that field names must match default field names
%   *PASSIVE DEFAULTS*
%    - .axial is axial resistivity (defaults to 100 Ohm*cm)
%    - .memRes is membrane resistivity (defaults to 20000 Ohm-cm^2)
%    - .memCap is specific capacitance (defaults to 1µF/cm^2)
%   *ACTIVE DEFAULTS*
%    - can add specific channels 
%
% -------------------------------------------------------------------------
% dx is maximum resolution - default is 10µm 
% -------------------------------------------------------------------------
% Andrew Landau September 2019

% handle inputs
if nargin<1 || isempty(cellMorph), cellMorph = cellMorphologyDefaults(); end
if nargin<2 || isempty(membranePhys), membranePhys = membraneParameterDefaults(0); end
if nargin<3 || isempty(dx), dx = 10; end


% construct table of cell structure parameters
cellTable = constructCellTable(cellMorph);

% check if cell structure is valid
[valid,cellStructureValid] = checkCellStructure(cellTable); % check validity
if valid~=1, disp(cellStructureValid); return, end % print message to screen if invalid

% build compartment structure
cStructure = buildCompartments(cellTable,membranePhys,dx); % build each compartment


% build ode model of cell
odeCellModel.physPrm = cell2mat(cStructure(2,:)');
odeCellModel.parents = cell2mat(cStructure(4,:)');

daughters = cellfun(@transpose, cStructure(5,:), 'uni', 0); % fuck u matlab
odeCellModel.daughters = [daughters{:}]';

end



%% build compartment structure for cell
function cStructure = buildCompartments(cellTable,membranePhys,dx)
    % {1,:} vector of compartment distances
    % {2,:} parameters of each compartment [Rm, Em, Cm, Ra]
    % {3,:} global compartment index
    % {4,:} idx to parent compartments (use parent's axial resistance)
    % {5,:} idx to daughter compartments (use daughter's axial resistance)

    standardPhys = [membranePhys.sMemRes, membranePhys.eMem, membranePhys.sMemCap, membranePhys.axialRes];

    % build each compartment
    NC = size(cellTable,2);
    cStructure = cell(5,NC);
    numCompartments = zeros(NC,1);
    for c = 1:NC
        switch cellTable(2,c)
            case 0 % soma
                numCompartments(c) = 1; % soma is 1 compartment

                somaArea = 4*pi*(cellTable(3,c)/2)^2;
                rSoma = membranePhys.sMemRes / somaArea;
                cSoma = membranePhys.sMemCap * somaArea;

                cStructure{1,c} = 0; % single compartment
                cStructure{2,c} = [rSoma, membranePhys.eMem, cSoma, 0];
                cStructure{3,c} = sum(numCompartments(1:c-1))+1; 

            case 1 % dendrite
                % discretize dendrite with a maximum segment length shorter
                % than dx, that still splits at all discretization points
                links = (cellTable(5,:)==cellTable(1,c)); % idx of compartments that link to this one
                linkloc = cellTable(6,links);
                iDisc = unique([0, cellTable(4,c), linkloc]); % initial discretization - unique and sorted
                iDistance = diff(iDisc);
                numSegments = ceil(iDistance/dx);
                segmentIdx = [0, cumsum(numSegments)]; 
                distance = zeros(sum(numSegments),1);
                for segment = 1:length(numSegments)
                    idx = segmentIdx(segment)+1:segmentIdx(segment+1);
                    distance(idx) = iDistance(segment)/numSegments(segment);
                end

                relativeDistance = distance / dx;
                numMinorCompartments = length(distance);
                numCompartments(c) = numMinorCompartments; % store number of compartments

                circumference = pi*cellTable(3,c);
                unitArea = circumference*dx;
                unitPhys = standardPhys .* [1/unitArea 1 unitArea dx/(pi*(cellTable(3,c)/2)^2)];

                minorCompartmentConversion = [1./relativeDistance, ones(numMinorCompartments,1), relativeDistance, relativeDistance];

                cStructure{1,c} = distance;
                cStructure{2,c} = minorCompartmentConversion .* unitPhys;
                cStructure{3,c} = sum(numCompartments(1:c-1))+(1:numMinorCompartments)';

            case 2 % spine
                % diameter field (cellTable(3,c) is head diameter
                % length field (cellTable(4,c) is neck resistance)

                numCompartments(c) = 1; % spine is 1 compartment

                spineArea = 4*pi*(cellTable(3,c)/2)^2;
                rSpine = membranePhys.sMemRes / spineArea;
                cSpine = membranePhys.sMemCap * spineArea;

                cStructure{1,c} = 0; % single compartment
                cStructure{2,c} = [rSpine, membranePhys.eMem, cSpine, cellTable(4,c)];
                cStructure{3,c} = sum(numCompartments(1:c-1))+1;
            
            otherwise
                error('Does not recognize compartment type, 0:soma, 1:dendrite, 2:spine');
        end
    end
    
    % Link neighboring compartments in tree
    for c = 1:NC
        switch cellTable(2,c)
            case 0 % soma
                idxNotSoma = cellTable(2,:)~=0; % the soma has a zero in the link field so we need to ignore itself
                links = (cellTable(5,:)==cellTable(1,c)); % idx of compartments that link to this one
                linksToSoma = (idxNotSoma & links); % non-somatic compartments linking to soma

                daughterIdx = cellfun(@(idx) idx(1), cStructure(3,linksToSoma),'uni',1);
                cStructure{4,c} = 0; % Soma has no parent compartment
                cStructure{5,c} = {daughterIdx}; % daughters link from first compartment

            case 1 % dendrite
                daughters = (cellTable(5,:)==cellTable(1,c)); % idx of compartments that link to this one
                position = round(10*cumsum(cStructure{1,c}))/10;

                numMinorCompartments = length(cStructure{3,c});
                cStructure{4,c} = zeros(numMinorCompartments,1);
                cStructure{5,c} = cell(numMinorCompartments,1);
                for nmc = 1:numMinorCompartments
                    % Parent Compartment
                    if nmc==1
                        % parent is on another major compartment
                        morphidParent = cellTable(5,c);
                        idxParentCompartment = cellTable(1,:)==morphidParent;
                        if cellTable(2,idxParentCompartment)==0
                            % parent is soma
                            cStructure{4,c}(nmc) = cStructure{3,idxParentCompartment};
                        else
                            % parent is dendrite
                            parentpos = round(10*cumsum(cStructure{1,idxParentCompartment}))/10;
                            parentloc = cellTable(6,c);
                            parentidx = (parentpos == parentloc);
                            if sum(parentidx)==0, error('couldn''t link to parent- discretization wrong'); end
                            cStructure{4,c}(nmc) = cStructure{3,idxParentCompartment}(parentidx);
                        end
                    else
                        % parent is previous minor compartment
                        cStructure{4,c}(nmc) = cStructure{3,c}(nmc-1);
                    end
                    % Daughter Compartments
                    if nmc<numMinorCompartments
                        % First daughter is neighbor on major compartment
                        cStructure{5,c}{nmc} = cStructure{3,c}(nmc+1);
                    end
                    linksToCompartment = cellTable(6,:)==position(nmc) & daughters;
                    if any(linksToCompartment)
                        daughterIdx = cellfun(@(idx) idx(1), cStructure(3,linksToCompartment),'uni',1);
                        cStructure{5,c}{nmc} = [cStructure{5,c}{nmc}; daughterIdx(:)];
                    end
                    if isempty(cStructure{5,c}{nmc}), cStructure{5,c}{nmc} = 0; end
                end

            case 2 % spine
                % Spine has parent and no daughters
                morphidParent = cellTable(5,c);
                idxParentCompartment = cellTable(1,:)==morphidParent;
                if cellTable(2,idxParentCompartment)==0
                    % parent is soma
                    cStructure{4,c} = cStructure{3,idxParentCompartment};
                else
                    % parent is dendrite
                    parentpos = round(10*cumsum(cStructure{1,idxParentCompartment}))/10;
                    parentloc = cellTable(6,c);
                    parentidx = (parentpos == parentloc);
                    if sum(parentidx)==0, error('couldn''t link to parent- discretization wrong'); end
                    cStructure{4,c} = cStructure{3,idxParentCompartment}(parentidx);
                end
                cStructure{5,c} = {0}; % no daughters
        end
    end
end


%% check if cell structure is valid
function [valid,cellStructureValid] = checkCellStructure(cellTable)
    valid = 1;
    cellStructureValid = 'valid';
    % More than 1 soma
    if sum(cellTable(2,:)==0) > 1
        valid = 0;
        cellStructureValid = [cellStructureValid{:}, {'more than one soma defined'}]; 
    end
    % Soma not morphid=0
    if cellTable(1,cellTable(2,:)==0)~=0
        valid = 0;
        cellStructureValid = [cellStructureValid{:}, {'soma needs to be morphid 0'}]; 
    end
    % Soma has parent
    if cellTable(5,cellTable(2,:)==0)~=0
        valid = 0;
        cellStructureValid = [cellStructureValid{:}, {'soma cannot have parent compartment'}];
    end
    % Linkage compartment doesn't exist in for at least one compartment
    idxDanglingLinkage = ~ismember(cellTable(5,:),cellTable(1,:));
    if any(idxDanglingLinkage)
        valid = 0;
        dangling = cellTable(1,idxDanglingLinkage);
        missing = cellTable(5,idxDanglingLinkage);
        msg = ['Compartment (',sprintf('%d ',dangling),') linked to nonexistent compartment (',sprintf('%d ',missing),')'];
        cellStructureValid = [cellStructureValid{:}, {msg}];
    end
    % Linkage too far
    idxLinkDendrite = cellTable(5,:)>0;
    linkTooFar = cellTable(6,idxLinkDendrite) > cellTable(4,cellTable(5,idxLinkDendrite)+1);
    if any(linkTooFar)
        valid = 0;
        msg = ['Compartments (',sprintf('%d ',find(linkTooFar)),') link too far on parent compartment'];
        cellStructureValid = [cellStructureValid{:}, {msg}];
    end
    % Spine is parent
    morphidSpines = cellTable(1,(cellTable(2,:)==2));
    if any(ismember(morphidSpines,cellTable(5,:)))
        valid = 0;
        msg = ['Compartments (',sprintf('%d ',find(ismember(cellTable(5,:),morphidSpines))),') link to spines'];
        cellStructureValid = [cellStructureValid{:}, {msg}];
    end
    % Disconnected 
    NC = size(cellTable,2);
    pairs = [cellTable(1,:)' cellTable(5,:)']+1;
    G = accumarray(cat(1, pairs, fliplr(pairs)),1,[NC NC]) + eye(NC);
    if ~checkc(G)
        valid = 0;
        cellStructureValid = [cellStructureValid{:}, {'Neural Morphology is not connected'}]; 
    end
    cellStructureValid = cellStructureValid'; % Transpose for display 
end

%% build cell table
function cellTable = constructCellTable(cellMorph)



    majorCompartments = fieldnames(cellMorph);
    numMajorCompartments = numel(majorCompartments);
    tableNames = {'morphid','type','diameter','length','link','location'};
    numTableNames = numel(tableNames);
    cellTable = zeros(numTableNames,numMajorCompartments);
    for nmc = 1:numMajorCompartments
        for ntn = 1:numTableNames
            if isfield(cellMorph.(majorCompartments{nmc}),tableNames{ntn})
                cellTable(ntn,nmc) = cellMorph.(majorCompartments{nmc}).(tableNames{ntn});
            end
        end
    end
    cellTable(6,:) = round(cellTable(6,:)); % need this for computer rounding errors
end

%% check if cell morphology linkage graph is connected
function conncheckie=checkc(g)
    % Thank you:
    % https://www.mathworks.com/matlabcentral/fileexchange/35347-connectivity-check-for-undirected-graphs
    n=length(g); I=zeros(1,n); I(1)=1; count=1;
    while sum(I)<n && count<n
        for i=1:n,if I(i)==1,for j=1:n,if g(i,j)~=0
            g(j,:)=g(j,:)+g(i,:); g(j,j)=0; I(j)=1; end, end, end, end
        count=count+1;
    end
    if sum(I)==n, conncheckie=1; else, conncheckie=0; end
end







































