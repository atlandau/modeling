function [odeCellModel,cellMorph,cellStructure,compStructure] = generateCell(cellMorph,membranePhys,dx)
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


% construct cell structure
cellStructure = constructCellStructure(cellMorph);

% check if cell structure is valid
[valid,cellStructureValid] = checkCellStructure(cellStructure); % check validity
if valid~=1, disp(cellStructureValid); return, end % print message to screen if invalid

% build compartment structure
compStructure = buildCompartments(cellStructure,membranePhys,dx); % build each compartment

% build ode model of cell
odeCellModel.physPrm = cell2mat(compStructure.parameters');
odeCellModel.parent = cell2mat(compStructure.parent');
odeCellModel.daughters = matricizeFromCell(compStructure.daughters,0);
odeCellModel.synapse = matricizeFromCell(compStructure.synapse,1);
odeCellModel.voltageClamp = matricizeFromCell(compStructure.voltageClamp,1);
odeCellModel.vcAccess = cell2mat(compStructure.vcAccess');
odeCellModel.currentClamp = matricizeFromCell(compStructure.currentClamp,1);

end



%% build compartment structure for cell
function compStructure = buildCompartments(cellStructure,membranePhys,dx)
    % compStructure.name-- cell of name of each compartment
    % compStructure.position -- vector of end position for each compartment
    % compStructure.parameters-- parameters each comp. [Rm, Em, Cm, Ra]
    % compStructure.gci-- global compartment index
    % compStructure.parents-- idx to parent compartments
    % compStructure.daughters-- idx to daughter compartments
    % compStructure.synapse-- cell of synaptic features
    % compStructure.voltageClamp-- cell of voltage clamp command
    % compStructure.currentClamp-- cell of current clamp command

    standardPhys = [membranePhys.sMemRes, membranePhys.eMem, membranePhys.sMemCap, membranePhys.axialRes];
    
    % build each compartment
    NC = numel(cellStructure.morphid);
    compStructure.name = cell(1,NC);
    compStructure.position = cell(1,NC);
    compStructure.parameters = cell(1,NC);
    compStructure.gci = cell(1,NC);
    compStructure.parent = cell(1,NC);
    compStructure.daughters = cell(1,NC);
    compStructure.synapse = cell(1,NC);
    compStructure.voltageClamp = cell(1,NC);
    compStructure.vcAccess = cell(1,NC);
    compStructure.currentClamp = cell(1,NC);
    
    numCompartments = zeros(1,NC);
    for nc = 1:NC
        compStructure.name{nc} = cellStructure.name{nc};
        switch cellStructure.type(nc)
            case 0 % soma
                numCompartments(nc) = 1; % soma is 1 compartment

                somaArea = 4*pi*(cellStructure.diameter(nc)/2)^2;
                rSoma = membranePhys.sMemRes / somaArea;
                cSoma = membranePhys.sMemCap * somaArea;
                
                compStructure.position{nc} = 0; % single compartment
                compStructure.parameters{nc} = [rSoma, membranePhys.eMem, cSoma, 0];
                compStructure.gci{nc} = sum(numCompartments(1:nc-1))+1;

            case 1 % dendrite
                % discretize dendrite with a maximum segment length shorter
                % than dx, that still splits at all discretization points
                links = (cellStructure.link==cellStructure.morphid(nc)); % idx of compartments that link to this one
                linkloc = cellStructure.location(links);
                iDisc = unique([0, cellStructure.length(nc), linkloc]); % initial discretization - unique and sorted
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
                numCompartments(nc) = numMinorCompartments; % store number of compartments
                
                radius = cellStructure.diameter(nc)/2;
                circumference = 2*pi*radius;
                unitArea = circumference*dx;
                unitPhys = standardPhys .* [1/unitArea 1 unitArea dx/(pi*radius^2)];

                minorCompartmentConversion = [1./relativeDistance, ones(numMinorCompartments,1), relativeDistance, relativeDistance];
                
                compStructure.position{nc} = round(10*cumsum(distance))/10;
                compStructure.parameters{nc} = minorCompartmentConversion.*unitPhys; % broadcasting
                compStructure.gci{nc} = sum(numCompartments(1:nc-1))+(1:numMinorCompartments)';

            case 2 % spine
                % diameter field (cellTable(3,c) is head diameter
                % length field (cellTable(4,c) is neck resistance)
                numCompartments(nc) = 1; % spine is 1 compartment

                spineArea = 4*pi*(cellStructure.diameter(nc)/2)^2;
                rSpine = membranePhys.sMemRes / spineArea;
                cSpine = membranePhys.sMemCap * spineArea;
                
                compStructure.position{nc} = 0; % single compartment
                compStructure.parameters{nc} = [rSpine, membranePhys.eMem, cSpine, cellStructure.rneck(nc)];
                compStructure.gci{nc} = sum(numCompartments(1:nc-1))+1;
                
            case {3,11,12} % synapse, voltage-clamp, current-clamp
                % Not a true compartment, will remove when vectorized
                numCompartments(nc) = 0; % not a compartment
            
            otherwise
                error('Does not recognize compartment type.\n%s',...
                    '0:soma, 1:dendrite, 2:spine, 3:synapse, 11:vclamp, 12:cclamp');
        end
    end
    
    % Link neighboring compartments in tree
    for nc = 1:NC
        switch cellStructure.type(nc)
            case 0 % soma
                % Compartment links
                compStructure.parent{nc} = nan;
                compStructure.synapse{nc}{1} = [];
                compStructure.voltageClamp{nc}{1} = [];
                compStructure.currentClamp{nc}{1} = [];
                
                links = cellStructure.link==cellStructure.morphid(nc) & (numCompartments>0);
                daughterIdx = cellfun(@(idx) idx(1), compStructure.gci(links), 'uni', 1);
                if ~isempty(daughterIdx)
                    compStructure.daughters{nc} = daughterIdx;
                else
                    compStructure.daughters{nc} = nan;
                end
                
                % Auxiliary (synaptic or pipette) links
                auxLinks = find(cellStructure.link==cellStructure.morphid(nc) & (numCompartments==0));
                for alink = auxLinks(:)'
                    switch cellStructure.type(alink)
                        case 3
                            compStructure.synapse{nc}{1}{end+1} = cellStructure.synapse(alink);
                        case 11
                            compStructure.voltageClamp{nc}{1}{end+1} = cellStructure.voltageClamp(alink);
                            compStructure.vcAccess{nc} = cellStructure.vcAccess(alink);
                        case 12
                            compStructure.currentClamp{nc}{1}{end+1} = cellStructure.currentClamp(alink);
                    end
                end
                if isempty(compStructure.vcAccess{nc}), compStructure.vcAccess{nc} = nan; end
                
            case 1 % dendrite
                % Compartment links
                daughters = cellStructure.link==cellStructure.morphid(nc) & (numCompartments>0);
                auxLinks = cellStructure.link==cellStructure.morphid(nc) & (numCompartments==0);
                
                compStructure.parent{nc} = nan(numCompartments(nc),1);
                compStructure.daughters{nc} = cell(numCompartments(nc),1);
                compStructure.synapse{nc} = cell(numCompartments(nc),1);
                compStructure.voltageClamp{nc} = cell(numCompartments(nc),1);
                compStructure.vcAccess{nc} = nan(numCompartments(nc),1);
                compStructure.currentClamp{nc} = cell(numCompartments(nc),1);
                for nmc = 1:numCompartments(nc)
                    % Parent compartments
                    if nmc==1
                        % parent is on another major compartment
                        morphidParent = cellStructure.link(nc);
                        idxParentCompartment = cellStructure.morphid==morphidParent;
                        if morphidParent==-1 % has no parent
                            compStructure.parent{nc}(nmc) = nan;
                        elseif cellStructure.type(idxParentCompartment)==0
                            % parent is soma
                            compStructure.parent{nc}(nmc) = compStructure.gci{idxParentCompartment};
                        else
                            % parent is dendrite
                            parentpos = compStructure.position{idxParentCompartment};
                            parentloc = cellStructure.location(nc);
                            parentidx = parentpos==parentloc;
                            if sum(parentidx)==0, error('couldn''t link to parent- discretization error'); end
                            compStructure.parent{nc}(nmc) = compStructure.gci{idxParentCompartment}(parentidx);
                        end
                    else
                        % parent is previous minor compartment
                        compStructure.parent{nc}(nmc) = compStructure.gci{nc}(nmc-1);
                    end
                    % Daughter Compartments
                    if nmc<numCompartments(nc)
                        % First daughter is neighbor on minor compartment
                        compStructure.daughters{nc}{nmc} = compStructure.gci{nc}(nmc+1);
                    end
                    linksToCompartment = cellStructure.location==compStructure.position{nc}(nmc) & daughters;
                    if any(linksToCompartment)
                        daughterIdx = cellfun(@(idx) idx(1), compStructure.gci(linksToCompartment), 'uni', 1);
                        compStructure.daughters{nc}{nmc} = [compStructure.daughters{nc}{nmc}, daughterIdx];
                    end
                    if isempty(compStructure.daughters{nc}{nmc}), compStructure.daughters{nc}{nmc}=nan; end
                    
                    alinksToCompartment = find(cellStructure.location==compStructure.position{nc}(nmc) & auxLinks);
                    for alink = alinksToCompartment(:)'
                        switch cellStructure.type(alink)
                            case 3
                                compStructure.synapse{nc}{nmc}{end+1} = cellStructure.synapse(alink);
                            case 11
                                compStructure.voltageClamp{nc}{nmc}{end+1} = cellStructure.voltageClamp(alink);
                                compStructure.vcAccess{nc}{nmc}(end+1) = cellStructure.vcAccess(alink);
                            case 12
                                compStructure.currentClamp{nc}{nmc}{end+1} = cellStructure.currentClamp(alink);
                        end
                    end
                end
                
            case 2 % spine
                compStructure.synapse{nc}{1} = [];
                compStructure.voltageClamp{nc}{1} = [];
                compStructure.currentClamp{nc}{1} = [];
                
                morphidParent = cellStructure.link(nc);
                idxParentCompartment = cellStructure.morphid==morphidParent;
                if cellStructure.type(idxParentCompartment)==0
                    % Parent is soma
                    compStructure.parent{nc} = compStructure.gci{idxParentCompartment};
                else
                    % Parent is dendrite
                    parentpos = compStructure.position{idxParentCompartment};
                    parentloc = cellStructure.location(nc);
                    parentidx = parentpos==parentloc;
                    if sum(parentidx)==0, error('couldn''t link to parent- discretization error'); end
                    compStructure.parent{nc} = compStructure.gci{idxParentCompartment}(parentidx);
                end
                compStructure.daughters{nc} = nan; % no daughters
                
                % Auxiliary Links
                auxLinks = find(cellStructure.link==cellStructure.morphid(nc) & (numCompartments==0));
                for alink = auxLinks(:)'
                    switch cellStructure.type(alink)
                        case 3
                            compStructure.synapse{nc}{1}{end+1} = cellStructure.synapse(alink);
                        case 11
                            compStructure.voltageClamp{nc}{1}{end+1} = cellStructure.voltageClamp(alink);
                            compStructure.vcAccess{nc} = cellStructure.vcAccess(alink);
                        case 12
                            compStructure.currentClamp{nc}{1}{end+1} = cellStructure.currentClamp(alink);
                    end
                end
                if isempty(compStructure.vcAccess{nc}), compStructure.vcAccess{nc} = nan; end
        end
    end
    
    trueCompartment = numCompartments>0;
    FN = fieldnames(compStructure);
    for fn = 1:numel(FN)
        compStructure.(FN{fn}) = compStructure.(FN{fn})(trueCompartment);
    end
end


%% check if cell structure is valid
function [valid,cellStructureValid] = checkCellStructure(cellStructure)
    NC = numel(cellStructure.morphid);
    valid = 1;
    cellStructureValid = {''};
    % Type not recognized
    if any(~ismember(cellStructure.type,[0 1 2 3 11 12]))
        valid = 0;
        cellStructureValid = [cellStructureValid{:}, {'some compartment types not recognized'}];
    end
    % More than 1 soma
    idxSoma = cellStructure.type==0;
    if sum(idxSoma) > 1
        valid = 0;
        cellStructureValid = [cellStructureValid{:}, {'more than one soma defined'}]; 
    end
    % Soma has parent
    if cellStructure.link(idxSoma)~=-1
        valid = 0;
        cellStructureValid = [cellStructureValid{:}, {'soma cannot have parent compartment'}];
    end
    % Errors in linkage, either: 
    % 1) linking to a non-existent compartment
    % 2) linking to non-soma or dendrite but not auxiliary component
    % 3) linking to dendrite but past the end of the dendrite
    hasLink = cellStructure.link ~= -1;
    allLinks = cellStructure.link(hasLink);
    [~,idxToLink] = ismember(allLinks,cellStructure.morphid);
    idxDangling = false(1,NC);
    idxDangling(hasLink) = ~ismember(allLinks,cellStructure.morphid);
    if any(idxDangling)
        valid = 0;
        dangling = cellStructure.morphid(idxDangling);
        missing = cellStructure.link(idxDangling);
        msg = ['Compartment (',sprintf('%d ',dangling),') linked to nonexistent compartment (',sprintf('%d ',missing),')'];
        cellStructureValid = [cellStructureValid{:}, {msg}];
    end
    
    
    
    
    
    
    % Dendrites have to link to soma or other dendrites
    % Spines have to link to soma or dendrites
    % Auxiliary Links can link to anything except each other
    
    
%     if any(cellStructure.type(idxToLink)>1)
%         valid = 0;
%         cellStructureValid = [cellStructureValid{:}, {'links must be made to soma and dendrite only'}];
%     end





    % Link on dendrite is too far
    linkToDendrite = false(1,NC);
    linkToDendrite(hasLink) = cellStructure.type(idxToLink)==1;
    [~,idxToDendrite] = ismember(cellStructure.link(linkToDendrite),cellStructure.morphid);
    dendLength = inf*ones(1,NC);
    dendLength(linkToDendrite) = cellStructure.length(idxToDendrite);
    linkDendLoc = zeros(1,NC);
    linkDendLoc(linkToDendrite) = cellStructure.location(linkToDendrite);
    linkTooFar = linkDendLoc>dendLength;
    if any(linkTooFar)
        valid = 0;
        msg = ['Compartments (',sprintf('%d ',find(linkTooFar)),') link too far on parent compartment'];
        cellStructureValid = [cellStructureValid{:}, {msg}];
    end
    % Disconnected 
    fprintf(2,'NOTE THAT the disconnected check might be buggy\n');
    nonNegativeLink = cellStructure.link;
    nonNegativeLink(nonNegativeLink<0)=cellStructure.morphid(nonNegativeLink<0);
    morphids = cellStructure.morphid;
    [~,idxID] = ismember(morphids,morphids);
    [~,nnLinkIdx] = ismember(nonNegativeLink,morphids);
    pairs = [idxID; idxID(nnLinkIdx)]';
    G = accumarray(cat(1, pairs, fliplr(pairs)),1,[NC NC]) + eye(NC);
    if ~checkc(G)
        valid = 0;
        cellStructureValid = [cellStructureValid{:}, {'Neural Morphology is not connected'}]; 
    end
    if valid==1, cellStructureValid='valid'; end
    cellStructureValid = cellStructureValid'; % Transpose for display 
end


%% build cell structure
function cellStructure = constructCellStructure(cellMorph)
    % First get vectors for morphid and type
    % Then load in each feature into the cell structure based on it's type    
    majorCompartments = fieldnames(cellMorph);
    NMC = numel(majorCompartments);
    
    % Morphid & Type
    default = -1*ones(1,NMC);
    defCell = cell(1,NMC);
    cellStructure.name = defCell;
    cellStructure.morphid = default;
    cellStructure.type = default;
    cellStructure.diameter = default;
    cellStructure.length = default;
    cellStructure.rneck = default;
    cellStructure.link = default;
    cellStructure.location = default;
    cellStructure.synapse = defCell;
    cellStructure.voltageClamp = defCell;
    cellStructure.vcAccess = default;
    cellStructure.currentClamp = defCell;
    
    for nmc = 1:NMC
        cComp = cellMorph.(majorCompartments{nmc});
        cellStructure.name{nmc} = majorCompartments{nmc};
        cellStructure.morphid(nmc) = cComp.morphid;
        cellStructure.type(nmc) = cComp.type;
        if isfield(cComp,'diameter') % Soma, dendrite, and spine
            cellStructure.diameter(nmc) = cComp.diameter;
        end
        if isfield(cComp,'link') % usually soma doesn't have a link because it is not a daughter to anything
            cellStructure.link(nmc) = cComp.link;
            if isfield(cComp,'location')
                % Don't put location if linking to soma or spine
                cellStructure.location(nmc) = round(cComp.location); % need this to avoid computer rounding errors
            end
        end
        switch cComp.type
            case 1 % dendrite
                cellStructure.length(nmc) = cComp.length;
            case 2 % spine
                cellStructure.rneck(nmc) = cComp.rneck;
            case 3 % synapse
                % 3 component vector giving:
                % [reversal potential, peak conductance, peak cond. time]
                cellStructure.synapse{nmc} = cComp.properties; 
                
            case 11 % voltage-clamp
                % inline function describing voltage clamp as vc(t)
                cellStructure.voltageClamp{nmc} = cComp.vcCommand;
                cellStructure.vcAccess(nmc) = cComp.vcAccess;
            case 12 % current-clamp
                % inline function describing current clamp as cc(t)
                cellStructure.currentClamp{nmc} = cComp.ccCommand;
        end
    end
end

%% Matricize cell
function mat = matricizeFromCell(C,keepCell)
    mat = cellfun(@(c) c(:)', C, 'uni', 0);
    mat = [mat{:}]';
    maxMat = max(cellfun(@length, mat, 'uni', 1));
    if keepCell
        out = cell(numel(mat),maxMat);
        for i = 1:numel(mat)
            cNum = length(mat{i});
            for c = 1:cNum
                out{i,c} = mat{i}{c}{1};
            end
            for c = cNum+1:maxMat
                out{i,c} = [];
            end
        end
        mat = out;
    else
        mat = cell2mat(cellfun(@(c) [c, nan(1, maxMat-length(c))], mat, 'uni', 0));
    end
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







































