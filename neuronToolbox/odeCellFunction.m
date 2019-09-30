function dv = odeCellFunction(~,v,odeCellModel)

NC = length(v);

% current from this component
iLeak = (-v + odeCellModel.physPrm(:,2))./odeCellModel.physPrm(:,1);

% current from parent
idxParent = odeCellModel.parents~=0;
iParent = zeros(NC,1);
iParent(idxParent) = (v(odeCellModel.parents(idxParent))-v(idxParent))./odeCellModel.physPrm(idxParent,4);

% current from daughters
ND = cellfun(@(c) sum(c>0), odeCellModel.daughters, 'uni', 1);
idxOneDaughter = ND>0;
idxFirstDaughter = cellfun(@(c) c(1), odeCellModel.daughters(idxOneDaughter), 'uni', 1);
iDaughters = zeros(NC,1);
iDaughters(idxOneDaughter) = -(v(idxOneDaughter)-v(idxFirstDaughter))./odeCellModel.physPrm(idxFirstDaughter,4);
for nd = find(ND(:)'>1)
    for d = 2:ND(nd)
        idxDaughter = odeCellModel.daughters{nd}(d);
        if idxDaughter==0, continue, end
        iDaughters(nd) = iDaughters(nd) - (v(nd) - v(idxDaughter))/odeCellModel.physPrm(idxDaughter,4);
    end
end

cdv = iParent + iLeak + iDaughters;
dv = cdv./odeCellModel.physPrm(:,3);


% Insert to add current in particular compartment
%idxHold = ceil(NC/2);
%cdv(idxHold) = cdv(idxHold) + 100e-12;

    
    