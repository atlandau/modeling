function dv = odeCellFunction(t,v,odeModel)
% persistent stepTime iteration
% fprintf(1,'%.2f\n',1000*t);

cdv = odeModel.allocateCurrent;

% Add Leak Current
cdv = cdv + (-v + odeModel.physPrm(:,2))./odeModel.physPrm(:,1);

% Add Parent Current
pidx = odeModel.parentIdx;
cdv(pidx) = cdv(pidx) + (v(odeModel.parent(pidx))-v(pidx))./odeModel.physPrm(pidx,4);

% Add daughter current for each daughter
for nd = 1:size(odeModel.daughtersIdx,2)
    didx = odeModel.daughtersIdx(:,nd);
    cdv(didx) = cdv(didx) - (v(didx)-v(odeModel.daughters{nd}))./odeModel.physPrm(odeModel.daughters{nd},4);
end

% Add synaptic current for each synapse
for ns = 1:size(odeModel.synapseIdx,2)
    sidx = odeModel.synapseIdx(:,ns);
    cdv(sidx) = cdv(sidx) + odeModel.synConductance{ns}(t).*(odeModel.synReversal{ns}-v(sidx));
end

% Add voltage-clamp current
for nvc = 1:size(odeModel.vcIdx,2)
    vidx = odeModel.vcIdx(:,nvc);
    cdv(vidx) = cdv(vidx) + (odeModel.voltageClamp{nvc}(t) - v(vidx))./odeModel.access{nvc};
end

% Add current-clamp current
for ncc = 1:size(odeModel.ccIdx,2)
    cidx = odeModel.ccIdx(:,ncc);
    cdv(cidx) = cdv(cidx) + odeModel.currentClamp{ncc}(t);
end

% Convert to dv
dv = cdv./odeModel.physPrm(:,3);

    

    
    