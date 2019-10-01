function dv = odeCellFunction(t,v,odeCellModel)

NC = length(v);

% leak current
iLeak = (-v + odeCellModel.physPrm(:,2))./odeCellModel.physPrm(:,1);

% current from parent
idxParent = ~isnan(odeCellModel.parent);
iParent = zeros(NC,1);
iParent(idxParent) = (v(odeCellModel.parent(idxParent))-v(idxParent))./odeCellModel.physPrm(idxParent,4);

% current to daughters
iDaughters = zeros(NC,1);
ND = size(odeCellModel.daughters,2);
for nd = 1:ND % iterate up to maximum number of daughters for any compartment
    idxDaughter = ~isnan(odeCellModel.daughters(:,nd)); % idx of compartments with nd daughters
    voltageDifference = v(idxDaughter)-v(odeCellModel.daughters(idxDaughter,nd)); % voltage difference of those pairs
    axialDaughter = odeCellModel.physPrm(odeCellModel.daughters(idxDaughter,nd),4); % axial of daughter
    iDaughters(idxDaughter) = iDaughters(idxDaughter) - voltageDifference./axialDaughter; % updated daughter current
end

% Synaptic Current
conductance = @(t,gpeak,tpeak,tstart) ((t-tstart)>0).* gpeak .*exp(1)/tpeak.*(t-tstart).*exp(-(t-tstart)./(tpeak-tstart));
iSynaptic = zeros(NC,1);
if any(cellfun(@(c) ~isempty(c),odeCellModel.synapse(:,1), 'uni', 1))
    NS = size(odeCellModel.synapse,2);
    for ns = 1:NS
        idxSynapse = cellfun(@(c) ~isempty(c), odeCellModel.synapse(:,ns), 'uni', 1); % idx of comp. with ns synapses
        synapseReversal = cellfun(@(c) c(1,1), odeCellModel.synapse(idxSynapse,ns), 'uni', 1);
        drivingForce = synapseReversal - v(idxSynapse);
        synapseConductance = cellfun(@(c) conductance(t,c(1,2),c(1,3),c(1,4)), odeCellModel.synapse(idxSynapse,ns), 'uni', 1);
        iSynaptic(idxSynapse) = iSynaptic(idxSynapse) + synapseConductance.*drivingForce;
    end
end

% Current from Voltage Clamp
iVoltageClamp = zeros(NC,1);
NVC = size(odeCellModel.voltageClamp,2);
for nvc = 1:NVC
    idxVC = cellfun(@(c) ~isempty(c), odeCellModel.voltageClamp(:,nvc), 'uni', 1); 
    vHold = cellfun(@(c) c(t), odeCellModel.voltageClamp(idxVC,nvc), 'uni', 1);
    vDiff = vHold - v(idxVC);
    iVoltageClamp(idxVC) = iVoltageClamp(idxVC) + vDiff./odeCellModel.vcAccess(idxVC);
end

% Current from current-clamp
iCurrentClamp = zeros(NC,1); 
NCC = size(odeCellModel.currentClamp,2);
for ncc = 1:NCC
    idxCC = cellfun(@(c) ~isempty(c), odeCellModel.currentClamp(:,ncc), 'uni', 1); 
    iCurrentClamp(idxCC) = iCurrentClamp(idxCC) + cellfun(@(c) c(t), odeCellModel.currentClamp(idxCC,ncc), 'uni', 1);
end

cdv = iLeak + iParent + iDaughters + iSynaptic + iVoltageClamp + iCurrentClamp;
dv = cdv./odeCellModel.physPrm(:,3);


% sv = v(402)*1000;
% ddv = v(201)*1000;
% sdC = 1e12 * 1e-3*(sv-ddv)/odeCellModel.physPrm(402,4);
% tC = cdv(201)*1e12;
% dC = iDaughters(201)*1e12;
% pC = iParent(201)*1e12;
% fprintf(1, ' S:%.1fmV, D:%.1fmV; C:%.1fpA, TC:%.1fpA, CPar:%.1fpA, CDau:%.1fpA',sv,ddv,sdC,tC,pC,dC);
% fprintf(1,'\n');

    
    