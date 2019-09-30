function dv = vcSpineDendrite(t,v,prm)
% spine has:
%   1) a depolarizing synaptic current
%   2) a hyperpolarizing leak to dendrite
%   3) a hyperpolarizing leak from leak channels
% dendrite has:
%   1) a depolarizing current from the spine
%   2) a hyperpolarizing leak channel current

% synaptic current is conductance times driving force, w/ Erev=0
iSynapse = -v(1)*prm.gsyn(t);

cdvs = iSynapse - (v(1)-v(2))/prm.rneck - (v(1)-prm.e)/prm.rhead;
cdvd = (v(1)-v(2))/prm.rneck - (v(2)-prm.e)/prm.rdend;
dv = [cdvs/prm.chead; cdvd/prm.cdend];


