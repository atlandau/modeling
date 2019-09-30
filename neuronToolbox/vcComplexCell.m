function dv = vcComplexCell(t,v,cellPrms,physPrms)
% physPrms is basic physiological parameters




cdvs = (prm.vc(t) - v(1))/prm.rp - (v(1) - prm.es)/prm.rs - (v(1)-v(2))/prm.ra;
cdvd = (v(1) - v(2))/prm.ra - (v(2) - prm.ed)/prm.rd;
dv = [cdvs/prm.cs; cdvd/prm.cs]; 
