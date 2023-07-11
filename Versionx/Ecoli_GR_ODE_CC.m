function [dxdt,growthRate,kelong,fracCharge,fracRactive,ppGppSynRate,ppGppDegRate] = ...
    Ecoli_GR_ODE_CC(t,x,nutr,cm,hp,ppGpp)

%   t: time
%   x: variables
%   nutr: concentration of growth-limiting nutrient
%   cm: concentration of chloramphenicol
%   hp: host cell parameters

AminoAcid   = abs(x(1));
Ribosome    = abs(x(2));

%   Charged and Uncharged tRNA
chargedTRNA     = hp.ratioTR*Ribosome*AminoAcid/(hp.thetaAaAc+AminoAcid);
unchargedTRNA   = hp.ratioTR*Ribosome*hp.thetaAaAc/(hp.thetaAaAc+AminoAcid);
if (chargedTRNA == 0)
    fracCharge = 0; %   0/0 will give NaN
else
    fracCharge = chargedTRNA/(chargedTRNA+unchargedTRNA);
end

%   RMF concentration is stimulated by ppGpp
rmf = hp.cRMF*ppGpp^2/(hp.thetaPpGppRMF^2+ppGpp^2);

%   Ribosomes are inactivated by RMF and CM
Ractive = (Ribosome-hp.KdRMF*(1+cm/hp.KdCM)-rmf+...
          sqrt((Ribosome-hp.KdRMF*(1+cm/hp.KdCM)-rmf)^2+...
          4*Ribosome*hp.KdRMF*(1+cm/hp.KdCM)))/2/(1+cm/hp.KdCM);
if (Ractive == 0)
    fracRactive = 0;    %   0/0 will give NaN
else
    fracRactive = Ractive/Ribosome;
end

%   Ribosomes bound by CM and RMF
%   Rib_cm = Ractive*cm/hp.KdCM;
%   Rib_rmf = Ribosome-Ractive-Rib_cm;

%   Peptide elongation rate
kelong  = hp.kelongMax*(chargedTRNA/hp.KdAcTRNA)./(1+chargedTRNA/hp.KdAcTRNA+unchargedTRNA/hp.KdUAcTRNA);

%   Concentration of stalled ribosomes
Rsr     = Ractive.*(unchargedTRNA/hp.KdUAcTRNA)./(1+chargedTRNA/hp.KdAcTRNA+unchargedTRNA/hp.KdUAcTRNA);

%   Total protein synthesis rate (in unit of amino acid/h)
kelongTot   = kelong*Ractive;

%   Nutrient uptake and amino acid synthesis
Eprot   = (hp.beta*hp.phiRMax-hp.massR*Ribosome)/hp.massE;
if (Eprot < 0)
    Eprot = 0;
end
JsAa    = nutr*Eprot*(hp.thetaAaTox/(AminoAcid+hp.thetaAaTox));

%   ppGpp synthesis rate
JsPpGpp = hp.ksPpGpp*hp.cRelA*Rsr;

%   Ribosome synthesis rate
JsR     = kelongTot*hp.phiRMax*(hp.thetaPpGppR./(hp.thetaPpGppR+ppGpp))/hp.massR;

%   Degradation rate of R-proteins
%   only inactive ribosomes are degraded!
JdR     = hp.kdR*(1-fracRactive)*Ribosome;

%   Degradation rate of ppGpp
JdPpGpp = hp.kdPpGpp*ppGpp;

%   Specific growth rate
growthRate = kelongTot/hp.beta-hp.mCost;
if (growthRate < 0)
    growthRate = 0;
end

%   Differential equations
%   they are written in form of (synthesis rate)-(active
%   degradation)-(dilution)
dxdt(1) = JsAa      -kelongTot*hp.fAa   -growthRate*AminoAcid;
dxdt(2) = JsR       -JdR                -growthRate*Ribosome;
%dxdt(3) = JsPpGpp   -JdPpGpp            -growthRate*ppGpp;
dxdt    = dxdt';

%   output
ppGppSynRate = JsPpGpp;
ppGppDegRate = JdPpGpp+growthRate*ppGpp;

end

