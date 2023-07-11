function [dxdt,growthRate,AminoAcidSynRate,AminoAcidDegRate] = ...
    Ecoli_GR_ODE_GA_unsatVer(t,x,nutr,cm,hp,model,eta,gamma,ppGpp,AminoAcid)

%   t: time
%   x: variables
%   nutr: concentration of growth-limiting nutrient
%   cm: concentration of chloramphenicol
%   hp: host cell parameters

Ribosome    = abs(x(1));

%   Charged and Uncharged tRNA
switch model
    case {0,1}
        chargedTRNA     = hp.ratioTR*Ribosome*AminoAcid/(hp.thetaAaAc+AminoAcid);
        unchargedTRNA   = hp.ratioTR*Ribosome*hp.thetaAaAc/(hp.thetaAaAc+AminoAcid);
    otherwise
        chargedTRNA     = hp.ratioTR*Ribosome*eta*AminoAcid/hp.thetaAaAc;
        unchargedTRNA   = hp.ratioTR*Ribosome*hp.thetaAaAc/(hp.thetaAaAc+AminoAcid);
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

%   Peptide elongation rate
switch model
    case {0,2}
        kelong  = hp.kelongMax*(chargedTRNA/hp.KdAcTRNA)./(1+chargedTRNA/hp.KdAcTRNA+unchargedTRNA/hp.KdUAcTRNA);
        Rsr     = Ractive.*(unchargedTRNA/hp.KdUAcTRNA)./(1+chargedTRNA/hp.KdAcTRNA+unchargedTRNA/hp.KdUAcTRNA);
    otherwise
        kelong  = hp.kelongMax*gamma*(chargedTRNA/hp.KdAcTRNA)./(1+unchargedTRNA/hp.KdUAcTRNA);
        Rsr     = Ractive.*(unchargedTRNA/hp.KdUAcTRNA)./(1+chargedTRNA/hp.KdAcTRNA+unchargedTRNA/hp.KdUAcTRNA);
end
%   Total protein synthesis rate (in unit of amino acid/h)
kelongTot   = kelong*Ractive;

%   Nutrient uptake and amino acid synthesis
Eprot   = (hp.beta*hp.phiRMax-hp.massR*Ribosome)/hp.massE;
if (Eprot < 0)
    Eprot = 0;
end
JsAa    = nutr*Eprot*(hp.thetaAaTox/(AminoAcid+hp.thetaAaTox));

%   Ribosome synthesis rate
JsR     = kelongTot*hp.phiRMax*(hp.thetaPpGppR./(hp.thetaPpGppR+ppGpp))/hp.massR;

%   Degradation rate of R-proteins
%   only inactive ribosomes are degraded!
JdR     = hp.kdR*(1-fracRactive)*Ribosome;

%   Specific growth rate
growthRate = kelongTot/hp.beta-hp.mCost;
if (growthRate < 0)
    growthRate = 0;
end

%   Differential equations
%   they are written in form of (synthesis rate)-(active
%   degradation)-(dilution)
dxdt(1) = JsR       -JdR                -growthRate*Ribosome;
dxdt    = dxdt';

%   output
AminoAcidSynRate = JsAa;
AminoAcidDegRate = kelongTot*hp.fAa + growthRate*AminoAcid;

end

