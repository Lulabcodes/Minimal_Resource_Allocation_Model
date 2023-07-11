function hostparams = readParameters()

%   # of amino acids in R-sector and affiliated proteins (aa)
%   biological constant
hostparams.('massR')            = 7336*1.6;

%   # of amino acids in E-sector protein (aa)
%   amino acid length of an averaged protein
hostparams.('massE')            = 325;

%   maximum peptide elongation rate (1/h)
%   close to 20 aa/s, the observed maximum elongation rate
hostparams.('kelongMax')        = 19*3600;

%   total amino acid concentration (uM)
%   Paper: Growth rate of Escherichia coli by A Marr (1991)
hostparams.('beta')             = 3.00E+6;

%   RelA concentration (uM)
%   Paper: Analysis of the relA gene product of Escherichia coli by F
%   Pedersen and N Kjeldgaard (1977)
hostparams.('cRelA')            = 0.21;

%   ppGpp synthesis rate per RelA (1/h)
%   Paper: Characterizaation of the tRNA and ribosome-dependent
%   pppGpp-synthesis by recombinant stringent factor from Escherichia coli
%   by K Jenvert et al. (2005)
hostparams.('ksPpGpp')          = 4461*60;

%   ppGpp decay rate (1/h)
%   Paper: Interaction of alleles of the relA, relC and spoT genes in
%   Escherichia coli: analysis of the interconversion of GTP, ppGpp, and
%   pppGpp by N Fiil et al. (1977)
hostparams.('kdPpGpp')          = 0.03*3600;

%   maximum proteome allocation for R sector
%   Close to observed 0.55 in paper
%   Interdependence of cell growth and gene expression: origins and
%   consequences by M Scott et al. (2010)
hostparams.('phiRMax')          = 0.5;

%   dissociation constant of chloramphenicol (uM)
%   Paper: How partially inhibitory concentrations of chloramphenicol
%   affect the growth of Escherichia coli by R Harvey and A Koch
hostparams.('KdCM')             = 2.04;

%   fraction of one particular amino acid
%   biological constant
hostparams.('fAa')              = 0.05;

%-------------------------------
%   Free parameters
%-------------------------------

%   tRNA:Ribosome concencetration ratio
hostparams.('ratioTR')          = 0.22;

%   dissociation constant of charged tRNA (uM)
hostparams.('KdAcTRNA')         = 3.00;

%   dissociation constant of uncharged tRNA (uM)
hostparams.('KdUAcTRNA')        = 30.00;

%   ppGpp activation threshold for ribosome synthesis (uM)
hostparams.('thetaPpGppR')      = 60;

%   amino acid threshold for tRNA aminoacylation (uM)
hostparams.('thetaAaAc')        = 10;

%   amino acid toxicity (uM)
hostparams.('thetaAaTox')       = 1e+4;

%   total RMF concentration (uM)
hostparams.('cRMF')             = 15;

%   dissociation constant of RMF (uM)
hostparams.('KdRMF')            = 1;

%   ppGpp activation threshold for RMF synthesis (uM)
hostparams.('thetaPpGppRMF')    = 232.65;

%   active degradation rate for inactive ribosomes (1/h)
hostparams.('kdR')              = 0.1;

%   mainteinance cost flux (1/h)
hostparams.('mCost')            = 0.1;

end

