%% Generalized Bond Models (not specific proteins)
kT = 4.114; % [pN*nm] at T = 298 K
koff0 = 1; % Intrinsic Off-Rate
Fb = 5; % Coefficent [pN] in Bell Model Slip Bond: koff = koff0*exp(F/Fb)
xb = kT/Fb; % Coefficent [nm] in Bell Model Slip Bond: koff = koff0*exp(F*xb/kT)

IdealBond_TwoPathBellModelParams.k01 = koff0;
IdealBond_TwoPathBellModelParams.xb1 = 0;
IdealBond_TwoPathBellModelParams.k02 = 0;
IdealBond_TwoPathBellModelParams.xb2 = 0;

SlipBond_TwoPathBellModelParams.k01 = koff0;
SlipBond_TwoPathBellModelParams.xb1 = xb;
SlipBond_TwoPathBellModelParams.k02 = 0;
SlipBond_TwoPathBellModelParams.xb2 = 0;

CatchSlipBond_TwoPathBellModelParams.k01 = koff0*.9;
CatchSlipBond_TwoPathBellModelParams.xb1 = -1*2*xb;
CatchSlipBond_TwoPathBellModelParams.k02 = koff0*.1;
CatchSlipBond_TwoPathBellModelParams.xb2 = xb;
