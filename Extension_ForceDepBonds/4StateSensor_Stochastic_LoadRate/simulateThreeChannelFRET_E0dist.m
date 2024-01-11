function PixelData = simulateThreeChannelFRET_E0dist(state_D1A1, state_D1A0, state_D0A1, state_D0A0, E0, Toggle_Noise, varargin)

% Computes three channel FRET measurements of a mixed population of n MTS
% each with a specified FRET Eff for the ideal state (E0, which should be
% based on some Eff vs Force LUT), and existing in one of four possible
% states based on the health of the donor and acceptor FP (D1A1, D1A0,D0A1,
% or D0A0).

% Assumptions: 
% (1) No intermolecular FRET
% (2) Acceptors in non-functional state (A0) cannot quench donor.

% Inputs
% E0: nx1 vector of E0 value for each MTS (FRET Eff for ideal state)
% state_D1A1: nx1 vector indicating if MTS is in D1A1 state
% state_D1A0: nx1 vector indicating if MTS is in D1A0 state
% state_D0A1: nx1 vector indicating if MTS is in D0A1 state
% state_D0A0: nx1 vector indicating if MTS is in D0A0 state
% Toggle_Noise: use to turn on/off imaging noise (1=On, 0=Off)

%% Ground truth calibration constants
G_true = 1.65;
k_true = 1.00;
gammaM_true = G_true;
betaX_true = 1/G_true/k_true;

%% Donor Bleedthru and Acceptor Direct Excitation constants
alphaBT_true = @(x) 0.75; % Actual donor bleedthru coefficient. Can be constant or function of intensity, x=Idd_TOT.
deltaDE_true = @(x) 0.25; % Actual acceptor direct excitation coefficient. Can be constant or function of intensity, x=Idd_TOT.
alphaBT_est = 0.75; % Estimated donor bleedthru coefficient using donor only experiments and intensity-independent assumption
deltaDE_est = 0.25; % Estimated acceptor direct excitation coefficient using acceptor only experiments and intensity-independent assumption

%% True Signal in Each Channel
% Check that each MTS has unique state defined
check = state_D1A1+state_D1A0+state_D0A1+state_D0A0;
if sum(check)~=length(check)
    error('Each MTS was not assigned a unique state.');
end

% Compute Total Channel Intensities for the general mixed population of
% sensor states based on numbers of sensors in each state, ideal FRET Eff
% of each sensor, and and the ground truth calibration constants.
k = 100; % mean signal in AA channel per acceptor fluorophore
Iaa_True = sum( (state_D1A1+state_D0A1) .* k );
Idd_True = sum( (state_D1A1.*(1-E0) + state_D1A0) .* k/gammaM_true/betaX_true );
Ida_True = sum( E0.*state_D1A1 .*k/betaX_true ) + alphaBT_true(Idd_True).*Idd_True + deltaDE_true(Iaa_True).*Iaa_True;

% Old Version Assuming Single Homogenous E0
% Iaa_True = (num_D0A1+num_D1A1)*k;
% Idd_True = ((1-E0)*num_D1A1 + num_D1A0)*k/gammaM_true/betaX_true;
% Ida_True = E0*num_D1A1*k/betaX_true + alphaBT_true(Idd_True).*Idd_True + deltaDE_true(Iaa_True).*Iaa_True;

%% Actual Signal with Additive Noise Models
% Actual Signal = True Signal + Read Noise + Photon Noise
% Read Noise: (Signal-Indep) Gaussian Read Noise ~ N(mu=0,sigma=sigma_G) [Mannam et al. Optical 2022]
% Photon Noise: (Signal-Dep) Poisson Photon Noise ~ N(mu=0,sigma=sigma_P*sqrt(T)) [Mannam et al. Optical 2022]
% Combine noise models assuming they are independent [Mannam et al. Optical 2022]
% Equivalent to Noise_MPG ~ N(0,sigma_G+a*sigma_P*sqrt(T)) [Mannam et al. Optical 2022]
if isempty(varargin)
    sigma_independent = 400;
    sigma_dependent = 0; 
else
    sigma_independent = varargin{1};
    sigma_dependent = varargin{2}; 
end
numPixels = 1;
Iaa = Iaa_True + Toggle_Noise * ( sigma_independent*randn(numPixels,1) + sigma_dependent.*sqrt(Iaa_True).*randn(numPixels,1) );
Idd = Idd_True + Toggle_Noise * ( sigma_independent*randn(numPixels,1) + sigma_dependent.*sqrt(Idd_True).*randn(numPixels,1) );
Ida = Ida_True + Toggle_Noise * ( sigma_independent*randn(numPixels,1) + sigma_dependent.*sqrt(Ida_True).*randn(numPixels,1) );

% Set negative intensities to zero
Iaa(Iaa<0)=0;
Idd(Idd<0)=0;
Ida(Ida<0)=0;

%% Corrected FRET Intensity
% Compute corrected FRET intensity (Fc or Ida_corr) using estimated
% coefficients for donor bleedthru and acceptor direct excitation
Fc = Ida - alphaBT_est*Idd - deltaDE_est*Iaa;

%% Apparent FRET Eff
% Compute Apparent FRET Eff (Eapp)
Eapp = Fc./(gammaM_true*Idd + Fc);

%% Apparent Stoichiometry
% Comptue Apparent Stoichiometry (Sapp)
Sapp = (gammaM_true*Idd + Fc)./(gammaM_true*Idd + Fc + Iaa/betaX_true);

%% Write Out Data
num_D1A1 = sum( state_D1A1 );
num_D1A0 = sum( state_D1A0 );
num_D0A1 = sum( state_D0A1 );
num_D0A0 = sum( state_D0A0 );

ratio_FreeAcceptorToIdeal = num_D0A1./num_D1A1;
ratio_FreeDonorToIdeal = num_D1A0./num_D1A1;
PixelData = table(num_D1A1, num_D1A0, num_D0A1, num_D0A0,...
    ratio_FreeAcceptorToIdeal, ratio_FreeDonorToIdeal,...
    Iaa, Idd, Ida, Fc, Eapp, Sapp); 

end