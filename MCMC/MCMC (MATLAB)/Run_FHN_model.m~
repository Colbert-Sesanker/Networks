% Run FHN: The FHN model and equations are from the Calderhead paper:

% Girolami, M. and Calderhead, B. (2011), Riemann manifold Langevin and Hamiltonian Monte Carlo methods. 
% Journal of the Royal Statistical Society: Series B (Statistical Methodology), 73: 123–214. 
% doi: 10.1111/j.1467-9868.2010.00765.x


% Add all dierectories in MCMC methods to path
addpath(genpath('./'))
% Close all Figures
close all;

Model.burnin                = 100;
Model.numPosteriorSamples   = 50;

Model.equationName                   = 'FHN'; % name for saving results
% function handle of model equations and sensitivity equations
Model.equations                      = @FitzHughNagumo;
Model.equationsSens                  = @FitzHughNagumoSens2; 
Model.equationsSens_1                = @FitzHughNagumoSens1;
Model.numSampledParams               = 3;

% indexes of observed species in state vector
Model.observedStates                 = [1 2];
Model.unobservedStates               = [];
% indexes of sampled parameters 
% All other parameters are held fixed
sampledParam_idxs                    = [1 2 3];
Model.numSampledParams               = length(sampledParam_idxs);
Model.sampledParam_idxs              = sampledParam_idxs;
% Noise
Model.addedNoise_SD                  = 0.5;
% The initial step size for parameter updates
Model.initialStepSize                = 0.6;
% The step size is adjusted online until acceptance ratio
% is in the range given by 'stepSizeRange'
Model.stepSizeRange                  = [50 70];
% epsilon is for finite differences
Model.epsilon                        = 5e-1; 
Model.zeroMetricTensorDerivatives    = true;
% If true plots trajectories for all proposed parameters
Model.plotProposedTrajectories       = true;

% For RMHMC
Model.numLeapFrogSteps               = 3;
Model.stepSize_RMHMC                 = 3 / Model.numLeapFrogSteps;
Model.maxFixedPointStepsMomentum     = 5;
Model.maxFixedPointStepsPosition     = 3;

% Choose sensitvity type
Model.isAnalyticSens                 = false;
sensitivityMethods                   = getSensitivityMethods();
% 1 = symbolic, 2= finite diff, 3= automatic diff
Model.sensitivityMethod              = sensitivityMethods{3};


% Model specific priors (function handles))
% The Prior Struct holds all information regarding priors

% Used in calculating prior probabilities of current and proposed parameters
Prior.prior                   = @ModelParameterPrior;
% Used in computing natural gradient of posterior
Prior.prior_derivative        = @ModelParameterLogPriorDerivative;
% Used in computing the metric tensor of posterior
Prior.prior_second_derivative = @FHN_prior_second_deriv;
% Used in computing the derivative of metric tensor for the Laplace
% Beltrami Operator
Prior.prior_third_derivative  = @FHN_PriorThirdDeriv;

Model.Prior                   = Prior;

% Model specific integration settings
numTimePts    = 200;
startTime     = 0;
endTime       = 20;
step          = (endTime - startTime) / numTimePts;
timePoints    = startTime: step :endTime;

% The zeros are the initial values for the sensitivitites
numParams  = length(Model.SampledParameter_idxs);
numSpecies = length(Model.ObservedSpecies);
initialValues   = [-1 1];

initialValuesSens   = [ -1 1                                 ... 
                         zeros(1, numParams*numSpecies)       ...
                         zeros(1, sum(1:numParams)*numSpecies)...
                      ];               

initialValuesSens_1 = [ -1 1                                 ... 
                         zeros(1, numParams*numSpecies)       ...                         
                      ]; 
               
               
% Cell array of all parameter names         
paramMap                   = {'a' 'b' 'c'};
params                     = [0.2 0.2 3];
Model.totalParameters      = params;
Model.numTotalParameters   = length(params);
Model.paramMap             = paramMap;
Model.initialValues        = initialValues;
Model.initialValues_Sens   = initialValuesSens;
Model.initialValues_Sens_1 = initialValuesSens_1;

% Integrate model equations
[timeData, speciesEstimates]  = ode45(@FitzHughNagumo,...
                                      timePoints,...
                                      initialValues,...      
                                      odeset('RelTol', 1e-6),...
                                      params);

Model.timeData         = timeData' ; % note the transposes here
speciesEstimates       = speciesEstimates' ;
Model.speciesEstimates = speciesEstimates;

[numStates, numTimePts] = size(speciesEstimates);

% Add Noise to trajectories
noiseSD         = Model.addedNoise_SD;
Model.NoisyData = speciesEstimates + ...
                  randn(numStates, numTimePts) .* noiseSD;


%%%%%%%%%%%%%%%%%%%%%%%%%        
% Call sampling routines %
%%%%%%%%%%%%%%%%%%%%%%%%%

RMHMC(Model);

MALA(Model);







