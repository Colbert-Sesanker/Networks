% Solve for MLE values of trajectories without noise
% if sensitivity equations supplied, 'trajectories' matrix 
% will hold analytic solutions to 1st/2nd order sensitivity 
% equations for each state variable (in the matrix columns)
% 'Equations' are typically just time derivatives, but may
% be sensitivity equations to compute analytic sensitivities 

function [timePointsCvode,  ...
          trajectories, ... 
          varargout] =  ...
                        ...                           
          integrateEquations(  equations,          timePoints,...
                               initialValues,      totalParameters,...
                               equations_AD,       numStates,...        
                               sampledParam_idxs,  numSampledParams... 
                            ) % include modelJacobian
if nargin > 4 
    
    % This uses the sundialsTB Cvode integrator to solve ODEs and compute
    % sensitivities using automatic differentiation see Cvode documentation
    % examples for usage. Note the transposes on InitialValues,
    % sensitivities and trajectories
    
    t0               = timePoints(1); % initial time point
    data.params      = totalParameters;
    data.paramList   = sampledParam_idxs;
    
    % Basic options for integrator
    options      = CVodeSetOptions('UserData',     data,  ...
                                   'RelTol',       1.e-5, ...
                                   'AbsTol',       1.e-8, ...
                                   'LinearSolver', 'Dense');
                           
    CVodeInit(equations_AD, 'BDF', 'Newton', ...
              t0, initialValues', options);
          
    sensVarsInitial = zeros(numStates, numSampledParams);
    
    % forward sensitivity analysis (FSA) options
    FSAoptions = CVodeSensSetOptions('method',     'Simultaneous',...
                                     'ErrControl',  true,...
                                     'ParamField', 'params',...
                                     'ParamList',   sampledParam_idxs,...
                                     'ParamScales', ...
                                      data.params(sampledParam_idxs));
          
    CVodeSensInit(numSampledParams, [],  sensVarsInitial, FSAoptions);
    
    [errorStatus,  timePointsCvode, ...
     trajectories, sensitivityMatrix] = CVode(timePoints(2:end), 'Normal');
    
    if errorStatus == -1
        error(errorStatus);
    end   
   
    % the sensitivities and trajectories are concantenated for the first
    % time point so output is the same as ODE45. Cvode does not include
    % values for the first time point
    
    sensAtFirstPoint = zeros(1, numStates);
    sensitivities    = cell(1, numSampledParams);
    for i = 1 : numSampledParams
        sensitivities{i} = sensitivityMatrix(:, i, :);
        sensitivities{i} = [sensAtFirstPoint ; sensitivities{i}(:,:)' ];  
    end    
    
    trajectories = [initialValues ; trajectories' ]; 
    varargout{1} = sensitivities;
    
    
    % free memory
    CVodeFree;
              
else               
                                                                                
    [timePointsCvode, trajectories] = ...
                                  ...
                  ode45(...
                         equations,... 
                         timePoints,... 
                         [initialValues ],...   
                         odeset('RelTol', 1e-6),...                   
                         totalParameters...
                       );                        
   
    % Ocassionally  trajectories may get imaginary artifacts                   
    if any(imag(trajectories(:)))
        error('Trajectory has imaginary components');
    end % if
    
end % if

end % function          
                  
