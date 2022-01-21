function [SV_nd,t_nd,varargout] = propagateTrajectory_CR3BP(SV0_nd,tf_nd,mu_star,N,varargin)
% =======================================================================
%        Circular Restricted Three-Body Problem (CR3BP) Propagator
% =======================================================================
%
% Author : Kevin Charls (jackcruose47)
%
% Last Update : 20-01-2022
%
% Format : [SV_nd,t_nd,varargout] = propagateTrajectory_CR3BP(SV0_nd,tf_nd,...
%                                  mu_star,N,varargin)
%
% -----------------------------------------------------------------------
%                               INPUTS
% -----------------------------------------------------------------------
% SV0_nd        : Initial Non-dimensional 3-Body State Vector [6x1]
% tf_nd         : Non-dimensional Integration time [1x1]
% mu_star       : 3-body constant [1x1]
% N             : Number of integration points [0:auto,>0:manual] [1x1]
% varargin      : additional options accepted as parameter-value pairs
% -----------------------------------------------------------------------
%
% -----------------------------------------------------------------------
%                              OUTPUTS
% -----------------------------------------------------------------------
% SV0_nd        : Matrix of MEE Column Vectors [6xN]
% t_nd          : Row vector of time [1xN]
% varargout     : additional outputs based on additional options
% -----------------------------------------------------------------------
%
% -----------------------------------------------------------------------
%                            CHANGE LOG
% -----------------------------------------------------------------------
% 10-10-2020 : Code Created
% 20-01-2022 : Updated Addon function to have dy outputs
%              - example mass derivative
% -----------------------------------------------------------------------

%% --- Pre-processing

% -- Input health check
Error = errorCheck(SV0_nd,tf_nd,N,varargin);
% - Final Error check
if Error > 0
    fprintf(2,'Solver Terminated with %g Error(s)!\n',Error);
    SV_nd = [];
    t_nd = [];
    return;
end

% -- Extracting parameter value pairs

% - Default options
RelTol = 2.5e-14;               % Relative solver tolerance
AbsTol = 1e-22;                 % Absolute solver tolerance
ode_solver = 'ode113';          % ODE solver to use
FunAccel = [];                  % Handle to acceleration computation function
FunEvents = [];                 % Handle to solver termination events function
FunAddon = [];                  % Handle to add-on function called every iteration
RecoverDerivatives = false;     % Flag to recover derivatives
ComputeSTM = false;             % Flag to compute STM
IntegrationDir = 1;             % Integration Direction (1 : Forward, -1 Backward)

% - Setting options based on additional inputs
for ipar = 1:2:length(varargin)
    switch varargin{ipar}
        case 'Solver'
            ode_solver = varargin{ipar+1};
        case 'RelTol'
            RelTol = varargin{ipar+1};
        case 'AbsTol'
            AbsTol = varargin{ipar+1};
        case 'Accel'
            FunAccel = varargin{ipar+1};
        case 'Events'
            FunEvents = varargin{ipar+1};
        case 'Addon'
            FunAddon = varargin{ipar+1};
        case 'RecoverDerivatives'
            RecoverDerivatives = varargin{ipar+1};
        case 'ComputeSTM'
            ComputeSTM = varargin{ipar+1};
        case 'InitialSTM'
            STM0 = varargin{ipar+1};
        case 'IntegrationDirection'
            IntegrationDir = varargin{ipar+1};
    end
end

% -- time span for integration
switch IntegrationDir
    case 1
        tspan = [0,tf_nd];
    case -1
        tspan = [tf_nd,0];
end
if N>0
    tspan = linspace(tspan(1),tspan(2),N);
end

% -- Case with STM to be computed
if ComputeSTM
    SV0_nd = [SV0_nd; STM0(1,1:6)';STM0(2,1:6)';STM0(3,1:6)';...
        STM0(4,1:6)';STM0(5,1:6)';STM0(6,1:6)'];
end

% -- setting up solver
options = odeset('RelTol',RelTol,'AbsTol',AbsTol,'Events',FunEvents);

%% --- Running solver

% -- Running solver based on type selection
switch ode_solver
    case 'ode45'
        [t_nd , SV_nd] = ode45(...
            @(t,y)dyn(t,y,mu_star,FunAccel,FunAddon,ComputeSTM)...
            ,tspan,SV0_nd,options);
    case 'ode113'
        [t_nd , SV_nd] = ode113(...
            @(t,y)dyn(t,y,mu_star,FunAccel,FunAddon,ComputeSTM)...
            ,tspan,SV0_nd,options);
    case 'ode15s'
        [t_nd , SV_nd] = ode15s(...
            @(t,y)dyn(t,y,mu_star,FunAccel,FunAddon,ComputeSTM)...
            ,tspan,SV0_nd,options);
end
% -- Transposing results
t_nd = t_nd';
SV_nd = SV_nd';

%% --- Post-processing

ptr_out = 1;

% - Outputing STM as 3rd output
if ComputeSTM
    STM = zeros(6,6,length(t_nd));
    STM(1,1:6,:) = SV_nd(7:12,:);
    STM(2,1:6,:) = SV_nd(13:18,:);
    STM(3,1:6,:) = SV_nd(19:24,:);
    STM(4,1:6,:) = SV_nd(25:30,:);
    STM(5,1:6,:) = SV_nd(31:36,:);
    STM(6,1:6,:) = SV_nd(37:42,:);
    SV_nd(7:42,:) = [];
    varargout{ptr_out} = STM;
    ptr_out = ptr_out + 1;
end

% Outputing derivative of states as 3rd output / der. of state and dSTM as
% 4th and 5th outputs
if RecoverDerivatives
    dSV_nd = zeros(6,length(t_nd));
    for n = 1:length(t_nd)
        dSV_nd(:,n) = dyn(t_nd(n),SV_nd(:,n),mu_star,FunAccel,FunAddon,ComputeSTM);
    end
    if ComputeSTM
        dSTM = zeros(6,6,length(t));
        dSTM(1,1:6,:) = dSV_nd(7:12,:);
        dSTM(2,1:6,:) = dSV_nd(13:18,:);
        dSTM(3,1:6,:) = dSV_nd(19:24,:);
        dSTM(4,1:6,:) = dSV_nd(25:30,:);
        dSTM(5,1:6,:) = dSV_nd(31:36,:);
        dSTM(6,1:6,:) = dSV_nd(37:42,:);
        dSV_nd(7:42,:) = [];
        varargout{ptr_out} = dSV_nd;
        ptr_out = ptr_out + 1;
        varargout{ptr_out} = dSTM;
    else
        varargout{ptr_out} = dSV_nd;
    end
end

%% --- Additional Functions

    function dy = dyn(t,y,mu,FunAccel,FunAddon,ComputeSTM)
        % Main dynamics function
        
        % -- State vector
        rx = y(1);
        ry = y(2);
        rz = y(3);
        vx = y(4);
        vy = y(5);
        vz = y(6);
        
        r13=sqrt( (rx+mu)^2 + ry^2 + rz^2 );
        r23=sqrt( (rx-(1-mu))^2 + ry^2 + rz^2 );
        
        % -- Computing acceleration (non-dimensional frame)
        if isempty(FunAccel)
            accel = [0;0;0];
        else
            accel = FunAccel(t,y);
        end
        
        % -- dynamics in synodic non-dimensional state vector
        r_dot = [vx;vy;vz];
        v_dot = [
            rx + 2*vy - (1-mu)*(rx+mu)/(r13^3) - mu*(rx-1+mu)/(r23^3);
            ry - 2*vx - (1-mu)*ry/(r13^3) - mu*ry/(r23^3);
            -(1-mu)*rz/(r13^3) - mu*rz/(r23^3)
            ] + accel;
        
        % -- Creating the output derivative vector
        dy = [r_dot;v_dot];
        
        % -- Computing Addon function
        if ~isempty(FunAddon)
            dy = FunAddon(t,y,dy);
        end
        
        % -- Computing STM function
        if ComputeSTM
            dy = calculatedSTM(y,dy,mu);
        end
        
    end

    function Error = errorCheck(SV0_nd,tf_nd,N,varargs)
        Error = 0;
        % - Initial SV vector
        if ~isnumeric(SV0_nd)
            fprintf('Error : Initial SV values must be double\n');
            Error = Error+1;
        end
        % - Initial time points
        if ~isnumeric(tf_nd)
            fprintf('Error : Final time points tf be double\n');
            Error = Error+1;
        end
        % - Number of points
        if ~isnumeric(N) || N<0
            fprintf('Error : Number of points => 0 of double type\n');
            Error = Error+1;
        end
        % - parameter-value pairs
        if ~isempty(varargs) && mod(length(varargs),2)~=0
            fprintf('Error : additional parameters mu_starst be paired with its values\n');
            Error = Error+1;
        end
    end

    function dy = calculatedSTM(y,dy,mu)
        % - Computes the STM at every time step
        
        % - Current STM
        phi = [y(7:12)';y(13:18)';y(19:24)'; ...
            y(25:30)';y(31:36)';y(37:42)'];
        
        % - Jacobian Matrix of CR3BP vector field
        A = getA_CR3BP(y(1:3),mu);
        
        % - Dynamics of STM
        phi_dot = A*phi;
        
        % - Adding dynamics to result vector
        dy(7:42,1) = zeros(36,1);
        dy(7:12,1) = phi_dot(1,1:6)';
        dy(13:18,1) = phi_dot(2,1:6)';
        dy(19:24,1) = phi_dot(3,1:6)';
        dy(25:30,1) = phi_dot(4,1:6)';
        dy(31:36,1) = phi_dot(5,1:6)';
        dy(37:42,1) = phi_dot(6,1:6)';
    end

end