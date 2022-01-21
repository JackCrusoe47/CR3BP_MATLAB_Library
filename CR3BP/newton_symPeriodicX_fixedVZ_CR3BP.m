function [SV0,tf] = newton_symPeriodicX_fixedVZ_CR3BP(SV0_0,tf_0,mu,tol)
% =======================================================================
%    Newton-Raphson Method To Estimate Initial Condtion For Periodic
%         Orbit Symmetric Across X-axis with Z-velocity Fixed
%                  Variable X-position and Y-velocity
% =======================================================================
%
% Author : Kevin Charls (jackcruose47)
%
% Last Update : 12-10-2020
%
% Format : [SV0,tf] = newton_symPeriodicX_fixedVZ_CR3BP(SV0_0,tf_0,mu,tol)
%
% -----------------------------------------------------------------------
%                               INPUTS
% -----------------------------------------------------------------------
% SV0_0         : Initial Guess State Vector in CR3BP [6x1]
% tf_0          : Initial Guess Orbit Time period [1x1]
% mu            : 3-body constant [1x1]
% tol           : Tolerance of Newton Raphson
% -----------------------------------------------------------------------
%
% -----------------------------------------------------------------------
%                              OUTPUTS
% -----------------------------------------------------------------------
% SV0           : Corrected Initial State Vector [6xN]
% tf            : Corrected Orbital Time Period [1x1]
% -----------------------------------------------------------------------
%
% -----------------------------------------------------------------------
%                            CHANGE LOG
% -----------------------------------------------------------------------
% 12-10-2020 : Code Created
% -----------------------------------------------------------------------

SV0 = SV0_0;
tf2 = tf_0/2;

% -- Variables Vector For Newton-Raphson
% - Selected variables (rx,vy,tf)
X = [SV0(1);SV0(5);tf2];

% -- Main loop of Newton-Raphson
while true
    
    % -- Current state vector
    SV0(1) = X(1);
    SV0(5) = X(2);
    tf2 = X(3);
    
    % -- Integrating trajectory in CR3BP with STM computation (half period)
    [SV,~,STM] = propagateTrajectory_CR3BP(SV0,tf2,mu,0,...
        'ComputeSTM',true,'InitialSTM',eye(6,6));
    % -- State vector at the end of computation
    SV = SV(:,end);
    % -- State trasition matrix at end of computation
    STM = STM(:,:,end);
    
    % -- Derivative state vector
    dSV = getSV_dot_CR3BP(SV,mu);
    
    % -- Building Constrain Vector
    % Final ry,rz,vx
    F = [SV(2); SV(3); SV(4)];
    
    % -- Matrix of Partial Derivatie of Constraints wrt Free-Variables
    dF = [STM(2,1), STM(2,5), dSV(2);
        STM(3,1), STM(3,5), dSV(3);
        STM(4,1), STM(4,5), dSV(4)];
    
    % -- Newton-Raphson Step
    X_new = X - dF\F;
    
    % -- Breaking Conditon
    if norm(F)<tol
        break
    end
    
    % -- Updating results
    X = X_new;
end

SV0(1) = X_new(1);
SV0(5) = X_new(2);
tf2 = X_new(3);
tf = tf2*2;
end