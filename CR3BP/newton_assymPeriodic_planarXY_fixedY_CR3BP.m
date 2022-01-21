function [SV0,tf] = newton_assymPeriodic_planarXY_fixedY_CR3BP(SV0_0,tf_0,mu,tol)
% =======================================================================
%    Newton-Raphson Method To Estimate Initial Condtion For Periodic
%     Orbit Assymetric and Planar in XY plane with fixed Y-velocity
%             Variable X-position, X-velocity and Y-velocity
%                      <Y is the reference axis>
% =======================================================================
%
% Author : Kevin Charls (jackcruose47)
%
% Last Update : 20-10-2020
%
% Format : [SV0,tf] = newton_assymPeriodic_planarXY_fixedVY_CR3BP(...
%                       SV0_0,tf_0,mu,tol)
%
% Ref : [1] Generating Periodic Orbits In The Circular Restricted Threebody
%           Problem With Applicaiton To Lunar South Pole Coverage
%           - Daniel J. Grebow
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
% 20-10-2020 : Code Created
% -----------------------------------------------------------------------

SV0 = SV0_0;
tf = tf_0;

% -- Variables Vector For Newton-Raphson
% - Selected variables (rx,vx,vy,tf)
X = [SV0(1);SV0(4);SV0(5);tf];

% -- Main loop of Newton-Raphson
while true
    
    % -- Current state vector
    SV0(1) = X(1);
    SV0(4) = X(2);
    SV0(5) = X(3);
    tf = X(4);
    
    % -- Integrating trajectory in CR3BP with STM computation (half period)
    [SV,~,STM] = propagateTrajectory_CR3BP(SV0,tf,mu,0,...
        'ComputeSTM',true,'InitialSTM',eye(6,6));
    
    % -- State vector at the end of computation
    SV = SV(:,end);
    % -- State trasition matrix at end of computation
    STM = STM(:,:,end);
    
    % -- Derivative state vector
    dSV = getSV_dot_CR3BP(SV,mu);
    
    % -- Modified STM (for given y)
    Q = STM - 1/dSV(2)*dSV*STM(2,:);
    
    % -- Building Constrain Vector
    % Final rx-rx0,vx-vx0
    F = [SV(1)-SV0(1); SV(4)-SV0(4); SV(5)-SV0(5)];
    
    % -- Matrix of Partial Derivatie of Constraints wrt Free-Variables
    dF = [Q(1,1)-1, Q(1,4), Q(1,5), dSV(1);
        Q(4,1), Q(4,4)-1, Q(4,5), dSV(4)
        Q(5,1), Q(5,4), Q(5,5)-1, dSV(5)];
    
    % -- Newton-Raphson Step
    X_new = X - dF'/(dF*dF')*F;
    
    % -- Breaking Conditon
    if norm(F)<tol
        break
    end
    
    % -- Updating results
    X = X_new;
end

SV0(1) = X(1);
SV0(4) = X(2);
SV0(5) = X(3);
tf = X(4);
end