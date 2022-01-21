function [F,dF] = computeFdF_symPeriodicPlanes_CR3BP(X,SV0,mu,plane)
% =======================================================================
%       Compute Constraint Vector and Partial Derivative Matrix for
%     Newton-Raphson Method To Estimate Initial Condtion For Periodic
%                   Orbit Symmetric Across Given Plane
% =======================================================================
%
% Author : Kevin Charls (jackcruose47)
%
% Last Update : 17-10-2020
%
% Format : [F,dF] = computeFdF_symPeriodicPlanes_CR3BP(X,SV0,mu,plane)
%
% -----------------------------------------------------------------------
%                               INPUTS
% -----------------------------------------------------------------------
% X             : Free-variable Vector for Newton-Raphson [3x1]
% SV0           : Initial State Vector in CR3BP [6x1]
% mu            : 3-body constant [1x1]
% plane         : Symmetric plane {1,2,3,12,13,23,4} [1x1]
% -----------------------------------------------------------------------
%
% -----------------------------------------------------------------------
%                              OUTPUTS
% -----------------------------------------------------------------------
% F             : Constraint Vector [3x1]
% dF            : Partial Derivative Matrix [3x3]
% -----------------------------------------------------------------------
%
% -----------------------------------------------------------------------
%                            CHANGE LOG
% -----------------------------------------------------------------------
% 17-10-2020 : Code Created
% -----------------------------------------------------------------------

% -- Recover States and Time of Flight
switch plane
    case 1 % X-axis
        % X = [rx;vy;vz,T/2]
        SV0(1) = X(1);
        SV0(5) = X(2);
        SV0(6) = X(3);
    case 12 % XY plane
        % X = [rx;ry;vz;T/2]
        SV0(1) = X(1);
        SV0(2) = X(2);
        SV0(6) = X(3);
    case 13 % XZ plane
        % X = [rx;rz;vy;T/2]
        SV0(1) = X(1);
        SV0(3) = X(2);
        SV0(5) = X(3);
    case 23 % YZ plane
        % X = [ry;rz;vx;T/2]
        SV0(2) = X(1);
        SV0(3) = X(2);
        SV0(4) = X(3);
end
tf2 = X(4);

% -- Integrating trajectory in CR3BP with STM computation (half period)
[SV,~,STM] = propagateTrajectory_CR3BP(SV0,tf2,mu,0,...
    'ComputeSTM',true,'InitialSTM',eye(6,6));

% -- State vector at the end of computation
SV = SV(:,end);

% -- State trasition matrix at end of computation
STM = STM(:,:,end);

% -- Derivative state vector
dSV = getSV_dot_CR3BP(SV,mu);

% -- Compute Constraint and Derivative Matrix based on symmetric plane
switch plane
    case 1 % X-axis 
        % F = [vx(f);ry(f);rz(f)];
        % -- Building Constrain Vector
        F = [SV(4); SV(2); SV(3)];
        % -- Matrix of Partial Derivatie of Constraints wrt Free-Variables
        dF = [STM(4,1)  STM(4,5), STM(4,6), dSV(4);
            STM(2,1), STM(2,5), STM(2,6), dSV(2);
            STM(3,1), STM(3,5), STM(3,6), dSV(3)];
    case 12 % XY plane
        % F = [vx(f);vy(f);rz(f)];
        % -- Building Constrain Vector
        F = [SV(4); SV(5); SV(3)];
        % -- Matrix of Partial Derivatie of Constraints wrt Free-Variables
        dF = [STM(4,1)  STM(4,2), STM(4,6), dSV(4);
            STM(5,1), STM(5,2), STM(5,6), dSV(5);
            STM(3,1), STM(3,2), STM(3,6), dSV(3)];
    case 13 % XZ plane
        % F = [vx(f);vz(f);ry(f)];
        % -- Building Constrain Vector
        F = [SV(4); SV(6); SV(2)];        
        % -- Matrix of Partial Derivatie of Constraints wrt Free-Variables
        dF = [STM(4,1)  STM(4,3), STM(4,5), dSV(4);
            STM(6,1), STM(6,3), STM(6,5), dSV(6);
            STM(2,1), STM(2,3), STM(2,5), dSV(2)];        
    case 23 % YZ plane
        % F = [vy(f);vz(f);rx(f)];
        % -- Building Constrain Vector
        F = [SV(5); SV(6); SV(1)];        
        % -- Matrix of Partial Derivatie of Constraints wrt Free-Variables
        dF = [STM(5,2)  STM(5,3), STM(5,4), dSV(5);
            STM(6,2), STM(6,3), STM(6,4), dSV(6);
            STM(1,2), STM(1,3), STM(1,4), dSV(1)];
end


end