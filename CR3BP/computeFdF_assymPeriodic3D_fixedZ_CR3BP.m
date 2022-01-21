function [F,dF] = computeFdF_assymPeriodic3D_fixedZ_CR3BP(X,SV0,mu)
% =======================================================================
%       Compute Constraint Vector and Partial Derivative Matrix for
%     Newton-Raphson Method To Estimate Initial Condtion For Periodic
%                       Assymmetric 3D orbit
%                         < Z0 reference >
% =======================================================================
%
% Author : Kevin Charls (jackcruose47)
%
% Last Update : 21-10-2020
%
% Format : [F,dF] = computeFdF_assymPeriodic3D_fixedZ_CR3BP(X,SV0,mu)
%
% Ref : [1] Generating Periodic Orbits In The Circular Restricted Threebody
%           Problem With Applicaiton To Lunar South Pole Coverage
%           - Daniel J. Grebow
%
% -----------------------------------------------------------------------
%                               INPUTS
% -----------------------------------------------------------------------
% X             : Free-variable Vector for Newton-Raphson [3x1]
% SV0           : Initial State Vector in CR3BP [6x1]
% mu            : 3-body constant [1x1]
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
% 21-10-2020 : Code Created
% -----------------------------------------------------------------------

% -- Recover States and Time of Flight
SV0 = X(1:6);
tf = X(7);

% -- Integrating trajectory in CR3BP with STM computation (half period)
[SV,~,STM] = propagateTrajectory_CR3BP(SV0,tf,mu,0,...
    'ComputeSTM',true,'InitialSTM',eye(6,6));

% -- State vector at the end of computation
SV = SV(:,end);

% -- State trasition matrix at end of computation
STM = STM(:,:,end);

% -- Derivative state vector
dSV = getSV_dot_CR3BP(SV,mu);

% -- Modified STM (for given z)
Q = STM - 1/dSV(3)*dSV*STM(3,:);

% -- Compute Constraint vector
% Final rx-rx0 ,rz-rz0,vx-vx0,vy-vy0,vz-vz0
F = [SV(1)-SV0(1); SV(2)-SV0(2); SV(3)-SV0(3); SV(4)-SV0(4); SV(5)-SV0(5); SV(6)-SV0(6)];

% -- Compute Derivative Matrix
dF = [Q(1,1)-1, Q(1,2), Q(1,3), Q(1,4), Q(1,5), Q(1,6), dSV(1);
    Q(2,1), Q(2,2)-1, Q(2,3), Q(2,4), Q(2,5), Q(2,6), dSV(2);
    Q(3,1), Q(3,2), Q(3,3)-1, Q(3,4), Q(3,5), Q(3,6), dSV(3);
    Q(4,1), Q(4,2), Q(4,3), Q(4,4)-1, Q(4,5), Q(4,6), dSV(4);
    Q(5,1), Q(5,2), Q(5,3), Q(5,4), Q(5,5)-1, Q(5,6), dSV(5);
    Q(6,1), Q(6,2), Q(6,3), Q(6,4), Q(6,5), Q(6,6)-1, dSV(6)];

end