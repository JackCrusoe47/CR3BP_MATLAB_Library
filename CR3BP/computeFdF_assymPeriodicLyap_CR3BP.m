function [F,dF] = computeFdF_assymPeriodicLyap_CR3BP(X,SV0,mu)
% =======================================================================
%       Compute Constraint Vector and Partial Derivative Matrix for
%     Newton-Raphson Method To Estimate Initial Condtion For Periodic
%                    Assymmetric Lyapanov orbit
%                          < Y0 reference > 
% =======================================================================
%
% Author : Kevin Charls (jackcruose47)
%
% Last Update : 20-10-2020
%
% Format : [F,dF] = computeFdF_assymPeriodicLyap_CR3BP(X,SV0,mu)
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
% 20-10-2020 : Code Created
% -----------------------------------------------------------------------

% -- Recover States and Time of Flight
SV0(1) = X(1);
SV0(2) = X(2);
SV0(4) = X(3);
SV0(5) = X(4);
tf = X(5);

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

% -- Compute Constraint vector
% Final rx-rx0, ry-ry0, vx-vx0, vy-vy0
F = [SV(1)-SV0(1); SV(2)-SV0(2); SV(4)-SV0(4); SV(5)-SV0(5)];

% -- Compute Derivative Matrix
dF = [Q(1,1)-1, Q(1,2), Q(1,4), Q(1,5), dSV(1);
    Q(2,1), Q(2,2)-1, Q(2,4), Q(2,5), dSV(2);
    Q(4,1), Q(4,2), Q(4,4)-1, Q(4,5), dSV(4);
    Q(5,1), Q(5,2), Q(5,4), Q(5,5)-1, dSV(5)];

end