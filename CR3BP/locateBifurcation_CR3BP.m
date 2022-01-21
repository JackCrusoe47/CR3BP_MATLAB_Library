function [SV0,tf,EVal,EVec,nu,k,tau,B] = locateBifurcation_CR3BP(SV0i,tfi,DeltaS,mu,type)
% =======================================================================
%         Bifurcation Location in Family with Bisection Method
%                Continuation with Pseudo Arc Length
% =======================================================================
%
% Author : Kevin Charls (jackcruose47)
%
% Last Update : 19-10-2020
%
% Format : [SV0,tf] = continuation_PAL(SV0i,tfi,DeltaS,N,mu,type,varargin)
%
% Ref : [1] Generating Periodic Orbits In The Circular Restricted Threebody
%           Problem With Applicaiton To Lunar South Pole Coverage
%           - Daniel J. Grebow
%
% -----------------------------------------------------------------------
%                               INPUTS
% -----------------------------------------------------------------------
% SV0i          : Initial Orbit Initial State Vector [6x1]
% tfi           : Initial Orbit Time Period [1x1]
% DeltaS        : Pseudo Arc-Length step size [1x1]
% N             : Number of orbits to compute [1x1]
% mu            : 3-body constant [1x1]
% type          : type of orbit [1x1]
%                 1     : Lyaponov
%                 2     : Halo
% varargin      : Addition parameters (parameter-value pairs)
%                 PALTolerance      : Tolerance of PAL Newton-Raphson [1x1]
%                 DiffCorrTolerance : Tolerance of Diff. Correction [1x1]
%                 IterMax           : Max iterations of PAL loop [1x1]
% -----------------------------------------------------------------------
%
% -----------------------------------------------------------------------
%                              OUTPUTS
% -----------------------------------------------------------------------
% SV0           : Initial State of Bifurcation Orbit [6x1]
% tf            : Time Period of Bifurcation Orbit [1x1]
% -----------------------------------------------------------------------
% 
% -----------------------------------------------------------------------
%                            CHANGE LOG
% -----------------------------------------------------------------------
% 19-10-2020 : Code Created
% -----------------------------------------------------------------------

IterMax = 100;
TolBisection = 1e-6;

% ------ INITIAL ORBIT -------

% -- Compute Initial Orbit Monodromy
[~,~,STM] = propagateTrajectory_CR3BP(SV0i,tfi,mu,0,...
        'ComputeSTM',true,'InitialSTM',eye(6,6));
Mi = STM(:,:,end);

% ------ FINAL ORBIT -------

% -- Compute Initial Orbit Eigen ID
[~,~,~,~,~,~,~,IDi] = eigenAnalysis_Monodromy_CR3BP(Mi);

% -- Find Final orbit after DeltaS shift
[SV0f,tff] = continuation_PAL_CR3BP(SV0i,tfi,DeltaS,1,mu,type);

% -- Compute Final Orbit Monodromy
[~,~,STM] = propagateTrajectory_CR3BP(SV0f,tff,mu,0,...
        'ComputeSTM',true,'InitialSTM',eye(6,6));
Mf = STM(:,:,end);

% -- Compute Final Orbit Eigen ID
[~,~,~,~,~,~,~,IDf] = eigenAnalysis_Monodromy_CR3BP(Mf);

% ------ BISECTION LOOP -------

if ~isequal(IDi,IDf)
    for iter = 1:IterMax
        
        % - Reduce shift by half
        DeltaS = DeltaS*0.5;
        
        % ------ MID ORBIT -------
        
        % - Find new orbit at DeltaS/2 shift
        [SV0m,tfm] = continuation_PAL_CR3BP(SV0i,tfi,DeltaS,1,mu,type);
        
        % -- Compute Mid Orbit Monodromy
        [~,~,STM] = propagateTrajectory_CR3BP(SV0m,tfm,mu,0,...
            'ComputeSTM',true,'InitialSTM',eye(6,6));
        Mm = STM(:,:,end);
        
        % -- Compute Unique ID of mid orbit
        [EVal,EVec,nu,k,tau,B,~,IDm] = eigenAnalysis_Monodromy_CR3BP(Mm);
        
        
        % ----- BREAKING CONDITION -----
        if B>0 && (norm(SV0i-SV0m)<TolBisection || norm(SV0f-SV0m)<TolBisection)
            SV0 = SV0m;
            tf = tfm;
            break;
        end
        
        % ------ BISECTION STEP -------
        
        % -- Compare with Initial and Final Orbit
        if isequal(IDm,IDi)
            % - bifurcation in section m-f
            SV0i = SV0m;
            tfi = tfm;
            IDi = IDm;
        else
            % - bifurcation in section i-m
            SV0f = SV0m;          
        end
    end
else
    fprintf(2,'WARNING : NO BISECTION DETECTED IN RANGE!');
    SV0 = SV0i;
    tf = tfi;
end

end