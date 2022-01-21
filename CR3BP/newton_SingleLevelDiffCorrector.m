function [SVp,tp,DV,DVmag] = newton_SingleLevelDiffCorrector(SVp_guess,tp_guess,mu,tol,periodic)
% =======================================================================
%        Newton-Raphson Method For Periodic and Aperiodic orbits
%    (Single Level Differential Corrector - only position continuity)
%       [Mostly finds Patched Set of Segments/Not Natural Orbits]
% =======================================================================
%
% Author : Kevin Charls (jackcruose47)
%
% Last Update : 13-10-2020
%
% Format : [SVp,tp,DV,DVmag] = newton_SingleLevelDiffCorrector(...
%                               SVp_guess,tp_guess,mu,tol,periodic)
%
% Ref : [1] Numerical determination of Lissajous trajectories in the
%           restricted three-body problem
%           - K.C.Howell and H.J.Pernicka
%       [2] Generating Periodic Orbits In The Circular Restricted Threebody
%           Problem With Applicaiton To Lunar South Pole Coverage
%           - Daniel J. Grebow
%
% -----------------------------------------------------------------------
%                               INPUTS
% -----------------------------------------------------------------------
% SVp_guess     : Initial guess set of patch point state vectors [6xN]
% tp_guess      : Initial guess time points for each patch point [1xN]
% mu            : 3-body constant [1x1]
% tol           : Tolerance of Newton Raphson [1x1]
% periodic      : Check Flag for case of Periodic Orbit [Logical 1x1]
%
% NOTE : For periodic case, the program expects tp vector to have N+1
% values (The last time point correspond to time point at return to first
% patch point). If not periodic, the number of time points must be equal to
% number of patch points. (Programe assumes last patch point is end of
% trajectory)
%
% -----------------------------------------------------------------------
%
% -----------------------------------------------------------------------
%                              OUTPUTS
% -----------------------------------------------------------------------
% SVp           : Corrected set of patch point state vectors [6xN]
% tp            : Corrected set of time points for each patch point [1xN]
% DV            : Delta-V vector (non-dim) [3xM]
% DVmag         : Magnitude of total Delta-V [1x1]
% -----------------------------------------------------------------------
%
% -----------------------------------------------------------------------
%                            CHANGE LOG
% -----------------------------------------------------------------------
% 13-10-2020 : Code Created
% -----------------------------------------------------------------------

% -- Initializing the final outputs
SVp = SVp_guess;
tp = tp_guess;

% -- Preprocessing Calculations
% - Number of patch points
M = size(tp,2) - periodic;
% - Number of Segments
N = M - (1-periodic);

% - Number of constraints for 1st level diff. corr
ncon = 3*N;
% - Number of free varaibles for 1st level diff. corr
nvar = 4*N;

% -- Computing the free-variable vector
if periodic
    V0 = SVp(4:6,:);
    tf = tp(2:end) - tp(1:end-1);
else
    V0 = SVp(4:6,1:end-1);
    tf = tp(2:end) - tp(1:end-1);
end
X = [reshape(V0,3*N,1);tf'];

% -- Main loop of Newton Raphson
while true
    
    % -- Initializing varaibles
    % - Constraint Vector
    F = zeros(ncon,1);
    % - Differential Matrix
    dF = zeros(ncon,nvar);
    % - Delta-V vector matrix
    DV = zeros(3,N);
    % - Delta-V magnitude
    DVmag = 0;
    
    % -- looping through segments
    for n = 1:N
        % - Extracting the velocities from Free-Variables for segment
        SV0 = [SVp(1:3,n);X(3*(n-1)+1:3*(n-1)+3)];
        % - Extracting the time of flight from Free-Variables for segment
        tf = X( 3*N+n );
        % - Integrating trajectory in CR3BP with STM computation
        [SV,~,STM] = propagateTrajectory_CR3BP(SV0,tf,mu,0,...
            'ComputeSTM',true,'InitialSTM',eye(6,6));
        % - State vector at the end of computation
        SV = SV(:,end);
        % - State trasition matrix at end of computation
        STM = STM(:,:,end);
        % - Computing elements of constraint vector
        % For all segments except last
        if n<N+(1-periodic)
            DP = SV(1:3) - SVp(1:3,n+1);
            DV(:,n) = SVp(4:6,n+1) - SV(4:6);
        else
            if periodic
                DP = SV(1:3) - SVp(1:3,1);
                DV(:,n) = SVp(4:6,1) - SV(4:6);
            end
        end
        F(3*(n-1)+1:3*(n-1)+3) = DP;
        DVmag = DVmag+norm(DV(:,n));
        % - Computing elements of derivative matrix
        % - For Partial of End-Position with Initial Velocity
        dF((3*(n-1)+1):(3*(n-1)+3),(3*(n-1)+1):(3*(n-1)+3)) = STM(1:3,4:6);
        % - For Partial of End-Position with Time of Flight
        dF((3*(n-1)+1):(3*(n-1)+3),3*N+n ) = SV(4:6);
    end
    
    % -- Newton-Raphson Step
    Xnew = X - dF\F;
    
    % -- Breaking Conditon
    if norm(F)<tol
        break
    end
    
    X = Xnew;
end

X = Xnew;

% -- Updating states with diff. corr.
if periodic
    tf = X(end-N+1:end)';
    SVp(4:6,1:end) = reshape(X(1:end-N),3,N);
else
    tf = X(end-N+1:end)';
    SVp(4:6,1:end-1) = reshape(X(1:end-N),3,N);
end
tp = cumsum([0,tf]);


end