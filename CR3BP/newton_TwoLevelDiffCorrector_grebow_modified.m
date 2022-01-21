function [SVp,tp,DV,Error] = newton_TwoLevelDiffCorrector_grebow_modified(SVp_guess,tp_guess,mu,varargin)
% =======================================================================
%        Newton-Raphson Method Two Level Differential Corrector
%                 (Position Error and Delta-V Reduction)
%             [Finds Natural Orbits close to initial Guess]
%                 < Daniel J. Grebow Modified Version >
%                 External Points are fixed in position
% =======================================================================
%
% Author : Kevin Charls (jackcruose47)
%
% Last Update : 06-11-2020
%
% Format : [SV0,tf,DV] = newton_TwoLevelDiffCorrector_grebow_modified(...
%                        SV0_guess,tp_guess,mu,periodic,varargin)
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
% varargin      : Addition parameters (parameter-value pairs)
%       DiffCorrTolerance1       : Tolerance of Diff. Correction 1 [1x1]
%       IterMax1                 : Max iterations Diff. Corr. 1 [1x1]
%       DiffCorrTolerance2       : Tolerance of Diff. Correction 2[1x1]
%       IterMax2                 : Max iterations Diff. Corr. 1 [1x1]
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
% 06-11-2020 : Code Created
% -----------------------------------------------------------------------

% -- Default settings of parameters
IterMax1 = 10;
IterMax2 = 100;
TolDiff1 = 1e-14;
TolDiff2 = 1e-7;
DVmax = 10;

% - Setting options based on additional inputs
for ipar = 1:2:length(varargin)
    switch varargin{ipar}
        case 'IterMax1'
            IterMax1 = varargin{ipar+1};
        case 'IterMax2'
            IterMax2 = varargin{ipar+1};
        case 'DiffCorrTolerance1'
            TolDiff1 = varargin{ipar+1};
        case 'DiffCorrTolerance2'
            TolDiff2 = varargin{ipar+1};
        case 'MaxDVThreshold'
            DVmax = varargin{ipar+1};
    end
end

% -- Initializing vectors
SVp = SVp_guess;
tp = tp_guess;

% -- Preprocessing Calculations
% - Number of patch points
M = size(tp,2);
% - Number of Segments
N = M - 1;

% - Initialzing error state
Error = 0;

% -- Performing 1st level differential correction on velocity and time
[SVp,tp,DV,DVmag] = diffCorr1(SVp,tp,mu,N,TolDiff1,IterMax1);

if DVmag>DVmax
    Error = 1;
    return;
end

% -- Main loop of Newton Raphson (2nd Level)
iter = 0;
while IterMax2==0 || iter<IterMax2
    
    % -- Perform Single Step of 2nd Level Diff. Corr.
    [SVp,tp,DVfinal,DRfinal] = diffCorr2_step(SVp,tp,mu,M,N);
    
    % -- Perform Complete 1st Level Diff. Corr.

    [SVp,tp,DV,DVmag] = diffCorr1(SVp,tp,mu,N,TolDiff1,IterMax1);
    
    DVmag = DVmag + DVfinal;
    
    % -- Check Delta-V magnitude (break)
    if DVmag<TolDiff2 && DRfinal<TolDiff2
        break;
    end   
    
    if DVmag>DVmax
        Error = 1;
        break;
    end
    
    iter = iter+1;
end


% --- Additional Functions

% - Function to performs complete 1st level diff. corr.
    function [SVp,tp,DV,DVmag] = diffCorr1(SVp,tp,mu,N,tol1,iterMax)
        
        % - Number of constraints for 1st level diff. corr
        ncon = 3*N;
        % - Number of free varaibles for 1st level diff. corr
        nvar = 4*N;
        
        % -- Computing the free-variable vector
        V0 = SVp(4:6,1:end-1);
        tf = tp(2:end) - tp(1:end-1);
        X = [reshape(V0,3*N,1);tf'];
        
        % -- Main loop of Newton Raphson
        k = 0;
        while iterMax==0 || k<iterMax
            
            % -- Initializing varaibles
            % - Constraint Vector
            F = zeros(ncon,1);
            % - Differential Matrix
            dF = zeros(ncon,nvar);
            % - Delta-V vector matrix
            DV = zeros(3,N-1);
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
                % Position error for all segments
                F(3*(n-1)+1:3*(n-1)+3) = SV(1:3) - SVp(1:3,n+1);
                % - Computing elements of derivative matrix
                % - For Partial of End-Position with Initial Velocity
                dF((3*(n-1)+1):(3*(n-1)+3),(3*(n-1)+1):(3*(n-1)+3)) = STM(1:3,4:6);
                % - For Partial of End-Position with Time of Flight
                dF((3*(n-1)+1):(3*(n-1)+3),3*N+n ) = SV(4:6);
                % -- Computing Delta-V
                if n<N
                    DV(:,n) = SVp(4:6,n+1) - SV(4:6);
                    DVmag = DVmag+norm(DV(:,n));
                end
            end
            
            % -- Newton-Raphson Step
            Xnew = X - dF'/(dF*dF')*F;
            
            % -- Breaking Conditon
            if norm(F)<tol1
                break
            end

            X = Xnew;
            
            k = k+1;
        end
        
        X = Xnew;
        
        % -- Updating states with diff. corr.
        tf = X(end-N+1:end)';
        SVp(4:6,1:end-1) = reshape(X(1:end-N),3,N);
        tp = cumsum([0,tf]);
        
    end

% - Function performs single step of 2nd level diff. corr.
    function [SVp,tp,DVfinal,DRfinal] = diffCorr2_step(SVp,tp,mu,M,N)
        
        % -- Initializing Cell array to save segment results
        SVi = cell(N,1);  tf = cell(N,1);  SVf = cell(N,1);
        STMif = cell(N,1);  dSVi = cell(N,1);  dSVf = cell(N,1);
        SVi(:) = {zeros(6,1)};  tf(:) = {0}; SVf(:) = {zeros(6,1)};
        STMif(:) = {zeros(6,6)}; dSVi(:) = {zeros(6,1)}; dSVf(:) = {zeros(6,1)};
        
        % -- Computing required data for all segment with integration
        % - Loop over all segments 
        for n = 1:N
            % - Extracting the velocities from Free-Variables for segment
            SVi{n} = SVp(:,n);
            % - Extracting the time of flight from Free-Variables for segment
            tf{n} = tp(n+1)-tp(n);
            % - Integrating trajectory in CR3BP with STM computation
            [SV,~,STM] = propagateTrajectory_CR3BP(SVi{n},tf{n},...
                mu,0,'ComputeSTM',true,'InitialSTM',eye(6,6));
            % - Final state vector
            SVf{n} = SV(:,end);
            % - State trasition matrix at final state
            STMif{n} = STM(:,:,end);
            % - Derivative of final state
            dSVf{n} = getSV_dot_CR3BP(SVf{n},mu);
            % - Derivative of initial state
            dSVi{n} = getSV_dot_CR3BP(SVi{n},mu);
        end
        
        nvar = M*4;
        ncon = (N-1)*3;
        
        % -- Initializing Free-variable Vector
        X = zeros(nvar,1);
        % -- Initializing Constraint Vector
        F = zeros(ncon,1);
        % -- Initializing Derivative Matrix
        dF = zeros(ncon,nvar);
        
        % -- contraint vector and differential matrix for general case
        % - Looping for N-1 Delta-V points
        for  n = 1: N-1
            % - Index of points of intrest
            n_o = n;   % Last point
            n_p = n+1; % Current point
            n_f = n+2; % Next point
            
            %  - Building up X-vector
            if n == 1
                X(1:12) = [SVp(1:3,n_o);tp(n_o);...
                    SVp(1:3,n_p);tp(n_p);...
                    SVp(1:3,n_f);tp(n_f)];
            else
                X(12+4*(n-2)+1:12+4*(n-2)+4) = [SVp(1:3,n_f);tp(n_f)];
            end
            
            % - Adding necessary constraints
            % Computing constraint vector and Diff. Matrix for each point
            [F_con,dF_con] = diffStep_point(n_o,n_p,SVi,SVf,STMif,dSVf,dSVi);
            % Adding to overall constraint vector
            F(3*(n-1)+1:3*(n-1)+3) = F_con;
            % Adding to overall diff. matrix
            dF(3*(n-1)+1:3*(n-1)+3 , 4*(n-1)+1:4*(n-1)+12) = dF_con;
        end
        
        % -- Removing first point and last point positions
        X([1:3,end-3:end-1]) = [];
        dF(:,[1:3,end-3:end-1]) = [];
        
        DVfinal = 0;
        DRfinal = 0;
        
        Xnew = X - dF'/(dF*dF')*F;
        
        % -- Recovering States from new free-variables
        tp(1) = Xnew(1);
        for  n = 2:M-1
           SVp(1:3,n) = Xnew(4*(n-2)+2:4*(n-2)+4);
           tp(n) = Xnew(4*(n-2)+5);
        end
        tp(M) = Xnew(end);
        
    end

% Function performs single step of 2nd level on single point
    function [F,dF] = diffStep_point(n_o,n_p,SVi,SVf,STMif,dSVf,dSVi)
        % - Constaint Vector
        F = -(SVf{n_o}(4:6)-SVi{n_p}(4:6)); % Delta-V
        % - Derivative Vector
        % STM from current to previus point
        STM_po = inv(STMif{n_o});
        Apo = STM_po(1:3,1:3);
        Bpo = STM_po(1:3,4:6);
        % STM from current to next point
        STM_pf = STMif{n_p};
        Apf = STM_pf(1:3,1:3);
        Bpf = STM_pf(1:3,4:6);
        % Computing intermediate matrices
        Mo  = -inv(Bpo);
        Mto = Bpo\SVi{n_o}(4:6);
        Mp  = Bpo\Apo-Bpf\Apf;
        Mtp = Bpf\Apf*SVi{n_p}(4:6) - Bpo\Apo*SVf{n_o}(4:6) ...
            + dSVi{n_p}(4:6) - dSVf{n_o}(4:6);
        Mf  = inv(Bpf);
        Mtf = -Bpf\SVf{n_p}(4:6);
        % Final derivative vector
        dF = [Mo,Mto,Mp,Mtp,Mf,Mtf];

    end

end