function [SV0,tf] = continuation_PAL_CR3BP(SV0i,tfi,DeltaS,N,mu,type,varargin)
% =======================================================================
%         Pseudo Arc-Length Continuation for Periodic Orbits
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
%                 1     : Lyaponov/Vertical/Butterfly Family
%                 2     : Halo
%                 3     : Axial
%                 4     : Lyaponov L4/L5 (planar asym. periodic y0 ref.)
%                 5     : Vertical L4/L5 (3D assym. periodic y0 ref.)
%                 6     : Axial L4/L5 (3D assym. periodic z0 ref.)
% varargin      : Addition parameters (parameter-value pairs)
%       PALTolerance            : Tolerance of PAL Newton-Raphson [1x1]
%       DiffCorrTolerance       : Tolerance of Diff. Correction [1x1]
%       IterMax                 : Max iterations of PAL loop [1x1]
%       DirectionalIncrement    : Allows automatic step direction based on
%                                 targeted vector and increment direction
%                                 (true/false)
%       TargetVector            : Index of targeted vector (Based on orbit
%                                 type. Look at X vector for each orbit 
%                                 family) eg: X = [rx,rz,vy,T/2] for lyap
%                                 and halo family. For improvement in
%                                 x-axis direction set value to (1)
%       TargetDirection         : Direction of improvement (+1/-1)
%       
% -----------------------------------------------------------------------
%
% -----------------------------------------------------------------------
%                              OUTPUTS
% -----------------------------------------------------------------------
% SV0           : Continued Orbit(s) Initial State Vector [6xN]
% tf            : Continued Orbit(s) Time Period [1xN]
% -----------------------------------------------------------------------
% 
% -----------------------------------------------------------------------
%                            CHANGE LOG
% -----------------------------------------------------------------------
% 19-10-2020 : Code Created
% 21-10-2020 : Added Assymetric periodic orbits (L4/L5)
% -----------------------------------------------------------------------

% -- Default Settings
TolPAL = 1e-6;
TolDiffCorr = 1e-6;
IterMax = 100;
DirIncrement = false;
TargetVector = 1;
TargetDirection = 1; 

% - Setting options based on additional inputs
for ipar = 1:2:length(varargin)
    switch varargin{ipar}
        case 'PALTolerance'
            TolPAL = varargin{ipar+1};
        case 'DiffCorrTolerance'
            TolDiffCorr = varargin{ipar+1};
        case 'IterMax'
            IterMax = varargin{ipar+1};
        case 'DirectionalIncrement'
            DirIncrement = varargin{ipar+1};
        case 'TargetVector'
            TargetVector = varargin{ipar+1};
        case 'TargetDirection'
            TargetDirection = varargin{ipar+1};
    end
end


% -- Computing Free-Variable Vector
switch type
    case 1 % Lyaponov/
        X = [SV0i(1);SV0i(3);SV0i(5);tfi/2];
        plane = 13;
    case 2 % Halo
        X = [SV0i(1);SV0i(3);SV0i(5);tfi/2];
        x0_last = SV0i(1);
        z0_last = SV0i(3);
        plane = 13;
    case 3 % axial
        X = [SV0i(1);SV0i(5);SV0i(6);tfi/2];
        vy0_last = SV0i(5);
        vz0_last = SV0i(6);
        plane = 1;
    case 4
        X = [SV0i(1);SV0i(2);SV0i(4);SV0i(5);tfi];
    case 5
        X = [SV0i;tfi];
    case 6
        X = [SV0i;tfi];
end

% -- Computing Jacobian of Constraints
if type<=3
    [~,dF] = computeFdF_symPeriodicPlanes_CR3BP(X,SV0i,mu,plane);
else
    switch type
        case 4
            [~,dF] = computeFdF_assymPeriodicLyap_CR3BP(X,SV0i,mu);
        case 5
            [~,dF] = computeFdF_assymPeriodic3D_fixedY_CR3BP(X,SV0i,mu); 
        case 6
            [~,dF] = computeFdF_assymPeriodic3D_fixedZ_CR3BP(X,SV0i,mu); 
    end
end

% -- Computing the taject vector (Null Space of Jacobian)
Xdot = null(dF);

% -- Initializing outputs
SV0 = zeros(6,N);
tf = zeros(1,N);

% -- Looping to compute each orbit
for n = 1:N
   
    % -- Initial Guess for New Free-Varaible Vector
    if DirIncrement
        DeltaX = DeltaS*Xdot;
        if TargetDirection*DeltaX(TargetVector)>0
            Xnew = X + DeltaS*Xdot;
        else
            Xnew = X - DeltaS*Xdot;
        end
    else
        Xnew = X + DeltaS*Xdot;
    end
    
    % -- Correcting Guess with Newton-Raphson (PAL contraint)
    for iter = 1:IterMax
        % - Computing new value of function and jacobian
        if type<=3
            [F,dF] = computeFdF_symPeriodicPlanes_CR3BP(X,SV0i,mu,plane);
        else
            switch type
                case 4
                    [F,dF] = computeFdF_assymPeriodicLyap_CR3BP(X,SV0i,mu);
                case 5
                    [F,dF] = computeFdF_assymPeriodic3D_fixedY_CR3BP(X,SV0i,mu); 
                case 6
                    [F,dF] = computeFdF_assymPeriodic3D_fixedZ_CR3BP(X,SV0i,mu); 
            end
        end
        % - Combained Constrain Vector
        G = [F;dot((Xnew-X),Xdot)-DeltaS];
        % - Combained Derivate Matrix
        dG = [dF;Xdot'];
        if norm(F)<TolPAL
            break
        end
        % - Newton Raphson Step
        Xnew = Xnew - dG\G;
    end
    
    % -- Updating values
    Xdot = null(dF);
    X = Xnew;
    
    % -- Recovering States and Differential Correction
    switch type
        
        case 1
            % - Rebuilding state vector
            SV0i(1) = X(1); SV0i(3) = X(2); SV0i(5) = X(3); tfi = 2*X(4);
            % - Periforming differential correction on state
            [SV0i,tfi] = newton_symPeriodicXZ_fixedX_CR3BP(SV0i,tfi,mu,TolDiffCorr);
            % - Rebuilding free-variable vector
            X = [SV0i(1);SV0i(3);SV0i(5);tfi/2];
            
        case 2
            % - Rebuilding state vector
            SV0i(1) = X(1); SV0i(3) = X(2); SV0i(5) = X(3); tfi = 2*X(4);
            % - Periforming differential correction on state
            x0 = SV0i(1);
            z0 = SV0i(3);
            if abs(x0-x0_last)>abs(z0-z0_last)
                [SV0i,tfi] = newton_symPeriodicXZ_fixedX_CR3BP(SV0i,tfi,mu,TolDiffCorr);
            else
                [SV0i,tfi] = newton_symPeriodicXZ_fixedZ_CR3BP(SV0i,tfi,mu,TolDiffCorr);
            end
            x0_last = SV0i(1);
            z0_last = SV0i(3);
            % - Rebuilding free-variable vector
            X = [SV0i(1);SV0i(3);SV0i(5);tfi/2];
            
        case 3
            % - Rebuilding state vector
            SV0i(1) = X(1); SV0i(5) = X(2); SV0i(6) = X(3); tfi = 2*X(4);
            % - Periforming differential correction on state
            vy0 = SV0i(5);
            vz0 = SV0i(6);
            if abs(vy0-vy0_last)>abs(vz0-vz0_last)
                [SV0i,tfi] = newton_symPeriodicX_fixedVY_CR3BP(SV0i,tfi,mu,TolDiffCorr);
            else
                [SV0i,tfi] = newton_symPeriodicX_fixedVZ_CR3BP(SV0i,tfi,mu,TolDiffCorr);
            end
            vy0_last = SV0i(5);
            vz0_last = SV0i(6);
            % - Rebuilding free-variable vector
            X = [SV0i(1);SV0i(5);SV0i(6);tfi/2];
            
        case 4
            % - Rebuilding state vector
            SV0i(1) = X(1); SV0i(2) = X(2); SV0i(4) = X(3); SV0i(5) = X(4); tfi = X(5);
            % - Periforming differential correction on state
            [SV0i,tfi] = newton_assymPeriodic_planarXY_fixedY_CR3BP(SV0i,tfi,mu,TolDiffCorr);
            % - Rebuilding free-variable vector
            X = [SV0i(1);SV0i(2);SV0i(4);SV0i(5);tfi];
            
        case 5
            % - Rebuilding state vector
            SV0i = X(1:6); tfi = X(7);
            % - Periforming differential correction on state
            [SV0i,tfi] = newton_assymPeriodic_3DXY_fixedY_CR3BP(SV0i,tfi,mu,TolDiffCorr);
            % - Rebuilding free-variable vector
            X = [SV0i;tfi];
            
        case 6
            % - Rebuilding state vector
            SV0i = X(1:6); tfi = X(7);
            % - Periforming differential correction on state
            [SV0i,tfi] = newton_assymPeriodic_3DXY_fixedZ_CR3BP(SV0i,tfi,mu,TolDiffCorr);
            % - Rebuilding free-variable vector
            X = [SV0i;tfi];
    end
    
    % -- Adding result to output vector
    SV0(:,n) = SV0i;
    tf(:,n) = tfi;
end

end