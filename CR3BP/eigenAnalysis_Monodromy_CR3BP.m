function [EVal,EVec,nu,k,tau,B,deltaEVP,ID,detM] = eigenAnalysis_Monodromy_CR3BP(M,varargin)
% =======================================================================
%              Performs Eigen Analysis on Mondoromy Matrix
% =======================================================================
%
% Author : Kevin Charls (jackcruose47)
%
% Last Update : 18-10-2020
%
% Format : [EVal,EVec,detM,nu] = eigenAnalysis_Monodromy_CR3BP(M)
%
% Ref : [1] Generating Periodic Orbits In The Circular Restricted Threebody
%           Problem With Applicaiton To Lunar South Pole Coverage
%           - Daniel J. Grebow
%       [2] Three-Dimensional, Periodic, ‘Halo’ Orbits - K. Howell
%       [3] LOW-ENERGY LUNAR TRAJECTORY DESIGN 
%           - Jeffrey S. Parker and Rodney L. Anderson
%
% -----------------------------------------------------------------------
%                               INPUTS
% -----------------------------------------------------------------------
% M             : Mondormy Matrix of Orbit [6x6]
% varargin      : Additional arguments (parameter value pairs)
%   'NumericTolerance'      : Numerical Tolerance for identifying pairs
%   'BifurcationTolerance'  : Tolerance to identify bifurcation
% -----------------------------------------------------------------------
%
% -----------------------------------------------------------------------
%                              OUTPUTS
% -----------------------------------------------------------------------
% EVal          : Eigen Values [6x1]
% EVec          : Corresponding Eigen Vectors [6x6]
% nu            : Stability index Type 1 [1x1] [1] [2]
%                 nu = 1/2(|lamda_U|+|lamda_S|)
% k             : Stability index Type 2 [1x1] [3]
%                 k>2 : unstable, k=2 neutral, k<2 stable
% tau           : Stability time fraction [1x1]
%                 perturbation doubling time if unstable orbit
%                 perturbation hald time if stable orbit
% B             : Flag to identify number of bifurcation [0,1,2]
% deltaEVP      : Vector of difference of Eigen value pairs [2x1]
% ID            : Unique ID based on Eigen structure [4x3]
%                 Excludes the Eigen value pair that is unity due to
%                 periodicity property. Thus the 4 other eigen values
%                 [Sign of Real Values, Is value complex, Sign of Complex]
% detM          : Determinant of Mondoromy Matrix [1x1]
%                 Used to estimate numerical quality on Monodromy
%                 Idealy it must be equal to 1
% -----------------------------------------------------------------------
%
% -----------------------------------------------------------------------
%              PROPERTIES OF MONODROMY MATRIX  EIGAN VALUES
% -----------------------------------------------------------------------
% General Properties
% det(M) = 1 ; prod(EigenValues) = 1;
% eig(M) = eig(M') = eig(inv(M));
% 2 values equal to 1
% For Halo Orbits
% - 2 values real and reciprocals (stable and unstable manifolds)
% - 2 values complex but are conjugate pairs and on unit circle
% For Lyap Orbits
% -----------------------------------------------------------------------
%
% -----------------------------------------------------------------------
%                            CHANGE LOG
% -----------------------------------------------------------------------
% 18-10-2020 : Code Created
% -----------------------------------------------------------------------

% - Settings parameters
TolNumeric = 1e-8;
TolBifurcation = 1e-6;

% - Setting options based on additional inputs
for ipar = 1:2:length(varargin)
    switch varargin{ipar}
        case 'BifurcationTolerance'
            TolBifurcation = varargin{ipar+1};
    end
end

% -- Compute Eigan value and vector
[EVec,diag_eig] = eig(M);
EVal = diag(diag_eig);

% -- Alternative method to identify eigan values
a = 2-trace(M);
b = (a^2-trace(M^2))/2+1;
p = (a+sqrt(a^2-4*b+8))/2;
q = (a-sqrt(a^2-4*b+8))/2;
% % - First four Eigen Values
EVal_a = zeros(4,1);
EVal_a(1) = (-p+sqrt(p^2-4))/2;
EVal_a(2) = (-p-sqrt(p^2-4))/2;
EVal_a(3) = (-q+sqrt(q^2-4))/2;
EVal_a(4) = (-q-sqrt(q^2-4))/2;
% % - Two values are 1 (property of periodic orbits)
% EVal_a(5) = 1; 
% EVal_a(6) = 1;

% - Stable and unstable Eigen Values


nu1 = 1/2*(abs(real(EVal_a(1)))+abs(real(EVal_a(2))));
nu2 = 1/2*(abs(real(EVal_a(3)))+abs(real(EVal_a(4))));
if nu1 > 1 || nu2 > 1
    if nu1>nu2
        lamda_U = max(abs(real(EVal_a(1:2))));
        lamda_S = min(abs(real(EVal_a(1:2))));
        nu = nu1;
    else
        lamda_U = max(abs(real(EVal_a(3:4))));
        lamda_S = min(abs(real(EVal_a(3:4))));
        nu = nu2;
    end
else
    nu = 1;
    lamda_U = 1;
end

% lamda_U = max(abs(real(EVal_a)));
% lamda_S = min(abs(real(EVal_a)));
% % - Stability index (type 1)
% nu = 1/2*(real(lamda_U)+real(lamda_S));

% - Stability index (type 2) [3]
k1 = -p;
k2 = -q;
k = max([abs(real(k1)),abs(real(k2))]);

% -- determinant of M
detM = det(M); % should be very close to 1

% -- Perturbation doubling time(unstable)/Perturbation half life(stable)
% Note : Value is in fraction, actual time is given when multiplied by the
% period of the orbit
if lamda_U>0
    tau = log(2)/log(lamda_U);
else
    tau = 0;
end

% -- Difference in eigen value pairs
deltaEVP = [real(EVal_a(1))-real(EVal_a(2));...
    real(EVal_a(3))-real(EVal_a(4))];

% -- Detecting bifurcation
B = 0;
if abs(deltaEVP(1))<TolBifurcation
    B = B+1;
end
if abs(deltaEVP(2))<TolBifurcation
    B = B+1;
end

% -- Creating unique ID from first 4 eigen values
idxReal = and(abs(real(EVal_a))>TolNumeric , abs(imag(EVal_a))<TolNumeric);
idxImaginary = and(abs(real(EVal_a))<TolNumeric , abs(imag(EVal_a))>TolNumeric);
idxComplex = and(abs(real(EVal_a))>TolNumeric , abs(imag(EVal_a))>TolNumeric);

ID = [idxReal,idxComplex,idxImaginary];

end