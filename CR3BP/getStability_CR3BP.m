function [nu,k,tau,BC,region,lamda] = getStability_CR3BP(M,varargin)
% =======================================================================
%            Computes Stability Indices For Mondoromy Matrix
% =======================================================================
%
% Author : Kevin Charls (jackcruose47)
%
% Last Update : 31-10-2020
%
% Format : [nu,k,tau,lamda] = getStability_CR3BP(M)
%
% Ref : [1] Generating Periodic Orbits In The Circular Restricted Threebody
%           Problem With Applicaiton To Lunar South Pole Coverage
%           - Daniel J. Grebow
%       [2] TRAJECTORY DESIGN IN THE EARTH-MOON SYSTEM AND LUNAR SOUTH POLE
%           COVERAGE - Daniel J. Grebow 
%       [3] Three-Dimensional, Periodic, ‘Halo’ Orbits - K. Howell
%       [4] LOW-ENERGY LUNAR TRAJECTORY DESIGN 
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
% nu            : Stability index Type 1 [1x1] [1] [2]
%                 nu = 1/2(|lamda_U|+|lamda_S|)
%                 nu=1 stable, nu>1 unstable
% k             : Stability index Type 2 [1x1] [3] [4]
%                 k>2 : unstable, k=2 neutral, k<2 stable
% tau           : Stability time fraction [1x1] [4]
%                 perturbation doubling time if unstable orbit
%                 perturbation half time if stable orbit
% BC            : Broucke’s Stability Coef [alpha;beta] [2x1]
% region        : Stability Region on Broucke’s Stability Diagram [1x1]
% lamda         : Analytical Eigen Values [6x1]
% -----------------------------------------------------------------------
%
% -----------------------------------------------------------------------
%                            CHANGE LOG
% -----------------------------------------------------------------------
% 31-10-2020 : Code Created
% -----------------------------------------------------------------------

% - Default prameters
stabilitydiagram = 0;
TolNumeric = 1e-8;

% - Setting options based on additional inputs
for ipar = 1:2:length(varargin)
    switch varargin{ipar}
        case 'BrouckeDiagram'
            stabilitydiagram = varargin{ipar+1};
    end
end

% - Broucke’s Stability coeff.
alpha = 2-trace(M);
beta = 1/2*(alpha^2+2-trace(M^2));

% - Discriminant
D = alpha^2-4*(beta-2);

% - Solution to quadratic formula [ sigma^2-alpha*sigma+(beta-2) = 0 ]
sigma1 = 1/2*(alpha+sqrt(D));
sigma2 = 1/2*(alpha-sqrt(D));

% - Eigen Values of M
lamda = zeros(6,1);
lamda(1) = -(sigma1+sqrt(sigma1^2-4))/2;
lamda(3) = -(sigma2+sqrt(sigma2^2-4))/2;
lamda(2) = 1/lamda(1);
lamda(4) = 1/lamda(3);
lamda(5) = 1;
lamda(6) = 1;

% - Eigan value with max modulus (Typically unstable, unless equal to 1)
lamda_max = max(abs(real(lamda)));

% - Stability Index of Type I
nu = 1/2*(lamda_max + 1/lamda_max);

% - Stability Index of Type II
p = (alpha+sqrt(alpha^2-4*beta+8))/2;
q = (alpha-sqrt(alpha^2-4*beta+8))/2;
k1 = -p;
k2 = -q;
k = max([abs(real(k1)),abs(real(k2))]);

% - Perturbation doubling time(unstable)/Perturbation half life(stable)
tau = log(2)/log(lamda_max);

% - Stability region
if imag(sigma1)<TolNumeric && imag(sigma2)<TolNumeric && abs(sigma1)<=2 && abs(sigma2)<=2
    region = 1; % Localy stable
elseif imag(sigma1)>TolNumeric && imag(sigma2)>TolNumeric && abs(conj(sigma1)-sigma2)<TolNumeric
    region = 2; % Complex Instability
elseif imag(sigma1)<TolNumeric && imag(sigma2)<TolNumeric && sigma1>2 && sigma2<-2
    region = 3; % Even-Odd Instability
elseif imag(sigma1)<TolNumeric && imag(sigma2)<TolNumeric && sigma1<-2 && sigma2<-2
    region = 4; % Even-Even Instability
elseif imag(sigma1)<TolNumeric && imag(sigma2)<TolNumeric && sigma1>2 && sigma2>2
    region = 5; % Odd-Odd Instability
elseif imag(sigma1)<TolNumeric && imag(sigma2)<TolNumeric && abs(sigma1)<2 && sigma2<-2
    region = 6; % Even Semi-Instability
elseif imag(sigma1)<TolNumeric && imag(sigma2)<TolNumeric && sigma1>2 && abs(sigma2)<2
    region = 7; % Odd Semi-Instability
end

BC = [alpha;beta];

if stabilitydiagram
   
    if abs(alpha)<20
        a = -20:0.1:20;
    else
        a = linspace(-abs(alpha)*1.2,abs(alpha)*1.2,1000);
    end
    % - boundary 1
    b1 = a.^2/4+2;
    % - boundary 2
    b2 = +2.*a - 2;
    % - boundary 3
    b3 = -2.*a - 2;
    
    figure();
    hold on;
    grid on;
    xlabel('\alpha');
    ylabel('\beta');
    if abs(alpha)<20
    xlim([-20,20]);
    ylim([-20,60]);
    else
        xlim([-abs(alpha)*1.2,abs(alpha)*1.2]);
        ylim([min(b2),max(b2)]);
    end
    title('Broucke’s Stability Diagram');
    plot(a,b1,'k');
    plot(a,b2,'k');
    plot(a,b3,'k');
    
    plot(alpha,beta,'*r');
end

end