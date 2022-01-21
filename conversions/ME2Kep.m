function Kep = ME2Kep(ME)
% =======================================================================
%      Modified Equinoctial Elements to Keplarian Elements Conversion
% =======================================================================
%
% Author : Kevin Charls (jackcruose47)
%
% Last Update : 29-09-2020
%
% Format : [Kep] = ME2Kep(ME,mu)
%
% -----------------------------------------------------------------------
%                               INPUTS
% -----------------------------------------------------------------------
% ME            : Modified Equinoctial Elements [6xN]
% mu            : Gravitational Parameter of primary [1x1]
% -----------------------------------------------------------------------
%
% -----------------------------------------------------------------------
%                              OUTPUTS
% -----------------------------------------------------------------------
% Kep           : Keplarian Elements [6xN]
% -----------------------------------------------------------------------
%
% -----------------------------------------------------------------------
%                            CHANGE LOG
% -----------------------------------------------------------------------
% 29-09-2020 : Code Vectorized
% -----------------------------------------------------------------------

% For SIMULINK CODE GENERATION
Kep = zeros(size(ME));

% - Extracting Modified Equinoctial Elements
p = ME(1,:);
f = ME(2,:);
g = ME(3,:);
h = ME(4,:);
k = ME(5,:);
L = ME(6,:);

% - Converting to Keplarian Elements
e = sqrt(f.^2+g.^2);
a = p./(1-e.^2);
inc = 2.*atan(sqrt(h.^2+k.^2));
omega = atan2((g.*h-f.*k),(f.*h+g.*k));
RA = atan2(k,h);
theta = L - (RA+omega);

% - Final Keplarian Element Vector Formulation
Kep = [a;e;inc;RA;omega;theta];

end