function ME = Kep2ME(Kep)
% =======================================================================
%      Keplarian Elements to Modified Equinoctial Elements Conversion
% =======================================================================
%
% Author : Kevin Charls (jackcruose47)
%
% Last Update : 29-09-2020
%
% Format : [ME] = Kep2ME(Kep,mu)
%
% -----------------------------------------------------------------------
%                               INPUTS
% -----------------------------------------------------------------------
% Kep           : Keplarian Elements [6xN]
% mu            : Gravitational Parameter of primary [1x1]
% -----------------------------------------------------------------------
%
% -----------------------------------------------------------------------
%                              OUTPUTS
% -----------------------------------------------------------------------
% ME            : Modified Equinoctial Elements [6xN]
% -----------------------------------------------------------------------
%
% -----------------------------------------------------------------------
%                            CHANGE LOG
% -----------------------------------------------------------------------
% 29-09-2020 : Code Vectorized
% -----------------------------------------------------------------------


% For SIMULINK CODE GENERATION
ME = zeros(size(Kep));

% - Extracting Keplarian Elements of Orbit
a = Kep(1,:);
e = Kep(2,:);
inc = Kep(3,:);
RA = Kep(4,:);
omega = Kep(5,:);
theta = Kep(6,:);

% - Converting to Equinoctial Elements
p = a.*(1-e.^2);
f = e.*cos(omega+RA);
g = e.*sin(omega+RA);
h = tan(inc./2).*cos(RA);
k = tan(inc./2).*sin(RA);
L = theta+omega+RA;

% - Final ME vector formulation
ME = [p;f;g;h;k;L];

end