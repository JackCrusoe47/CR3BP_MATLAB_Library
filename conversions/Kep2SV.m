function SV = Kep2SV(Kep,mu)
% =======================================================================
%         Keplarian Elements to Cartesian State Vector Conversion
% =======================================================================
%
% Author : Kevin Charls (jackcruose47)
%
% Last Update : 29-09-2020
%
% Format : [SV] = Kep2SV(Kep,mu)
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
% SV            : Cartesian State Vector [6xN]
% -----------------------------------------------------------------------
%
% -----------------------------------------------------------------------
%                            CHANGE LOG
% -----------------------------------------------------------------------
% 29-09-2020 : Code Vectorized
% -----------------------------------------------------------------------


% - Initializing for use in Simulink
SV = zeros(size(Kep));
% - Extracting keplarian elements
a = Kep(1,:);
e = Kep(2,:);
inc = Kep(3,:);
RA = Kep(4,:);
omega = Kep(5,:);
theta = Kep(6,:);

% - Converting true anomaly to be positive
theta(theta<0) = 2*pi + theta(theta<0);

% - Computing orbital parameters
ra_orbit = a.*(1+e); % Apoapsis of orbit
rp_orbit = a.*(1-e); % Periapsis of orbit
p_orbit = 2.*ra_orbit.*rp_orbit./(ra_orbit+rp_orbit);
h_orbit = sqrt( p_orbit .* mu );

% - Checking true anomaly of hyperbola/parabola with in limit
theta_inf = zeros(1,length(e));
theta_inf(e>=1) = acos(-1./e(e>=1));
theta0_check = theta;
theta0_check(theta>pi) = 2*pi- theta(theta>pi);
if sum(theta0_check(e>=1)>theta_inf(e>=1))>1
    fprintf(2,'\tERROR : Initial point outside computable limits of hyperbola!\n');
        return;
end

% - Computing co-ordinates in perifocal frame
r_peri_mag = (h_orbit.^2)./mu * 1./( 1 + e .* cos(theta) ); % Radius magnitude in perifocal reference frame
r_peri = [ r_peri_mag.*cos(theta) ; r_peri_mag.*sin(theta) ]; % radius vector in perifocal reference frame
v_radial = mu./h_orbit .* e .* sin(theta); % Magnitude of radial velocity
v_trans = mu./h_orbit .* (1 + e .* cos(theta));% Magnitude of transversal velocity
v_peri = [ (v_radial.*cos(theta) - v_trans.*sin(theta)) ; (v_radial.*sin(theta) + v_trans.*cos(theta)) ]; % velocity in perifocal refernce frame
r_peri = [r_peri(1,:);r_peri(2,:);zeros(1,length(e))]; 
v_peri = [v_peri(1,:);v_peri(2,:);zeros(1,length(e))]; 

% - Converting from perifocal to inertial reference frame
r = Perifocal2Inertial(r_peri,inc,RA,omega);
v = Perifocal2Inertial(v_peri,inc,RA,omega);

% - Final cartesian state vectors
SV = [r;v];

end