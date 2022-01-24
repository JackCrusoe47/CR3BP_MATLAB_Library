function Kep  = SV2Kep(SV,mu)
% =======================================================================
%         Cartesian State Vector to Keplarian Elements Conversion
% =======================================================================
%
% Author : Kevin Charls (jackcruose47)
%
% Last Update : 30-01-2021
%
% Format : [Kep] = SV2Kep(SV,mu)
%
% -----------------------------------------------------------------------
%                               INPUTS
% -----------------------------------------------------------------------
% SV            : Cartesian State Vector [6xN]
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
% 30-01-2021 : Vectorized Correction on RA and NODELINE
% 24-01-2022 : Updated condition for equitorial orbit vectorized
% -----------------------------------------------------------------------

% - Initializing for use in Simulink
Kep = zeros(size(SV));
% - Extracting position and velocity
r = SV(1:3,:);
v = SV(4:6,:);

% - Position magnitude
r_mag=vecnorm(r);
% - Angular momentum vector
h_vector = cross(r,v);
% - Eccentricity vector
e_vector = 1/mu .* cross(v,h_vector) - r./r_mag;

% -- Magnitudes of vectors and semi-major axis
e = vecnorm(e_vector);
h = vecnorm(h_vector);
a = h.^2/mu .* 1./(1-e.^2);

% -- Unit vectors
k_unit = [0; 0; 1].*ones(3,length(e));
i_unit = [1; 0; 0].*ones(3,length(e));

% -- The true anomaly of the orbit
theta = acos( dot(r,e_vector) ./ (r_mag.*e)  );
theta(dot(r,v) < 0) = 2*pi - theta(dot(r,v) < 0);
% - Special case for complex theta 
theta = real(theta);

% -- The Inclination of the orbit
inc = acos ( dot(h_vector,k_unit) ./ h );

% - Cross product of Z-axis and Angular momentum vector
MVECTOR = cross(k_unit,h_vector);

% -- Computing Nodeline vector
% - Initializing Nodeline vector
NODELINE = zeros(3,length(e));
% - Finding index where MVECTOR has nonzero norm
idx = ~(vecnorm(MVECTOR)==0);
NODELINE(:,idx) = MVECTOR(:,idx)./vecnorm( MVECTOR(:,idx) );

% -- Computing Right Ascension
idx1 = NODELINE(2,:)>=0;
idx2 = ~idx1;
RA = zeros(1,length(e));
RA(idx1) = acos( NODELINE(1,idx1) ./ vecnorm(NODELINE(:,idx1)) ); 
RA(idx2) = 2*pi - acos( NODELINE(1,idx2) ./ vecnorm(NODELINE(:,idx2)) );
% - Special condition check for undefined nodeline
RA(isnan(RA)) = 0;


% -- The argument of periapsis
% - Computing condition index
idx1 = dot(NODELINE,e_vector) == 0;
idx2 = dot(NODELINE,e_vector) ~= 0;
% - Computing argument based on conditions
omega = zeros(1,length(e));
omega(idx1) = acos( dot(i_unit(:,idx1),e_vector(:,idx1)) ./ e(idx1) );
omega( e_vector(2,idx1) < 0 ) = 2*pi - omega( e_vector(2,idx1) < 0 );
if any(idx2)
    omega(idx2) = acos( dot(NODELINE(:,idx2),e_vector(:,idx2)) ./ (vecnorm(NODELINE(:,idx2).*e(idx2)) ));
    omega( e_vector(3,idx2) < 0 ) = 2*pi - omega( e_vector(3,idx2) < 0 );
end

% -- Final Keplarian Elements
Kep = [a;e;inc;RA;omega;theta];

end