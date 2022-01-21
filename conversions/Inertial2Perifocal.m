function [perifocal] = Inertial2Perifocal(inertial,inc,RA,omega)
% This function converts any vector in inertial to perifocal.
%
% INPUTS OF FUNCTION
% inertial     : Inertial fram vector
% inclination  : Inclination of orbit
% RA           : Right Accension of accent node
% omega        : The argument of perigee
%
% OUTPUTS OF THE FUNCTION
% perifocal    : Perifocal frame vector

% Creating individual the rotating matrix definitions
% R1 = [ 1 0 0; 0 cos(x) sin(x); 0 -sin(x) cos(x) ];
% R3 = [ cos(x) sin(x) 0; -sin(x) cos(x) 0; 0 0 1];

perifocal = zeros(size(inertial));

N = size(inertial,2);

if length(RA) == 1 && N>1
    RA = ones(1,N)*RA;
    inc = ones(1,N)*inc;
    omega = ones(1,N)*omega;
end

for i = 1:N
    R3_omega = [ cos(omega(i)) sin(omega(i)) 0; -sin(omega(i)) cos(omega(i)) 0; 0 0 1];
    R1_inclination = [ 1 0 0; 0 cos(inc(i)) sin(inc(i)); 0 -sin(inc(i)) cos(inc(i)) ];
    R3_RA = [ cos(RA(i)) sin(RA(i)) 0; -sin(RA(i)) cos(RA(i)) 0; 0 0 1];
    
    % The final overall rotating matrix for Inertial to Perifocal rotation
    Q = R3_omega*R1_inclination*R3_RA;
    
    % Converting Inertial vector to Perifocal vector
    perifocal(:,i) = Q * inertial(:,i);
end

end