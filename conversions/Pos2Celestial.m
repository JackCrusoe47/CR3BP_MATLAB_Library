function [RA,dec] = Pos2Celestial(pos)
% Inertial position to Celestial Co-ordinates conversion
%
% Author : KEVIN CHARLS (jackcrusoe47) 
%
% Keplarian Elements to Inertial State Vector Conversion
%
% INPUTS FOR THE FUNCTION
% Pos    : The cartesian position vector or unit vector
%
% OUTPUTS OF THE FUNCTION
% RA    : Right ascension of body
% dec   : Declination of body

l = pos(1)/norm(pos);
m = pos(2)/norm(pos);
n = pos(3)/norm(pos);
% - Right Ascension and Declination of Sun
dec = asin(n);
if m>0
    RA = acos(l/cos(dec));
else
    RA = 2*pi - acos(l/cos(dec));
end

end