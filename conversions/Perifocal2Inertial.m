function [inertial] = Perifocal2Inertial(perifocal,inc,RA,omega)
% =======================================================================
%            Perifocal Frame to Inertial Frame Conversion
% =======================================================================
%
% Author : Kevin Charls (jackcruose47)
%
% Last Update : 07-10-2020
%
% Format : [inertial] = Perifocal2Inertial(perifocal,inc,RA,omega)
%
% -----------------------------------------------------------------------
%                               INPUTS
% -----------------------------------------------------------------------
% perifocal     : Perifocal Vector [3xN]
% inc           : Inclination list [1xN]
% RA            : Right accension list [1xN]
% omega         : Argument of Periapsis [1xN]
% -----------------------------------------------------------------------
%
% -----------------------------------------------------------------------
%                              OUTPUTS
% -----------------------------------------------------------------------
% inertial      : Inertial Vector [3xN]
% -----------------------------------------------------------------------
%
% -----------------------------------------------------------------------
%                            CHANGE LOG
% -----------------------------------------------------------------------
% 22-10-2019 : Code Created
% 07-10-2020 : Code Vectorized
% -----------------------------------------------------------------------

inertial = zeros(size(perifocal));

N = size(perifocal,2);

if length(RA) == 1 && N>1
    RA = ones(1,N)*RA;
    inc = ones(1,N)*inc;
    omega = ones(1,N)*omega;
end

for i = 1:N
    R3_RA = [ cos(-RA(i)) sin(-RA(i)) 0; -sin(-RA(i)) cos(-RA(i)) 0; 0 0 1];
    R1_inc = [ 1 0 0; 0 cos(-inc(i)) sin(-inc(i)); 0 -sin(-inc(i)) cos(-inc(i)) ];
    R3_omega = [ cos(-omega(i)) sin(-omega(i)) 0; -sin(-omega(i)) cos(-omega(i)) 0; 0 0 1];
    
    % The final overall rotating matrix for Perifocal to Inertial rotation
    Q = R3_RA*R1_inc*R3_omega;
    
    % Converting perifocal vector to inertial vector
    inertial(:,i) = Q * perifocal(:,i);
end

end