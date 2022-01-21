function dSV = getSV_dot_CR3BP(SV,mu)
% =======================================================================
%            Derivative of Non-dimensional State Vector in CR3BP
% =======================================================================
%
% Author : Kevin Charls (jackcruose47)
%
% Last Update : 10-10-2020
%
% Format : dSV = getSV_dot_CR3BP(SV,mu)
%
% -----------------------------------------------------------------------
%                               INPUTS
% -----------------------------------------------------------------------
% SV            : Non-dimensional state vector in CR3BP [6x1]
% mu            : 3-body constant [1x1]
% -----------------------------------------------------------------------
%
% -----------------------------------------------------------------------
%                              OUTPUTS
% -----------------------------------------------------------------------
% dSV           : Derivative of non-dimensional state vector [6x6]
% -----------------------------------------------------------------------
% 
% -----------------------------------------------------------------------
%                            CHANGE LOG
% -----------------------------------------------------------------------
% 10-10-2020 : Code Created (Non-forced)
% -----------------------------------------------------------------------

% -- State vector
rx = SV(1);
ry = SV(2);
rz = SV(3);
vx = SV(4);
vy = SV(5);
vz = SV(6);

r13=sqrt( (rx+mu)^2 + ry^2 + rz^2 );
r23=sqrt( (rx-(1-mu))^2 + ry^2 + rz^2 );


% -- dynamics in synodic non-dimensional state vector
r_dot = [vx;vy;vz];
v_dot = [
    rx + 2*vy - (1-mu)*(rx+mu)/(r13^3) - mu*(rx-1+mu)/(r23^3);
    ry - 2*vx - (1-mu)*ry/(r13^3) - mu*ry/(r23^3);
    -(1-mu)*rz/(r13^3) - mu*rz/(r23^3)
    ];

% -- Creating the output derivative vector
dSV = [r_dot;v_dot];

end