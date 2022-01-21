function [SV_uvw,t,T] = analyticalFirstOrder_positionUVW(pos0_uvw,L,mu,tf,N,type)
% =======================================================================
%     First Order Analytical Solution to CR3BP Near Colinear Points
%          [Uses initial position U,V and W to define orbit]
%                     [Lissajous/Halo Orbit Case]
% =======================================================================
%
% Author : Kevin Charls (jackcruose47)
%
% Last Update : 15-10-2020
%
% Format : [SV_uvw,t] = analyticalFirstOrder_positionUVW(pos0_uvw,...
%                       L,mu,tf,N)
%
% Ref : [1] Generating Periodic Orbits In The Circular Restricted Threebody
%           Problem With Applicaiton To Lunar South Pole Coverage
%           - Daniel J. Grebow
%
% -----------------------------------------------------------------------
%                               INPUTS
% -----------------------------------------------------------------------
% pos0_uvw      : Initial Position in UVW [1x1]
% L             : Required colinear Lagrange point [1x1] (1,2,3)
% mu            : Three-body constant [1x1]
% tf            : Final time point for computation [1x1]
% N             : Number of points to compute
% type          : Type of orbit (0 for Lissajous, 1 for Halo)
% -----------------------------------------------------------------------
%
% -----------------------------------------------------------------------
%                              OUTPUTS
% -----------------------------------------------------------------------
% SV_uvw        : State Vectors in UVW relative coordinates [6xN]
% t             : Time points [1xN]
% T             : Time period of motion [1x1]
% -----------------------------------------------------------------------
%
% -----------------------------------------------------------------------
%                            CHANGE LOG
% -----------------------------------------------------------------------
% 15-10-2020 : Code Created
% -----------------------------------------------------------------------

% -- Compute Hessian of Potential at given equilibrium point
% - Compute Lagrange points
LP = getLagrangePoints_CR3BP(mu);
% - Selecting Lagrange point of intrest
if L>=1 && L<=3 && round(L)==L
    r_eq = LP(:,L);
else
    fprintf(2,'ERROR : L must be an integer between 1 and 3 !\n');
    fprintf('Program only works for 3 Colinear Points !\n');
    return
end
% - Compute Udiff matrix
Udiff = getUdiff_CR3BP(r_eq,mu);

% -- Computing preliminary parameters
% - z-axis oscilation frequency
nu = sqrt(abs(Udiff(3,3))); 
% - coef. of characteristic equation (1,2)
b_1 = 2 - (Udiff(1,1)+Udiff(2,2))/2;
b_2 = sqrt(- Udiff(1,1)*Udiff(2,2));
% - xy oscillation frequency
s = sqrt( b_1 + sqrt( b_1^2 + b_2^2  ) );
% - coef. of characteristic equation (3)
k = (s^2 + Udiff(1,1))/(2*s);
% - Time period of motion
T = 2*pi/s;

% - Check for halo
if type==1
    % Valid only if Delta = nu-s << (nu + s)/2
    if abs(nu-s)>0.5*(nu+s)*0.01
        fprintf(2,'WARNING : In plane frequency very different from out of plane !\n');
    end
    nu = s; % all axis frequency equal for halo orbit case
    
end

% -- Descritizing time points
t = linspace(0,tf,N);

% -- Extracting initial positions
u0 = pos0_uvw(1);
v0 = pos0_uvw(2);
w0 = pos0_uvw(3);

% -- Computing analytical solutions
% - Position
u = u0.*cos(s.*t)+v0/k.*sin(s.*t);
v = v0.*cos(s.*t)+u0*k.*sin(s.*t);
w = w0.*cos(nu.*t)+w0/nu.*sin(nu.*t);
% - Velocity
du = -u0*s.*sin(s.*t)+v0/k*s.*cos(s.*t);
dv = -v0*s.*sin(s.*t)+u0*k*s.*cos(s.*t);
dw = -w0*nu.*sin(nu.*t)+w0.*cos(nu.*t);
% -- Final UVW state vector
SV_uvw = [u;v;w;du;dv;dw];

end