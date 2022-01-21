function Udiff = getUdiff_CR3BP(r_nd,mu)
% =======================================================================
%         Compute Hessian of Pseudo Potential Function in CR3BP
% =======================================================================
%
% Author : Kevin Charls (jackcruose47)
%
% Last Update : 10-10-2020
%
% Format : A = getA_CR3BP(r_nd,mu)
%
% -----------------------------------------------------------------------
%                               INPUTS
% -----------------------------------------------------------------------
% r_nd          : Non-dimensional position in CR3BP [3x1]
% mu            : 3-body constant [1x1]
% -----------------------------------------------------------------------
%
% -----------------------------------------------------------------------
%                              OUTPUTS
% -----------------------------------------------------------------------
% Udiff         : Hessian of Potential [3x3]
% -----------------------------------------------------------------------
% 
% -----------------------------------------------------------------------
%                            CHANGE LOG
% -----------------------------------------------------------------------
% 10-10-2020 : Code Created (No-forcing terms)
% -----------------------------------------------------------------------

rx=r_nd(1);
ry=r_nd(2);
rz=r_nd(3);

% -- computing the necessary distances to primary and secondary
r1=sqrt((rx+mu)^2+ry^2+rz^2);
r2=sqrt((rx-(1-mu))^2+ry^2+rz^2);

% -- computing the double derivates with postion
Uxx = 1 - (1-mu)/r1^3 - mu/r2^3 ...
    + 3*(1-mu)*(rx+mu)^2/r1^5 + 3*mu*(rx+mu-1)^2/r2^5;
Uyy = 1 - (1-mu)/r1^3 - mu/r2^3 ...
    + 3*(1-mu)*ry^2/r1^5 + 3*mu*ry^2/r2^5;
Uzz = -(1-mu)/r1^3 - mu/r2^3 ...
    + 3*(1-mu)*rz^2/r1^5 + 3*mu*rz^2/r2^5;
Uxy = 3*ry*(1-mu)*(rx+mu)/r1^5 + 3*ry*mu*(rx-(1-mu))/r2^5;
Uxz = 3*rz*(1-mu)*(rx+mu)/r1^5 + 3*rz*mu*(rx-(1-mu))/r2^5;
Uyz = 3*ry*rz*(1-mu)/r1^5 + 3*ry*rz*mu/r2^5;

% -- exploiting the symmetry
Uyx = Uxy;
Uzx = Uxz;
Uzy = Uyz;

% -- Final second order differential matrix
Udiff = [Uxx,Uxy,Uxz;Uyx,Uyy,Uyz;Uzx,Uzy,Uzz];

end