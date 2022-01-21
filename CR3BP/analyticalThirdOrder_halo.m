function [SV_uvw,t,T] = analyticalThirdOrder_halo(Au,Aw,phi,L,mu,tf,N,type)
% =======================================================================
%     First Order Analytical Solution to CR3BP Near Colinear Points
%          [Uses Amplitude in U and W motion to define orbit]
%                [Periodic Halo Orbit (Class I/Class II)]
% =======================================================================
%
% Author : Kevin Charls (jackcruose47)
%
% Last Update : 15-10-2020
%
% Format : [SV_uvw,t,T] = analyticalThirdOrder_Halo(Au,Aw,phi,L,...
%                         mu,tf,N,type)
%
% Ref : [1] Analytic Construction of Periodic Orbits about the Collinear 
%           Points - David L. Richardson
%       [2] Periodic and Quasi Periodic Halo Orbits in the Earth-Sun / 
%           Earth-Moon Systems - David L. Richardson
%
% -----------------------------------------------------------------------
%                               INPUTS
% -----------------------------------------------------------------------
% Au            : Amplitude in U-axis [1x1]
% Aw            : Amplitude in W-axis [1x1]
% L             : Required colinear Lagrange point [1x1] (1,2,3)
% mu            : Three-body constant [1x1]
% tf            : Final time point for computation [1x1]
% N             : Number of points to compute [1x1]
% type          : Type of orbit (0 for Class I, 1 for Class II)
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

% -- Time vector
t = linspace(0,tf,N);

% -- Compute Lagrange points
LP = getLagrangePoints_CR3BP(mu);

% -- Selecting Lagrange point of intrest
if L>=1 && L<=3 && round(L)==L
    r_eq = LP(:,L);
else
    fprintf(2,'ERROR : L must be an integer between 1 and 3 !\n');
    fprintf('Program only works for 3 Colinear Points !\n');
    return
end

% -- Compute Udiff matrix
Udiff = getUdiff_CR3BP(r_eq,mu);

% -- Computing all necessary coef.
% - coef. of characteristic equation (1,2)
b_1 = 2 - (Udiff(1,1)+Udiff(2,2))/2;
b_2 = sqrt(- Udiff(1,1)*Udiff(2,2));
% - Linearized frequncy
lamda = sqrt( b_1 + sqrt( b_1^2 + b_2^2  ) );
% - coef. of characteristic equation (3)
k = (lamda^2 + Udiff(1,1))/(2*lamda); % also called beta-3
% - Time period
T = 2*pi/lamda;

% - Computing gamma parameter 
% gamma = distance between nearest primary and lagrange point
if L<3
    g = abs((1-mu)-r_eq(1));
else
    g = -mu-r_eq(1);
end
% - Computing c-coef. from Legendre Polynomial
c = zeros(4:1);
for i = 2:4
    switch L
        case 1
            c(i) = 1/g^3 * ( mu + ((-1)^i)*( (1-mu)*g^(i+1) )/( (1-g)^(1+i) ) );
        case 2 
            c(i) = 1/g^3 * ( ((-1)^i)*mu + ((-1)^i)*( (1-mu)*g^(i+1) )/( (1+g)^(1+i) ) );
        case 3
            c(i) = 1/g^3 * ( 1 - mu + ( mu*g^(i+1) )/( (1+g)^(1+i) ) );
    end
end
c2 = c(2);
c3 = c(3);
c4 = c(4);
% - Computing Delta
Delta = lamda^2 - c(2);
% - Additional coef. aij,bij,dij,di
a21 = 3*c3*(k^2-2)/(4*(1+2*c2));
a22 = 3*c3/(4*(1+2*c2));
d1 = 3*lamda^2/k*(k*(6*lamda^2-1)-2*lamda);
d2 = 8*lamda^2/k*(k*(11*lamda^2-1)-2*lamda);
a23 = -3*c3*lamda/(4*k*d1) * (3*k^3*lamda-6*k*(k-lamda)+4);
a24 = -3*c3*lamda/(4*k*d1) *(2 + 3*k*lamda);
b21 = -3*c3*lamda/(2*d1) * (3*k*lamda-4);
b22 = 3*c3*lamda/d1;
d21 = -c3/(2*lamda^2);
a31 = - 9*lamda/(4*d2)*(4*c3*(k*a23-b21) + k*c4*(4+k^2)) + ...
    ((9*lamda^2+1-c2)/(2*d2))*(3*c3*(2*a23-k*b21)+c4*(2+3*k^2));
a32 = -1/d2*( 9*lamda/4*(4*c3*(k*a24-b22)+k*c4) + ...
     3/2*(9*lamda^2+1-c2)*(c3*(k*b22+d21-2*a24)-c4)  );
b31 = 3/(8*d2)*( 8*lamda*(3*c3*(k*b21-2*a23)-c4*(2+3*k^2)) + ...
    (9*lamda^2+1+2*c2)*(4*c3*(k*a23-b21)+k*c4*(4+k^2)) );
b32 = 1/d2*( 9*lamda*(c3*(k*b22+d21-2*a24)-c4) + ...
    3/8*(9*lamda^2+1+2*c2)*(4*c3*(k*a24-b22)+k*c4) );
d31 = 3/(64*lamda^2)*(4*c3*a24+c4); 
d32 = 3/(64*lamda^2)*(4*c3*(a23-d21)+c4*(4+k^2));
% - Frequency correction ceof.
s1 = 1/(2*lamda*(lamda*(1+k^2)-2*k)) * ( 3/2*c3*(2*a21*(k^2-2)-...
    a23*(k^2+2)-2*k*b21) - 3/8*c4*(3*k^4-8*k^2+8));
s2 = 1/(2*lamda*(lamda*(1+k^2)-2*k)) * ( 3/2*c3*(2*a22*(k^2-2)+...
    a24*(k^2+2)+2*k*b22+5*d21) + 3/8*c4*(12-k^2));
% - amplitude-constraint relationship
a1 = -3/2*c3*(2*a21+a23+5*d21) - 3/8*c4*(12-k^2);
a2 = 3/2*c3*(a24-2*a22)+9/8*c4;
l1 = a1 + 2*lamda^2*s1;
l2 = a2 + 2*lamda^2*s2;

% - Alternative Equations for lamda and k [Approximation with legendre poly]
% k = 2*lamda/(lamda^2+1-c(2))
% lamda^4 + (c(2)-2)*lamda^2-(c(2)-1)*(1+2*c(2)) = 0


% - Converting to correct normalization
% Au originaly is defined normalized to r12 
% Here Au is normalized w.r.t. distance to nearest primary
Au = Au/g;
Aw = Aw/g;


% -- Computing Ax or Az if any are not given
% From amplitude relationship l1*Ax^2+l2*Az^2+Delta = 0
if Aw == 0 && Au>0
    Aw = sqrt(-(Delta+l1*Au^2)/l2);
elseif Au == 0 && Aw>0
    Au = sqrt(-(Delta+l2*Aw^2)/l1);
elseif Au == 0 && Aw == 0
    Au = sqrt(abs(Delta/l1));
else
    fprintf(2,'ERROR : Ax and Az must be non negetive!\n');
    return
end

% -- Overall phase (tau = frequncy*time + phase_shift)
tau = lamda.*t + phi;

% -- Computing Analytical Third Order Position
% - Switching function for Halo-orbit class
switch type
    case 0
        delta = 1;
    case 1
        delta = -1;
end
% - Position in each axis (normalized w.r.t. distance to primary)
u = a21*Au^2 + a22*Aw^2 - Au.*cos(tau) + ...
    (a23*Au^2-a24*Aw^2).*cos(2*tau) + ...
    (a31*Au^3-a32*Au*Aw^2).*cos(3*tau);
v = k*Au.*sin(tau) + (b21*Au^2-b22*Aw^2).*sin(2*tau) + ...
    (b31*Au^3 - b32*Au*Aw^2).*sin(3*tau);
w = delta * (Aw.*cos(tau) + d21*Au*Aw.*(cos(2*tau)-3) + ...
    (d32*Aw*Au^2-d31*Aw^3).*cos(3*tau));
% - Velocity in each axis (normalized w.r.t. distance to primary)
du = Au*lamda.*sin(tau) - (a23*Au^2-a24*Aw^2)*2*lamda.*sin(2*tau) - ...
    (a31*Au^3-a32*Au*Aw^2)*3*lamda.*sin(3*tau);
dv = k*Au*lamda.*cos(tau) + (b21*Au^2-b22*Aw^2)*2*lamda.*cos(2*tau) + ...
    (b31*Au^3 - b32*Au*Aw^2)*3*lamda.*cos(3*tau);
dw = delta * (-Aw*lamda.*sin(tau) - d21*Au*Aw*2*lamda.*sin(2*tau) - ...
    (d32*Aw*Au^2-d31*Aw^3)*3*lamda.*sin(3*tau));


% -- Final UVW state vector (normalized w.r.t. distance to primary)
SV_uvw = [u;v;w;du;dv;dw];

% - Converting values from normalized in g to normalized in r12
SV_uvw = SV_uvw .* g;

end