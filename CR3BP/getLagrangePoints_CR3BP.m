function LP = getLagrangePoints_CR3BP(mu)
% =======================================================================
%                   CR3BP Lagrange Points Computation
% =======================================================================
%
% Code Author : Kevin Charls (jackcruose47) [Copy from Web Source)
% Original Code : Matlab-Monkey
%
% Last Update : 10-10-2020
%
% Format : LP = getLagrangePoints_CR3BP(mu_star)
%
% Ref : [1] https://matlab-monkey.com/celestialMechanics/CRTBP/...
%           LagrangePoints/lagrangePoint
%
% -----------------------------------------------------------------------
%                               INPUTS
% -----------------------------------------------------------------------
% mu            : 3-body constant [1x1]
% -----------------------------------------------------------------------
%
% -----------------------------------------------------------------------
%                              OUTPUTS
% -----------------------------------------------------------------------
% LP            : Matrix of Lagrange Point Vectors [3x5]
% -----------------------------------------------------------------------
%
% -----------------------------------------------------------------------
%                            CHANGE LOG
% -----------------------------------------------------------------------
% 10-10-2020 : Code Created
% -----------------------------------------------------------------------

eps=1-mu;

% -- Initializing Matrix of Lagrange Point Vectors
LP = zeros(3,5);

% -- L-1 Co-ordinates
p_L1=[1, 2*(mu-eps), eps^2-4*eps*mu+mu^2, 2*mu*eps*(eps-mu)+mu-eps, mu^2*eps^2+2*(eps^2+mu^2), mu^3-eps^3];
L1roots=roots(p_L1);
% Estimating correct roots of L1
L1=0;
for i=1:5
    if (L1roots(i) > -mu) && (L1roots(i) < eps)
        L1=L1roots(i);
    end
end
LP(1,1) = L1;


% -- L-2 Co-ordinates
p_L2=[1, 2*(mu-eps), eps^2-4*eps*mu+mu^2, 2*mu*eps*(eps-mu)-(mu+eps), mu^2*eps^2+2*(eps^2-mu^2), -(mu^3+eps^3)];
L2roots=roots(p_L2);
% Estimating correct roots of L2
L2=0;
for i=1:5
    if (L2roots(i) > -mu) && (L2roots(i) > eps)
        L2=L2roots(i);
    end
end
LP(1,2) = L2;


% -- L-3 Co-ordinates
p_L3=[1, 2*(mu-eps), eps^2-4*mu*eps+mu^2, 2*mu*eps*(eps-mu)+(eps+mu), mu^2*eps^2+2*(mu^2-eps^2), eps^3+mu^3];
L3roots=roots(p_L3);
% Estimating correct roots of L3
L3=0;
for i=1:5
    if L3roots(i) < -mu
        L3=L3roots(i);
    end
end
LP(1,3) = L3;


% -- L-4 Co-ordinates
LP(1,4) = 0.5 - mu;
LP(2,4) = sqrt(3)/2;


% -- L-5 Co-ordinates
LP(1,5) = 0.5 - mu;
LP(2,5) = -sqrt(3)/2;

end