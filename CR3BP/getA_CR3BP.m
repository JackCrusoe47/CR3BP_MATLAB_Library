function A = getA_CR3BP(r_nd,mu)
% =======================================================================
%           Jacobian Matrix of CR3BP vector field [STM_dot = A*STM]
% =======================================================================
%
% Author : Kevin Charls (jackcruose47)
%
% Last Update : 10-10-2020
%
% Format : A = getA_CR3BP(r_nd,mu_star)
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
% A             : State Matrix for STM dynamics [6x6]
% -----------------------------------------------------------------------
% 
% -----------------------------------------------------------------------
%                            CHANGE LOG
% -----------------------------------------------------------------------
% 10-10-2020 : Code Created (No-forcing terms)
% -----------------------------------------------------------------------

% -- Computing Hetian of Potential(U)
Udiff = getUdiff_CR3BP(r_nd,mu);

% -- Omega matrix
Omega = [0,2,0; -2,0,0; 0,0,0];

% -- Top rows of G-Matrix
O = zeros(3,3);
I = eye(3,3);

% -- Filling G Matrix
A = [ O,     I;
      Udiff, Omega];

end