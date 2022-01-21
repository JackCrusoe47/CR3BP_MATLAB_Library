function dF = jacobian_FFD(f,X,eps)
% =======================================================================
%    Numerical Method To Find Jacobian With Forward Finite Difference
% =======================================================================
%
% Author : Kevin Charls (jackcruose47)
%
% Last Update : 17-10-2020
%
% Format : dF = jacobian_FFD(f,X,eps)
%
% NOTE : Less accurate than CFD but has twice as less function evaluations
%
% Ref : [1] Iterative Methods for Linear and Nonlinear Equations 
%           - C. T. Kelley
%       [2] Iterative Solution of Nonlinear Equations in Several Variables
%           - J. M Ortega, W. C. Rheinboldt
%
% -----------------------------------------------------------------------
%                               INPUTS
% -----------------------------------------------------------------------
% f             : function handle to function [function handle] [Mx1]
% X             : Free-variable field [Nx1]
% eps           : perturbation coef. [1x1]
% -----------------------------------------------------------------------
%
% -----------------------------------------------------------------------
%                              OUTPUTS
% -----------------------------------------------------------------------
% dF            : Jacobian of function [MxN]
% -----------------------------------------------------------------------
% 
% -----------------------------------------------------------------------
%                            CHANGE LOG
% -----------------------------------------------------------------------
% 17-10-2020 : Code Created
% -----------------------------------------------------------------------
% - Function value at current point 
Fx = f(X);

% - Number of variables in free-variable field
nvar = length(X);

% - Number of constraints 
ncon = length(Fx);

% - Initialzing Jacobian Matrix
dF = zeros(ncon,nvar);

% - Evaluation perturbation of state
deltaX = eps*X;
deltaX(deltaX==0) = eps;

for n=1:nvar
    % - Initializing pertubed State
    Xperturb = X;
    % - Perturbing selected element
    Xperturb(n) = Xperturb(n) + deltaX(n);
    % - Forward finite difference
    dF(:,n) = (f(Xperturb)-Fx)/deltaX(n);
end

end