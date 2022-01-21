function dF = jacobian_CTSE(f,X,eps)
% =======================================================================
%    Numerical Method To Find Jacobian With Complex Taylors Series 
%                       Expansion (CTSE) Method
% =======================================================================
%
% Author : Kevin Charls (jackcruose47)
%
% Last Update : 17-10-2020
%
% Format : dF = jacobian_CTSE(f,X,eps)
%
% NOTE : More accurate than CFD but also more expensive. (No cancellation
% errors typically found in finite difference methods and thus perturbation
% coef. can be as small as machine level accuracy would allow it)
%
% Ref : [1] Using Complex Variables to Estimate Derivatives of Real 
%           Functions - William Squire and George Trapp
%       [2] Multidisciplinary Sensitivity Derivatives Using Complex 
%           Variables - J.C.Newman, III, W.K.Anderson, D.L.Whitfield
%
% -----------------------------------------------------------------------
%                               INPUTS
% -----------------------------------------------------------------------
% f             : function handle to function [function handle] [Mx1]
% X             : Free-variable field [Nx1]
% eps           : perturbation coef. [1x1]
%
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
    % - Perturbing selected element (complex perturbation)
    Xperturb(n) = Xperturb(n) + sqrt(-1)*deltaX(n);
    % - CTSE approximation
    dF(:,n) = -imag(f(Xperturb))/deltaX(n);
end

end