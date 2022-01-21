function [SV_cr3bp] = Perifocal2CR3BPRot(SV_peri,mu,l_star,v_star,n_star,theta)
% =======================================================================
%         Dimensional Perifocal to CR3BP Non-Dimensional Rotational
%                           Frame Conversion 
%    (Primary body Centered Perifocal Frame of Secondary Body's Orbit)
% =======================================================================
%
% Author : Kevin Charls (jackcruose47)
%
% Last Update : 30-01-2021
%
% Format : [SV_cr3bp] = Perifocal2CR3BPRot(SV_peri,mu,l_star,...
%                           v_star,n_star,theta)
%
% Ref : [1] Small satellite earth-to-moon direct transfer trajectories 
%           using the CR3BP 
%           - Garrett Levi McMillan
%
% -----------------------------------------------------------------------
%                               INPUTS
% -----------------------------------------------------------------------
% SV_peri       : Dimensional inertial state vector [6xN]
% mu            : 3-body constant [1x1]
% l_star        : Characteristic length [1x1]
% v_star        : Characteristic velocity [1x1]
% n_star        : Characteristic angular velocity [1x1]
% theta         : Current ture anomaly in Perifocal frame [1xN]
% -----------------------------------------------------------------------
%
% -----------------------------------------------------------------------
%                              OUTPUTS
% -----------------------------------------------------------------------
% SV_cr3bp      : Non-dimensional rotational state vector in CR3BP [6xN]
% -----------------------------------------------------------------------
% 
% -----------------------------------------------------------------------
%                            CHANGE LOG
% -----------------------------------------------------------------------
% 27-01-2021 : Code Created
% 30-01-2021 : Velocity corrected
% -----------------------------------------------------------------------

% -- Magnitude of position and velocity of barry center
r_barry_mag = l_star*mu;
v_barry_mag = n_star*r_barry_mag;

% -- position and velocity of body 3 in perifocal frame (dimensional)
r_peri = SV_peri(1:3,:);
v_peri = SV_peri(4:6,:);

% -- angular velocity vector
omega = [0; 0; n_star];

% -- Initializing perifocal state vector
SV_cr3bp = zeros(size(SV_peri));

for i = 1:length(theta)
    
    % -- Current angle
    theta_c = theta(i);
    
    % -- Position and velocity vector of barrycenter
    r_barry = r_barry_mag*[cos(theta_c); sin(theta_c); 0];
    v_barry = v_barry_mag*[-sin(theta_c); cos(theta_c); 0];
    
    % -- Position and velocity vector in barrycenter inertial
    r_peri_barry = r_peri - r_barry;
    v_peri_barry = v_peri - v_barry;
    
    % -- Perifocal to rotatinal frame rotation matrix
    x_ref = [cos(theta_c);    sin(theta_c);     0];
    y_ref = [-sin(theta_c);   cos(theta_c);     0];
    z_ref = [0;             0;              1];
    A_ref = [x_ref';y_ref';z_ref'];
    
    % -- converting rotational to perifocal frame
    r_rot = A_ref*r_peri_barry(:,i);
    v_rot = A_ref*v_peri_barry(:,i) - cross(omega,r_rot);
    
    % -- converting to non-dimensinal units
    SV_cr3bp(1:3,i) = r_rot/l_star;
    SV_cr3bp(4:6,i) = v_rot/v_star;
end

end