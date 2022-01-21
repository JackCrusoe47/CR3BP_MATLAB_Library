function [SV_peri] = CR3BPRot2Perifocal(SV_cr3bp,mu,l_star,v_star,n_star,theta)
% =======================================================================
%        CR3BP Non-Dimensional Rotational to Dimensional Perifocal
%                           Frame Conversion 
%    (Primary body Centered Perifocal Frame of Secondary Body's Orbit)
% =======================================================================
%
% Author : Kevin Charls (jackcruose47)
%
% Last Update : 11-10-2020
%
% Format : [SV_peri] = CR3BPRot2Perifocal(SV_cr3bp,mu,l_star,v_star,...
%                      n_star,theta)
%
% Ref : [1] Small satellite earth-to-moon direct transfer trajectories 
%           using the CR3BP 
%           - Garrett Levi McMillan
%
% -----------------------------------------------------------------------
%                               INPUTS
% -----------------------------------------------------------------------
% SV_cr3bp      : Non-dimensional rotational state vector in CR3BP [6xN]
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
% SV_peri       : Dimensional inertial state vector [6xN]
% -----------------------------------------------------------------------
% 
% -----------------------------------------------------------------------
%                            CHANGE LOG
% -----------------------------------------------------------------------
% 11-10-2020 : Code Created
% -----------------------------------------------------------------------

% -- Magnitude of position and velocity of barry center
r_barry_mag = l_star*mu;
v_barry_mag = n_star*r_barry_mag;

% -- position and velocity of body 3 in rotating frame (non-dimensional)
r_rot_nd = SV_cr3bp(1:3,:);
v_rot_nd = SV_cr3bp(4:6,:);

% -- position and velocity of body 3 in rotating frame (dimensional)
r_rot = r_rot_nd*l_star;
v_rot = v_rot_nd*v_star;

% -- angular velocity vector
omega = [0; 0; n_star];


% -- Initializing perifocal state vector
SV_peri = zeros(size(SV_cr3bp));

for i = 1:length(theta)
    
    % -- Current angle
    theta_c = theta(i);
    
    % -- Position and velocity vector of barrycenter
    r_barry = r_barry_mag*[cos(theta_c); sin(theta_c); 0];
    v_barry = v_barry_mag*[-sin(theta_c); cos(theta_c); 0];
    
    % -- Perifocal to rotatinal frame rotation matrix
    x_ref = [cos(theta_c);    sin(theta_c);     0];
    y_ref = [-sin(theta_c);   cos(theta_c);     0];
    z_ref = [0;             0;              1];
    A_ref = [x_ref';y_ref';z_ref'];
    
    % -- converting rotational to perifocal frame
    r_peri_barry = A_ref'*r_rot(:,i);
    v_peri_barry = A_ref'*v_rot(:,i) + cross(omega,r_peri_barry);
    
    % -- Converting to primary body centered frame
    SV_peri(1:3,i) = r_peri_barry + r_barry;
    SV_peri(4:6,i) = v_peri_barry + v_barry;
end

end