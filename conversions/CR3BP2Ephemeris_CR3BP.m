function [SV_ephem,t_ephem,et_ephem] = CR3BP2Ephemeris_CR3BP(SV_cr3bp,t_cr3bp,et0,primary,secondary,frame)
% =======================================================================
%          Non-dimensional CR3BP Rotational to Ephemeris Frame
% =======================================================================
%
% Author : Kevin Charls (jackcruose47)
%
% Last Update : 30-01-2021
%
% Format : [SV_ephem,t_ephem,et_ephem] = CR3BP2Ephemeris_CR3BP(SV_cr3bp,...
%               t_cr3bp,et0,primary,secondary,frame)
%
%
% -----------------------------------------------------------------------
%                               INPUTS
% -----------------------------------------------------------------------
%
% -----------------------------------------------------------------------
%
% -----------------------------------------------------------------------
%                              OUTPUTS
% -----------------------------------------------------------------------
%
% -----------------------------------------------------------------------
% 
% -----------------------------------------------------------------------
%                            CHANGE LOG
% -----------------------------------------------------------------------
% 30-01-2021 : Code Created
% -----------------------------------------------------------------------

G = phy_const("C01");
abcorr = 'none'; % Light time correction

% -- Recovering parameters 

Gm1 = cspice_bodvrd( primary, 'GM', 1 );
Gm2 = cspice_bodvrd( secondary, 'GM', 1 );

m1 = Gm1/G;
m2 = Gm2/G;

m_star = m1+m2;
mu = m2/m_star;

N = length(t_cr3bp);
SV_ephem = zeros(6,N);
t_ephem = zeros(1,N);
et_ephem = zeros(1,N);


et_last = et0;
t_cr3bp_last = 0;
t_ephem_last = 0; 

for n = 1:N
    
    % -- computing secondary body's states at et
    SV2 = cspice_spkezr( secondary , et_last , frame, ...
        abcorr, primary);
    % -- Extracting position and velocity vector
    r2 = SV2(1:3);
    v2 = SV2(4:6);
    % -- Extracting current distance from primary
    r12 = norm(r2);
    
    % -- Compute instantaneous characteristic values
    l_star = r12;
    t_star = sqrt(l_star^3/(G*m_star));
    v_star = l_star/t_star;
    

    % -- position and velocity of body 3 in rotating frame (non-dimensional)
    r_rot_nd = SV_cr3bp(1:3,n);
    v_rot_nd = SV_cr3bp(4:6,n);
    
    % -- primary body position and velocity in rotating frame (non-dim)
    r_p = [-mu;0;0];
    v_p = [0;0;0];
    
    % -- converting to primary centered frame (non-dim)
    r_pc_nd = r_rot_nd + r_p;
    v_pc_nd = v_rot_nd + v_p;
    
    % -- position and velocity  in primary centered frame (dimensional)
    r_pc = r_pc_nd*l_star;
    v_pc = v_pc_nd*v_star;
    
    % -- Attitude matrix between inertial and rotating frame
    X_ref = r2/norm(r2);
    Z_ref = cross(r2,v2)/norm(cross(r2,v2));
    Y_ref = cross(Z_ref,X_ref);
    A_ref = [X_ref';Y_ref';Z_ref']';
    
    % -- Instantaneous angular velocity
    Om = norm(cross(r2,v2))/(norm(r2)^2);
    
    % -- Extracting terms of Attitude matrix
    C11 = A_ref(1,1); C12 = A_ref(1,2); C13 = A_ref(1,3);
    C21 = A_ref(2,1); C22 = A_ref(2,2); C23 = A_ref(2,3);
    C31 = A_ref(3,1); C32 = A_ref(3,2); C33 = A_ref(3,3);
    
    % -- Creating B matrix
    B_ref = [Om*C12, -Om*C11, 0;
        Om*C22, -Om*C21, 0;
        Om*C32, -Om*C31, 0];
    
    % -- Zero matrix
    O_ref = zeros(3,3);
    
    % -- Full transformation matrix
    A_full = [A_ref , O_ref;
        B_ref, A_ref];
    
    % -- Final conversion to inertial ephemeris frame
    SV_ephem(:,n) = A_full * [r_pc;v_pc];
    
    % -- Adding time point [dimensional]
    t_ephem(n) = t_ephem_last + (t_cr3bp(n)-t_cr3bp_last)*t_star;
    
    % -- Adding ephemeris time 
    et_ephem(n) = et0 + t_ephem(n);
    
    
    et_last = et_ephem(n);
    t_cr3bp_last = t_cr3bp(n);
    t_ephem_last = t_ephem(n);
end

end