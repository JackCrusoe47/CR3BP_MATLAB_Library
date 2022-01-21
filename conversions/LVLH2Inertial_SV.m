function [SV_Inertial] = LVLH2Inertial_SV(SV_LVLH,SV)

r = SV(1:3,:);
v = SV(4:6,:);

pos_LVLH = SV_LVLH(1:3);
vel_LVLH = SV_LVLH(4:6);

omega = cross(r,v)/norm(r)^2;

% deriving unit vectors in LVLH frame
i_unit = r/norm(r);
k_unit = cross(r,v)/norm(cross(r,v));
j_unit = cross(k_unit,i_unit);

A_LVLH = [ i_unit' ; j_unit'  ; k_unit' ];

pos_Inertial = A_LVLH' * pos_LVLH;
vel_Inertial = A_LVLH' * vel_LVLH + cross(omega,pos_Inertial);

SV_Inertial = [pos_Inertial;vel_Inertial];

end