function [SV_LVLH] = Inertial2LVLH_SV(SV_Inertial,SV)

r = SV(1:3,:);
v = SV(4:6,:);

pos_Inertial = SV_Inertial(1:3);
vel_Inertial = SV_Inertial(4:6);

omega = cross(r,v)/norm(r)^2;

% deriving unit vectors in LVLH frame
i_unit = r/norm(r);
k_unit = cross(r,v)/norm(cross(r,v));
j_unit = cross(k_unit,i_unit);

A_LVLH = [ i_unit' ; j_unit'  ; k_unit' ];

pos_LVLH = A_LVLH * pos_Inertial;
vel_LVLH = A_LVLH * (vel_Inertial - cross(omega,pos_Inertial));

SV_LVLH = [pos_LVLH;vel_LVLH];

end