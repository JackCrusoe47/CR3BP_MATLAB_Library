function [Inertial] = LVLH2Inertial(LVLH,SV)

r = SV(1:3,:);
v = SV(4:6,:);

% deriving unit vectors in LVLH frame
i_unit = r/norm(r);
k_unit = cross(r,v)/norm(cross(r,v));
j_unit = cross(k_unit,i_unit);

A_LVLH = [ i_unit' ; j_unit'  ; k_unit' ];

Inertial = A_LVLH' * LVLH;

end