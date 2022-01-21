function [LVLH] = Inertial2LVLH(Inertial,SV)

r = SV(1:3,:);
v = SV(4:6,:);

% deriving unit vectors in LVLH frame
i_unit = r/norm(r);
k_unit = cross(r,v)/norm(cross(r,v));
j_unit = cross(k_unit,i_unit);

A_LVLH = [ i_unit' ; j_unit'  ; k_unit' ];

LVLH = A_LVLH * Inertial;

end