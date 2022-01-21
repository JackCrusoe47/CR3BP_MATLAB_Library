function [TNW] = Inertial2TNW(Inertial,SV)

% TNW frame
% x-axis : along velocity 
% y-axis : along cross of velocity and angular momentum
% z-axis : along anglular momentum

r = SV(1:3,:);
v = SV(4:6,:);

% deriving unit vectors in TNW frame
i_unit = v/norm(v);
k_unit = cross(r,v)/norm(cross(r,v));
j_unit = cross(i_unit,k_unit);

A_TNW = [ i_unit' ; j_unit'  ; k_unit' ];

TNW = A_TNW * Inertial;

end