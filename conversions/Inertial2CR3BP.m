function SV_cr3bp = Inertial2CR3BP(SV_inert,mu,l_star,v_star,n_star,theta,inc,RA,omega)


% -- Converting inertial frame to perifocal frame
SV_peri = zeros(size(SV_inert));
% - Converting position
SV_peri(1:3,:) = Inertial2Perifocal(SV_inert(1:3,:),inc,RA,omega);
% - Converting velocity
SV_peri(4:6,:) = Inertial2Perifocal(SV_inert(4:6,:),inc,RA,omega);

% -- Converting secondary perifocal to CR3BP rotational frame
[SV_cr3bp] = Perifocal2CR3BPRot(SV_peri,mu,l_star,v_star,n_star,theta);

end

