function SV_inert = CR3BP2Inertial(SV_cr3bp,mu,l_star,v_star,n_star,theta,inc,RA,omega)

% -- Converting CR3BP to secondary perifocal frame
[SV_peri] = CR3BPRot2Perifocal(SV_cr3bp,mu,l_star,v_star,n_star,theta);


% -- Converting perifocal to inertial frame
SV_inert = zeros(size(SV_peri));
% - Converting position
SV_inert(1:3,:) = Perifocal2Inertial(SV_peri(1:3,:),inc,RA,omega);
% - Converting velocity
SV_inert(4:6,:) = Perifocal2Inertial(SV_peri(4:6,:),inc,RA,omega);

end

