function [M,E] = Theta2Mean(Theta,e)
%Function to calculate the mean and eccentric anomaly from the knowledge of
%the true anomaly and eccentricity of the orbit
%21/10/2018 - IVAN BARCELONA MORENO

tol_e = 1e-8; % requried tolerance for eccentricity

% -  Hyperbolic
if e>(1+tol_e) %For hyperbolic orbits
    E = 2* atanh( sqrt((e-1)/(e+1)) * tan(Theta/2) );
    if E<0
        E = 2*pi + E;
    end
    M = e*sinh(E) - E;

% - Parabolic 
elseif abs(e-1)<=tol_e % Parabolic orbits
    M = 1/2*tan(Theta/2) + 1/6*(tan(Theta/2))^3;

% - Elipthical
elseif e>(tol_e) && e<(1-tol_e) %For elipthical orbits 
    E = 2* atan( sqrt((1-e)/(1+e)) * tan(Theta/2) );
    if E<0
        E = 2*pi + E;
    end
    M = E - e*sin(E);

% - Circular
else % For circular orbits
    M = Theta;
end


end