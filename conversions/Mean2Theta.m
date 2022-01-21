function [theta] = Mean2Theta(M,e)
% Estimates True Anomaly for corresponding Mean Anomaly
% Ref : Fundamentals of Astrodynamics and Applications - Vallado
% Ref : Orbital Mechanics for Engineering Students - H. Curtis

tol_e = 1e-8; % requried tolerance for eccentricity
tol = 1e-10; %required tolerance of result at convergence

% - Elipthical orbits
if e>(tol_e) && e<(1-tol_e) %For elipthical orbits
   %intial estimate of E
   if M<pi
       E = M + e/2;
   elseif M>pi
       E = M - e/2;
   else
       E = M;
   end

   while 1
       f1 = E - e*sin(E) - M;
       f2 = 1 - e*cos(E);
       ratio = f1/f2;
       if abs(ratio)<tol
           break;
       else
           E = E - ratio;
       end
   end
   
   %Calculation of Elipthical Trua anomally
   theta = 2* atan( sqrt((1+e)/(1-e)) * tan(E/2) ); 
   
   
% - Hyperbolic orbits   

elseif e>(1+tol_e) %For hyperbolic orbits
    
    % - Initial Guess for F based on e and M
    if e<1.6
        if ( (-pi<M)&&(M<0) ) || M>pi
            F = M-e;
        else
            F = M+e;
        end
    else
        if e<3.6 && abs(M)>pi
            F = M - sign(M)*e;
        else
            F = M/(e-1);
        end
    end
    
    
    while 1
       f1 = e*sinh(F) - F - M;
       f2 = e*cosh(F) - 1;
       ratio = f1/f2;
       if abs(ratio)<tol
           break;
       else
           F = F - ratio;
       end
    end
 
   %Calculation of Hyperbolic Trua anomally
   theta = 2* atan( sqrt((e+1)/(e-1)) * tanh(F/2) ); 

   
% - Parabolic orbit
elseif abs(e-1)<=tol %For parabolic orbits
    % Calcualtion of true anomally of a parabolic orbit
    theta = 2 * atan( ( 3*M + sqrt((3*M)^2 + 1))^(1/3) - ( 3*M + sqrt((3*M)^2 + 1))^(-1/3)  );
    

% - For circular orbit
else 
    theta = M;
end

% - Converting True anomaly to positive
if theta<0
       theta = 2*pi + theta;
end

end