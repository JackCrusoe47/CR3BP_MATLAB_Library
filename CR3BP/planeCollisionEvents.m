function [position,isterminal,direction] = planeCollisionEvents(t,y,plane)

switch plane
    case 13
        position = y(2); % The value that we want to be zero
    case 12
        position = y(3); % The value that we want to be zero
    case 23
        position = y(1); % The value that we want to be zero
end

isterminal = 1;  % Halt integration 
direction = 0;   % The zero can be approached from either direction

end