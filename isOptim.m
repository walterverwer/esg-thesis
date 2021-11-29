function [position,isterminal,direction] = isOptim(z,y)
  position      = y(2); % The value that we want to be zero
  isterminal    = 0;  % Halt integration 
  direction     = 0;   % The zero can be approached from either direction
end

