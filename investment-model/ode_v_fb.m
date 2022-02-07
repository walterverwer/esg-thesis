function dy = ode_v_fb(z,y,r,A_G, A_B ,delta,theta,expr)
%ODE_V_FB Function to solve for the ODE of the first best value of V
%   write the second order ODE as a system of first-order equations. Define
%   y(1) := V(z) and y(2) := V'(z)
% if z<1e-1
%     z=1e-1;
% elseif z>1-1e-1
%     z=1-1e-1;
% end


i_G = 1/theta*( 1- 1/(y(1)-z*y(2)) ); % optimal i_G
i_B = 1/theta*( 1-1/(y(1)+(1-z)*y(2)) ); % optimal i_B

dy = zeros(2,1);  
dy(1) = y(2);

dy(2) = ( expr(z) ) * ( r*y(1) - (1-z)*(A_B-i_B) ...
    - z*(A_G-i_G) - (1-z)*(y(1) - z*y(2))*(i_B - theta/2*i_B^2 - delta ) ...
    - z*(y(1) + (1-z)*y(2))*(i_G - theta/2*i_G^2 - delta) )  ;


% % Taylor series at 0.5:
% if z<0.2
%     dy(2) = ( subs(expr,z) ) * ( r*y(1) - (1-z)*(A_B-i_B) ...
%         - z*(A_G-i_G) - (1-z)*(y(1) - z*y(2))*(i_B - theta/2*i_B^2 - delta ) ...
%         - z*(y(1) + (1-z)*y(2))*(i_G - theta/2*i_G^2 - delta) )  ;
% 
% % Taylor series at 0.5:
% elseif z>0.8
%     dy(2) = ( subs(expr,z) ) * ( r*y(1) - (1-z)*(A_B-i_B) ...
%     - z*(A_G-i_G) - (1-z)*(y(1) - z*y(2))*(i_B - theta/2*i_B^2 - delta ) ...
%     - z*(y(1) + (1-z)*y(2))*(i_G - theta/2*i_G^2 - delta) )  ;
% 
% % Original equation:
% else
%     dy(2) = 2/((sigma_B^2+sigma_G^2)*z^2*(1-z)^2) * ( r*y(1) - (1-z)*(A_B-i_B) ...
%     - z*(A_G-i_G) - (1-z)*(y(1) - z*y(2))*(i_B - theta/2*i_B^2 - delta ) ...
%     - z*(y(1) + (1-z)*y(2))*(i_G - theta/2*i_G^2 - delta) )  ;
% end



% % Taylor series at 0.1, 0.9, and original equation respectively:
% if z<0.2
%     dy(2) = ( 32/(sigma_B^2+sigma_G^2) + 256*(z-0.5)^2/(sigma_B^2+sigma_G^2) + 1536*(z-0.5)^4/(sigma_B^2+sigma_G^2) ) * ( r*y(1) - (1-z)*(A_B-i_B) ...
%         - z*(A_G-i_G) - (1-z)*(y(1) - z*y(2))*(i_B - theta/2*i_B^2 - delta ) ...
%         - z*(y(1) + (1-z)*y(2))*(i_G - theta/2*i_G^2 - delta) )  ;
% 
% % Taylor series at 0.5:
% elseif z>0.8
%     dy(2) = ( 32/(sigma_B^2+sigma_G^2) + 256*(z-0.5)^2/(sigma_B^2+sigma_G^2) + 1536*(z-0.5)^4/(sigma_B^2+sigma_G^2) ) * ( r*y(1) - (1-z)*(A_B-i_B) ...
%     - z*(A_G-i_G) - (1-z)*(y(1) - z*y(2))*(i_B - theta/2*i_B^2 - delta ) ...
%     - z*(y(1) + (1-z)*y(2))*(i_G - theta/2*i_G^2 - delta) )  ;
% 
% % Original equation:
% else
%     dy(2) = 2/((sigma_B^2+sigma_G^2)*z^2*(1-z)^2) * ( r*y(1) - (1-z)*(A_B-i_B) ...
%     - z*(A_G-i_G) - (1-z)*(y(1) - z*y(2))*(i_B - theta/2*i_B^2 - delta ) ...
%     - z*(y(1) + (1-z)*y(2))*(i_G - theta/2*i_G^2 - delta) )  ;
% end
% 
% 
% 
% if z<0.2
%     dy(2) = ( 32/(sigma_B^2+sigma_G^2) + 256*(z-0.5)^2/(sigma_B^2+sigma_G^2) + 1536*(z-0.5)^4/(sigma_B^2+sigma_G^2) ) * ( r*y(1) - (1-z)*(A_B-i_B) ...
%         - z*(A_G-i_G) - (1-z)*(y(1) - z*y(2))*(i_B - theta/2*i_B^2 - delta ) ...
%         - z*(y(1) + (1-z)*y(2))*(i_G - theta/2*i_G^2 - delta) )  ;
% 
% % Taylor series at 0.5:
% elseif z>0.8
%     dy(2) = ( 32/(sigma_B^2+sigma_G^2) + 256*(z-0.5)^2/(sigma_B^2+sigma_G^2) + 1536*(z-0.5)^4/(sigma_B^2+sigma_G^2) ) * ( r*y(1) - (1-z)*(A_B-i_B) ...
%     - z*(A_G-i_G) - (1-z)*(y(1) - z*y(2))*(i_B - theta/2*i_B^2 - delta ) ...
%     - z*(y(1) + (1-z)*y(2))*(i_G - theta/2*i_G^2 - delta) )  ;
% 
% % Original equation:
% else
%     dy(2) = 2/((sigma_B^2+sigma_G^2)*z^2*(1-z)^2) * ( r*y(1) - (1-z)*(A_B-i_B) ...
%     - z*(A_G-i_G) - (1-z)*(y(1) - z*y(2))*(i_B - theta/2*i_B^2 - delta ) ...
%     - z*(y(1) + (1-z)*y(2))*(i_G - theta/2*i_G^2 - delta) )  ;
% end

end

