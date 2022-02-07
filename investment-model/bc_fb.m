function res = bc_fb(ya, yb, p_B, p_G)
% Write boundary conditions in the form g(ya, yb) = 0, where ya denotes
% V(0) and yb denotes V(1).
res = [ya(1)- p_B
       yb(1)- p_G];
end

