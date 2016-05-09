function [P_r, roh_r,T_r, A_r]= static_cond(M, gamma)

P_r = (1 + (gamma-1)*M^2 /2)^( (gamma)/(gamma-1) ) ;
roh_r =(1 + (gamma-1)*M^2 /2)^( 1/(gamma-1) ) ;
T_r = 1 + (gamma-1)*M^2 /2;
A_r = sqrt(1/M^2*( 2/(gamma+1) * (1+(gamma-1)*M^2/2))^((gamma+1)/(gamma-1)));

end