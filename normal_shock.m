function [M_2, p_s, T_s, po2_o1, po1_1, po2_1]=  normal_shock(M, gamma)

M_2 = sqrt( ((1 + (gamma-1)/2*M^2 )/( (gamma*M^2 -(gamma-1)/2))));
p_s = 1 + (2*gamma/(gamma+1))*(M^2-1);
T_s = p_s*( (2 + (gamma -1)*M^2)/( (gamma+1)*M^2));
ss_1 = (((gamma + 1)^2*M^2)/( 4*gamma*M^2 - 2*(gamma -1)) ) ^(gamma/(gamma -1));
ss_2 = (1-gamma+2*gamma*M^2)/(gamma+1);
po2_1 = ss_1*ss_2;  % P_02/P_1
[po1_1, roh, T, A] = static_cond(M, gamma); % P_01/P_1
po2_o1 = po2_1/po1_1; % P_02/P_01
end