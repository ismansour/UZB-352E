%% The second problem
% 

clc
clear all
close all
format long
p_inf = 10332.27555; % mmH2O
Me = [1.39224743 1.778409538 1.925918898];


gamma = 1.4;
p_o = p_inf; %atm
alpha = 0; %deg
angle = 6; %deg

%%      
%       For region 2:
for n=1:length(Me)
 
M_e = Me(n);
M_1 = Me(n);
theta = -alpha + angle;
[po1_p1] = static_cond(M_e, gamma);
[M_2, beta, po2_o1, p2_p1] = oblique_shock(M_e, theta);
A_1 = sqrt((gamma+1)/(gamma-1));
A_2 = sqrt(M_2^2 -1);
v_2 = ((A_1)*atan( (A_2/A_1) ) - atan((A_2)) )*180/pi;
[po2_p2] = static_cond(M_2, gamma);
fprintf('\n For region 2:') 
fprintf('\n theta     = %g', theta)
fprintf('\n beta      = %g', beta)
fprintf('\n P_o1/P_1  = %g', po1_p1)
fprintf('\n P_o2/P_o1 = %g', po2_o1)
fprintf('\n M_2       = %g', M_2)
fprintf('\n P_2/P_1   = %g', p2_p1)
fprintf('\n v_2       = %g', v_2)
fprintf('\n P_o4/P_4  = %g', po2_p2)
%%      
%       For Region 3:
v_3 = v_2 + 2*angle;
syms M_3 
A_1 = sqrt((gamma+1)/(gamma-1));
A_2 = sqrt(M_3^2 -1);
equ = v_3 == ((A_1)*atan( (A_2/A_1) ) - atan((A_2)) )*180/pi;

sol_M_3 = vpasolve(equ, M_3,M_2);
M_3 = double(sol_M_3);
[po3_p3] = static_cond(M_3, gamma);
fprintf('\n For region 3:') 
fprintf('\n theta     = %g', 2*angle)
fprintf('\n v_2       = %g', v_3)
fprintf('\n M_3       = %g', M_3)
fprintf('\n P_o1/P_1  = %g', po3_p3)

%%
%       For region 4:
% 

theta = alpha + angle;

[M_4, beta, po4_o1, p4_p1] = oblique_shock(M_e, theta);
A_1 = sqrt((gamma+1)/(gamma-1));
A_2 = sqrt(M_4^2 -1);
v_4 = ((A_1)*atan( (A_2/A_1) ) - atan((A_2)) )*180/pi;
[po4_p4] = static_cond(M_4, gamma);
fprintf('\n For region 4:') 
fprintf('\n theta     = %g', theta)
fprintf('\n beta      = %g', beta)
fprintf('\n P_o4/P_o1 = %g', po4_o1)
fprintf('\n M_4       = %g', M_4)
fprintf('\n P_4/P_1   = %g', p4_p1)
fprintf('\n v_5       = %g', v_4)
fprintf('\n P_o4/P_4  = %g', po4_p4)
%%
%       For region 5:

v_5 = v_4+ 2*angle;
syms M_5 
A_1 = sqrt((gamma+1)/(gamma-1));
A_2 = sqrt(M_5^2 -1);
equ = v_5 == ((A_1)*atan( (A_2/A_1) ) - atan((A_2)) )*180/pi;

sol_M_5 = vpasolve(equ, M_5,M_4);
M_5 = double(sol_M_5);
[po5_p5] = static_cond(M_5, gamma);
fprintf('\n For region 5:') 
fprintf('\n theta     = %g', 2*angle)
fprintf('\n v_5       = %g', v_5)
fprintf('\n M_5       = %g', M_5)
fprintf('\n P_o5/P_5  = %g', po5_p5)
%%
%       Pressure ratios
p2_p1 = p2_p1; % p_2/p_1 = p_4/p_1
p3_p1 = 1/po3_p3 * po2_o1*po1_p1; % p_3/p_1 = p_1/p_1 p_3/p_o3 p_o3/p_o1 p_o2/p_2 
p4_p1 = p4_p1; % p_4/p_1
p5_p1 = 1/po5_p5 * po4_o1*po1_p1; % p_5/p_1 = p_5/p_o5 p_o5/p_o4 p_o4/p_o1
% p_o1/p_1
fprintf('\n Pressure ratios:') 
fprintf('\n P_2/P_1   = %g', p2_p1)
fprintf('\n P_3/P_1   = %g', p3_p1)
fprintf('\n P_4/P_1   = %g', p4_p1)
fprintf('\n P_5/P_1   = %g', p5_p1)
%%
%       The Drag and Lift Coefficients
lc = 1/(2*cosd(angle));
c_l = 2*lc*( (p4_p1-p3_p1)*cosd(alpha+angle) + (p5_p1-p2_p1)*cosd(alpha-angle) )/(gamma*M_e^2);
c_d = 2*lc*( (p4_p1-p3_p1)*sind(alpha+angle) + (p5_p1-p2_p1)*sind(alpha-angle) )/(gamma*M_e^2);
fprintf('\n the drag and lift coefficients:') 
fprintf('\n C_L   = %g', c_l)
fprintf('\n C_D   = %g', c_d)

end