clc;
clear all;
close all;
[x,y,z,m] = DRAW_API_DATE("52");
DRWA_API_TEMPLATE(x,y,z,m(1,:),m(2,:),m(3,:),m(4,:),m(5,:),m(6,:),m(7,:),m(8,:),m(9,:),m(10,:));
% [Lt_standard,mv_standard,phi_standard,offset_standard] = SMI_API_STANDARD_AUTO(Lt_reconstruct,200, 100, 125000, 3000)