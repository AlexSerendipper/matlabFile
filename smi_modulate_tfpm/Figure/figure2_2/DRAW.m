clc;
clear all;
close all;
[x,y,z,m] = DRAW_API_DATE("12");
DRWA_API_TEMPLATE(x,y,z,m(1,:),m(2,:));