clc;
clear all;
close all;
[x,y,z,m] = DRAW_API_DATE("22");
x = 0.348;
y = 0.3079;
z = [0.06233357972933,0.696375359406102,0];
DRWA_API_TEMPLATE(x,y,z,m(1,:),m(2,:),m(3,:),m(4,:));