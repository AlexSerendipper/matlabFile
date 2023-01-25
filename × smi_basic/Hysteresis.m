 % figure out the link of phi0 and phiF under 1<C<4.6
 f = @(phiF,C,a) phiF+C.*sin(phiF+atan(a));  %solve the phi0；
 C = 3; 
 a = 3;  %a means the alpha
 phiF=0:0.1:20;
 phi0=f(phiF,C,a);
 figure(1);
 plot(phi0,phiF);
 xlabel('φ0');
 ylabel("φF");
 hold;
 x1 = [6.59415,6.59415];
 y1 = [9.4,5.4,];
 x2 = [9.76892,9.76892];
 y2 = [7,10.9];
 x3 = [12.8727,12.8727];
 y3 = [15.6,11.7];
 x4 = [16.0491,16.0491];
 y4 = [17.2,13.3];
 plot(x1,y1,':'); 
 plot(x2,y2,':');
 plot(x3,y3,':');
 plot(x4,y4,':');
 text(6.59415,9.4,"a");
 text(6.59415,5.4,"c");
 text(9.76892,10.9,"d");
 text(9.76892,7,"b");
 text(12.8727,15.6,"A");
 text(12.8727,11.7,"C");
 text(16.0491,17.2,"B");
 text(16.0491,13.3,"D");
 
 % figure out the link of phi0 and G(phiF) under 1<C<4.6
 figure(2);
 Gf = cos(phiF);
 phi0=f(phiF,C,a);
 plot(phi0,Gf);
 xlabel('φ0');
 ylabel("G(φF)");
 hold;
 x1 = [6.59415,6.59415];
 y1 = [0.70867,-0.99222];
 x2 = [9.77035,9.77035];
 y2 = [0.815725,-0.0954289];
 x3 = [12.8727,12.8727];
 y3 = [0.65607,-0.994178];
 x4 = [16.0553,16.0553];
 y4 = [0.805884,-0.0786782];
 plot(x1,y1,':'); 
 plot(x2,y2,':');
 plot(x3,y3,':');
 plot(x4,y4,':');
 text(6.59415,0.70867,"a");
 text(6.59415,-0.99222,"c");
 text(9.77035,0.815725,"b");
 text(9.77035,-0.0954289,"d");
 text(12.8727,0.65607,"A");
 text(12.8727,-0.994178,"C");
 text(16.0553,0.805884,"B");
 text(16.0553,-0.0786782,"D");
 