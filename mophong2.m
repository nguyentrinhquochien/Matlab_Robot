%tich phan thich nghi
close all
clear all
%thong so bo dieu khien
Q1=25;Q2=25;Q3=25;Q=[Q1 0 0;0 Q2 0;0 0 Q3];

%Q1=10;Q2=10;Q3=10;Q=[Q1 0 0;0 Q2 0;0 0 Q3];nho khong tot

P1=1.5;P2=1.5;P3=1.5;P=[P1 0 0;0 P2 0;0 0 P3];

%Teta1=0.3; Teta2=0.3; Teta3=0.3;  %cung lon

Teta1=0.1; Teta2=0.1; Teta3=0.1;  %cung lon

% load D:\mophmoi\ref_ISMC xr yr phir dxr dyr d2xr d2yr dwr wr vr n
%%%%%%%%% qu? ??o 1
% dt=0.01;   
% Vr=0.05;
% Rr=6240*Vr*dt*(4/(2*pi)); 
% Lr=4000*Vr*dt;             
% 
% x1=0.25;					y1=0.4;
% x2=x1+Lr;   				y2=y1;      
% x3=x2+Rr;					y3=y2+Rr;   
% x4=x3-Rr;                   y4=y3+Rr;   
% x5=x4-Rr;	                y5=y4+Rr;   
% x6=x5+Rr;                   y6=y5+Rr;   
% 
% xr1(1)=x1;
% yr1(1)=y1;
% phir1(1)=0;
% wr1(1)=0;
% 
% n1=round(Lr/(Vr*dt));
% for i=2:n1
%    xr1(i)=xr1(i-1)+Vr*dt;
%    yr1(i)=yr1(i-1); % ngang 1
%    phir1(i)=0;
%    wr1(i)=0;
% end;
% 
% xr2(1)=x2;
% yr2(1)=y2;
% phir2(1)=0;
% wr2(1)=Vr/Rr;
% 
% n2=round(pi*Rr/(2*Vr*dt));   
% for i=2:n2
%    if i==2 i23=0; end;
%    xr2(i)=x2+Rr*sin(Vr*i23*dt/Rr);
%    yr2(i)=y2+Rr*(1-cos(Vr*i23*dt/Rr));
%    phir2(i)=Vr*i23*dt/Rr;    
%    wr2(i)=Vr/Rr;
%    i23=i23+1;
% end;
% 
% xr3(1)=x3;
% yr3(1)=y3;
% phir3(1)= pi/2;
% wr3(1)= Vr/Rr;
% 
% n3=round(pi*Rr/(2*Vr*dt));
% for i=2:n3
%    if i==2 i34=0; end;
%    xr3(i)=x3-Rr*(1-cos(Vr*i34*dt/Rr));
%    yr3(i)=y3+Rr*sin(Vr*i34*dt/Rr);
%    phir3(i)= pi/2 + Vr*i34*dt/Rr;    
%    wr3(i)=Vr/Rr;
%    i34=i34+1;
% end;
% 
% xr4(1)=x4;
% yr4(1)=y4;
% phir4(1)= pi;
% wr4(1)= -Vr/Rr;
% 
% n4=round(pi*Rr/(2*Vr*dt));
% for i=2:n4
%    if i==2 i45=0; end;
%    xr4(i)=x4-Rr*sin(Vr*i45*dt/Rr);
%    yr4(i)=y4+Rr*(1-cos(Vr*i45*dt/Rr));
%    phir4(i)= pi-Vr*i45*dt/Rr;
%    wr4(i)= -Vr/Rr;   
%    i45=i45+1;   
% end;
% 
% xr5(1)=x5;
% yr5(1)=y5;
% phir5(1)=pi/2;
% wr5(1)=-Vr/Rr;
% %n5=400;
% n5=round(pi*Rr/(2*Vr*dt));
% for i=2:n5
%    if i==2 i56=0; end;
%    xr5(i)=x5+Rr*(1-cos(Vr*i56*dt/Rr));
%    yr5(i)=y5 + Rr*sin(Vr*i56*dt/Rr);     %ngang 4
%    phir5(i)=pi/2-Vr*i56*dt/Rr;
%    wr5(i)=-Vr/Rr;
%    i56=i56+1;
% end;
% 
% 
% xr=[xr1,xr2,xr3,xr4,xr5];
% yr=[yr1,yr2,yr3,yr4,yr5];
% 
% phir=[phir1,phir2,phir3,phir4,phir5];
% wr=[wr1,wr2,wr3,wr4,wr5];
% 
% n=n1+n2+n3+n4+n5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%so 8
dt=0.01;   
Vr=0.05;

Rr=6240*Vr*dt*(4/(2*pi)); 
Lr=4000*Vr*dt;             

x1=0.25;					y1=0.4;
x2=x1+Rr;   				y2=y1+Rr;     
x3=x2-Rr;					y3=y2+Rr;   
x4=x3-Rr;					y4=y3+Rr;   
x5=x4+Rr;                    y5=y4+Rr;
x6=x5+Rr;                    y6=y5-Rr;
x7=x6-Rr;                    y7=y6-Rr;
x8=x7-Rr;                    y8=y7-Rr;
x9=x8+Rr;                    y9=y8-Rr;


xr1(1)=x1;
yr1(1)=y1;
phir1(1)=0;
wr1(1)=0;

n1=round(pi*Rr/(2*Vr*dt));
for i=2:n1
  if i==2 i12=0; end;
   xr1(i)=x1+Rr*sin(Vr*i12*dt/Rr);
   yr1(i)=y1+Rr*(1-cos(Vr*i12*dt/Rr));
   phir1(i)= Vr*i12*dt/Rr;
   wr1(i)= Vr/Rr;
   i12=i12+1;
end;

xr2(1)=x2;
yr2(1)=y2;
phir2(1)=pi/2;
wr2(1)=Vr/Rr;

n2=round(pi*Rr/(2*Vr*dt));   
for i=2:n2
   if i==2 i23=0; end;
   xr2(i)=x2-Rr*(1-cos(Vr*i23*dt/Rr));
   yr2(i)=y2+Rr*sin(Vr*i23*dt/Rr);
   phir2(i)=pi/2 + Vr*i23*dt/Rr;
   wr2(i)=Vr/Rr;
   i23=i23+1;
end;

xr3(1)=x3;
yr3(1)=y3;
phir3(1)=pi;
wr3(1)=Vr/Rr;
n3=round(pi*Rr/(2*Vr*dt));
for i=2:n3
   if i==2 i34=0; end;
   xr3(i)=x3-Rr*sin(Vr*i34*dt/Rr);
   yr3(i)=y3+Rr*(1-cos(Vr*i34*dt/Rr));
   phir3(i)=pi - Vr*i34*dt/Rr;
   wr3(i)=Vr/Rr;
   i34=i34+1;
end;

xr4(1)=x4;
yr4(1)=y4;
phir4(1)=pi/2;
wr4(1)=Vr/Rr;

n4=round(pi*Rr/(2*Vr*dt));
for i=2:n4
   if i==2 i45=0; end;
   xr4(i)=x4+Rr*(1-cos(Vr*i45*dt/Rr));
   yr4(i)=y4+Rr*sin(Vr*i45*dt/Rr);
   phir4(i)=pi/2-Vr*i45*dt/Rr;
   wr4(i)=Vr/Rr;   
   i45=i45+1;   
end;

xr5(1)=x5;
yr5(1)=y5;
phir5(1)=0;
wr5(1)=Vr/Rr;

n5=round(pi*Rr/(2*Vr*dt));
for i=2:n5
   if i==2 i56=0; end;
   xr5(i)=x5+Rr*sin(Vr*i56*dt/Rr);
   yr5(i)=y5-Rr*(1-cos(Vr*i56*dt/Rr));
   phir5(i)=-Vr*i56*dt/Rr;
   wr5(i)=Vr/Rr;   
   i56=i56+1;   
end;

xr6(1)=x6;
yr6(1)=y6;
phir6(1)=-pi/2;
wr6(1)=Vr/Rr;

n6=round(pi*Rr/(2*Vr*dt));
for i=2:n6
   if i==2 i67=0; end;
   xr6(i)=x6-Rr*(1-cos(Vr*i67*dt/Rr));
   yr6(i)=y6-Rr*sin(Vr*i67*dt/Rr);
   phir6(i)=-pi/2-Vr*i67*dt/Rr;
   wr6(i)=Vr/Rr;   
   i67=i67+1;   
end;

xr7(1)=x7;
yr7(1)=y7;
phir7(1)= - pi;
wr7(1)=Vr/Rr;

n7=round(pi*Rr/(2*Vr*dt));
for i=2:n7
   if i==2 i78=0; end;
   xr7(i)=x7-Rr*sin(Vr*i78*dt/Rr);
   yr7(i)=y7-Rr*(1-cos(Vr*i78*dt/Rr));
   phir7(i)=-pi+Vr*i78*dt/Rr;
   wr7(i)=Vr/Rr;   
   i78=i78+1;   
end;

xr8(1)=x8;
yr8(1)=y8;
phir8(1)=-pi/2;
wr8(1)=Vr/Rr;

n8=round(pi*Rr/(2*Vr*dt));
for i=2:n8
   if i==2 i89=0; end;
   xr8(i)=x8+Rr*(1-cos(Vr*i89*dt/Rr));
   yr8(i)=y8-Rr*sin(Vr*i89*dt/Rr);
   phir8(i)=-pi/2+Vr*i89*dt/Rr;
   wr8(i)=Vr/Rr;   
   i89=i89+1;   
end;
xr=[xr1,xr2,xr3,xr4,xr5,xr6,xr7,xr8];
yr=[yr1,yr2,yr3,yr4,yr5,yr6,yr7,yr8];

phir=[phir1,phir2,phir3,phir4,phir5,phir6,phir7,phir8];
wr=[wr1,wr2,wr3,wr4,wr5,wr6,wr7,wr8];

n=n1+n2+n3+n4+n5+n6+n7+n8;

%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%% qu? ??o tròn
% dt=0.01;   
% Vr=0.05;
% 
% Rr=6240*Vr*dt*(4/(2*pi)); 
% Lr=4000*Vr*dt;             
% 
% x1=0.25;					y1=0.4;
% x2=x1+Rr;   				y2=y1+Rr;      % ngang 1
% x3=x2-Rr;					y3=y2+Rr;   % cong trai 1
% x4=x3-Rr;					y4=y3-Rr;   % dung thang
% x5=x4+Rr;                    y5=y4-Rr;   % cong phai 2
% 
% xr1(1)=x1;
% yr1(1)=y1;
% phir1(1)=0;
% wr1(1)=0;
% 
% n1=round(pi*Rr/(2*Vr*dt));
% for i=2:n1
%   if i==2 i12=0; end;
%    xr1(i)=x1+Rr*sin(Vr*i12*dt/Rr);
%    yr1(i)=y1+Rr*(1-cos(Vr*i12*dt/Rr));
%    phir1(i)= Vr*i12*dt/Rr;
%    wr1(i)= Vr/Rr;
%    i12=i12+1;
% end;
% 
% xr2(1)=x2;
% yr2(1)=y2;
% phir2(1)=pi/2;
% wr2(1)=Vr/Rr;
% 
% n2=round(pi*Rr/(2*Vr*dt));   
% for i=2:n2
%    if i==2 i23=0; end;
%    xr2(i)=x2-Rr*(1-cos(Vr*i23*dt/Rr));
%    yr2(i)=y2+Rr*sin(Vr*i23*dt/Rr);
%    phir2(i)=pi/2 + Vr*i23*dt/Rr;
%    wr2(i)=Vr/Rr;
%    i23=i23+1;
% end;
% 
% xr3(1)=x3;
% yr3(1)=y3;
% phir3(1)=pi;
% wr3(1)=Vr/Rr;
% n3=round(pi*Rr/(2*Vr*dt));
% for i=2:n3
%    if i==2 i34=0; end;
%    xr3(i)=x3-Rr*sin(Vr*i34*dt/Rr);
%    yr3(i)=y3-Rr*(1-cos(Vr*i34*dt/Rr));
%    phir3(i)=pi + Vr*i34*dt/Rr;
%    wr3(i)=Vr/Rr;
%    i34=i34+1;
% end;
% 
% xr4(1)=x4;
% yr4(1)=y4;
% phir4(1)=3*pi/2;
% wr4(1)=Vr/Rr;
% 
% n4=round(pi*Rr/(2*Vr*dt));
% for i=2:n4
%    if i==2 i45=0; end;
%    xr4(i)=x4+Rr*(1-cos(Vr*i45*dt/Rr));
%    yr4(i)=y4-Rr*sin(Vr*i45*dt/Rr);
%    phir4(i)=3*pi/2+Vr*i45*dt/Rr;
%    wr4(i)=Vr/Rr;   
%    i45=i45+1;   
% end;
% xr=[xr1,xr2,xr3,xr4];
% yr=[yr1,yr2,yr3,yr4];
% 
% phir=[phir1,phir2,phir3,phir4];
% wr=[wr1,wr2,wr3,wr4];
% 
% n=n1+n2+n3+n4;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dxr(1)=0;dyr(1)=0;dphir(1)=0;dwr(1)=0;
for i=2:n
   dxr(i)=(xr(i)-xr(i-1))/dt;
   dyr(i)=(yr(i)-yr(i-1))/dt;
   dphir(i)=(phir(i)-phir(i-1))/dt;
   dwr(i)=(wr(i)-wr(i-1))/dt;
end;

dxr(1)=0;  %bac mot
dyr(1)=0;  %bac mot
d2xr(1)=0; %bac hai
d2yr(1)=0; %bac hai
dwr(1)=0;  %bac hai
for i=2:n
   %dxr(i)=(xr(i)-xr(i-1))/dt;
   d2xr(i)=(dxr(i)-dxr(i-1))/dt;
   %dyr(i)=(yr(i)-yr(i-1))/dt;
   d2yr(i)=(dyr(i)-dyr(i-1))/dt;
   dwr(i)=(wr(i)-wr(i-1))/dt;
end;

vr=Vr;

k1=23; k2=23; k3=23;K=diag([k1,k2,k3]);  
%k1=10; k2=10; k3=10;K=diag([k1,k2,k3]);

kv1=15; kv2=15; kv3=15;Kv=diag([kv1,kv2,kv3]);
        
dt=0.01;

%****************************************
%Main Program
%****************************************
%thong so cua OMP
L=0.18;	    	 
r=0.04;          
m=4.5;		    
J=0.12;		    
anpha=0.2;
peta=1.5;

%cac ma tran
M=[m 0 0;0 m 0;0 0 J];
V=[1.5*peta 0 0;0 1.5*peta 0;0 0 3*peta*L*L];
%vi tri OMP luc dau
x(1)=1;
y(1)=1;
phi(1)=101*pi/180; 

%van toc OMR luc dau
dx(1)=0; 
dy(1)=0;
dphi(1)=0;
v(1)=0;
w(1)=0;
dq=[dx(1);dy(1);w(1)];
dqr=[dxr(1);dyr(1);wr(1)];
d2qr=[d2xr(1);d2yr(1);dwr(1)];

%bien do cua nhieu
fM1=2;fM2=2;fM3=2;fA1=1.5;fA2=1.5;fA3=1.5;

%gia toc OMP luc dau
d2x(1)=0; 
d2y(1)=0; 
dw(1)=0; 

%matran H^-1 luc dau
H_inv=[-sin(phi(1)) cos(phi(1)) L
-sin((pi/3)-phi(1)) -cos((pi/3)-phi(1)) L
sin((pi/3)+phi(1)) -cos((pi/3)+phi(1)) L];

H_inv_dot=[-w(1)*cos(phi(1)) -w(1)*sin(phi(1)) 0
w(1)*cos((pi/3)-phi(1)) -w(1)*sin((pi/3)-phi(1)) 0
w(1)*cos((pi/3)+phi(1)) w(1)*sin((pi/3)+phi(1)) 0];

H=inv(H_inv);

%tinh z (3.9)
z=(1/r)*H_inv*dq;

%cac sai so luc dau
e1(1)=xr(1)-x(1);
e2(1)=yr(1)-y(1);
e3(1)=phir(1)-phi(1);
e=[e1(1);e2(1);e3(1)];

%dao ham cac sai so luc dau
%de1(1)=0; %de2(1)=0; %de3(1)=0; 
de1(1)=dxr(1)-dx(1);
de2(1)=dyr(1)-dy(1);
de3(1)=wr(1)-w(1);
de=[de1(1);de2(1);de3(1)];

%tin hieu van toc 
zd=(1/r)*H_inv*(K*e+dqr);
wd1(1)=zd(1);
wd2(1)=zd(2);
wd3(1)=zd(3);

dwd1(1)=0;
dwd2(1)=0;
dwd3(1)=0;
dzd=[dwd1(1);dwd2(1);dwd3(1)];

%sai so van toc
ev=z-zd;
ev1(1)=ev(1);
ev2(1)=ev(2);
ev3(1)=ev(3);

%mat phang truot luc dau
iev1(1)=0;
iev2(1)=0;
iev3(1)=0;
iev=[iev1(1);iev2(1);iev3(1)];

S=ev+Kv*iev;
S1(1)=S(1);      
S2(1)=S(2);
S3(1)=S(3);

%signS=[sign(S1(1));sign(S2(1));sign(S3(1))];
satS=[sign(S1(1)/Teta1);sign(S2(1)/Teta2);sign(S3(1)/Teta3)];

%bom nhieu ngau nhien luc dau
f1d(1)=-fM1*sin(phi(1))-fM2*sin(pi/3-phi(1))+fM3*sin(pi/3+phi(1))+fA1*cos(phi(1))+fA2*cos(2*pi/3+phi(1))+fA3*cos(4*pi/3+phi(1));
f2d(1)=fM1*cos(phi(1))-fM2*cos(pi/3-phi(1))-fM3*cos(pi/3+phi(1))+fA1*sin(phi(1))+fA2*sin(2*pi/3+phi(1))+fA3*sin(4*pi/3+phi(1));
f3d(1)=L*(fM1+fM2-fM3);
fd=[f1d(1);f2d(1);f3d(1)];

%tin hieu dieu khien luc dau
M_bar=(1/anpha)*H'*M;
V_bar=(1/anpha)*H'*V;
ud=(1/anpha)*H'*fd;
u=M_bar*H;
u= u*(-(H_inv_dot+Kv*H_inv+H_inv*V_bar)*dq+(H_inv_dot*K+Kv*H_inv*K)*e+H_inv*K*de+(H_inv_dot+Kv*H_inv)*dqr+H_inv*d2qr+H_inv*ud)-r*Q*S-r*P*satS  
u1(1)=u(1); 
u2(1)=u(2);
u3(1)=u(3);

%van toc banh xe luc dau 
w1(1)=z(1);w2(1)=z(2); w3(1)=z(3);
%w1(1)=z(1)*(180/pi);%w2(1)=z(2)*(180/pi); %w3(1)=z(3)*(180/pi); 

%Main Loop 
for i=2:n
t(i)=(i-1)*dt;

%tinh gia toc cua xe
d2q=M_bar^-1*(u-ud-V_bar*dq);
d2x(i)=d2q(1); 
d2y(i)=d2q(2); 
dw(i)=d2q(3); 

%tinh van toc cua xe
dx(i)=dx(i-1)+ d2x(i)*dt; 
dy(i)=dy(i-1)+ d2y(i)*dt;
w(i)=w(i-1)+ dw(i)*dt;
dq=[dx(i);dy(i);w(i)];
dqr=[dxr(i);dyr(i);wr(i)];

%tinh vi tri cua xe
x(i)=x(i-1)+ dx(i)*dt; 
y(i)=y(i-1)+ dy(i)*dt;
phi(i)=phi(i-1)+ w(i)*dt;

%cac sai so 
e1(i)=xr(i)-x(i);
e2(i)=yr(i)-y(i);
e3(i)=phir(i)-phi(i);
e=[e1(i);e2(i);e3(i)];

%dao ham cac sai so
%de1(i)=(e1(i)-e1(i-1))/dt; 
%de2(i)=(e2(i)-e2(i-1))/dt; 
%de3(i)=(e3(i)-e3(i-1))/dt; 
de1(i)=dxr(i)-dx(i);
de2(i)=dyr(i)-dy(i);
de3(i)=wr(i)-w(i);

%matran H luc dau
H_inv=[-sin(phi(i)) cos(phi(i)) L
-sin((pi/3)-phi(i)) -cos((pi/3)-phi(i)) L
sin((pi/3)+phi(i)) -cos((pi/3)+phi(i)) L];

H_inv_dot=[-w(i)*cos(phi(i)) -w(i)*sin(phi(i)) 0
w(i)*cos((pi/3)-phi(i)) -w(i)*sin((pi/3)-phi(i)) 0
w(i)*cos((pi/3)+phi(i)) w(i)*sin((pi/3)+phi(i)) 0];
H=inv(H_inv);

%tinh z
z=(1/r)*H_inv*dq;

%zd 
zd=(1/r)*H_inv*(K*e+dqr);
wd1(i)=zd(1);
wd2(i)=zd(2);
wd3(i)=zd(3);

dwd1(i)=(wd1(i)-wd1(i-1))/dt;
dwd2(i)=(wd2(i)-wd2(i-1))/dt;
dwd3(i)=(wd3(i)-wd3(i-1))/dt;
dzd=[dwd1(i);dwd2(i);dwd3(i)];

%sai so van toc
ev=z-zd;
ev1(i)=ev(1);
ev2(i)=ev(2);
ev3(i)=ev(3);

%mat phang truot luc dau
iev1(i)=iev1(i-1)+ev1(i)*dt;
iev2(i)=iev2(i-1)+ev2(i)*dt;
iev3(i)=iev3(i-1)+ev3(i)*dt;
iev=[iev1(i);iev2(i);iev3(i)];

S=ev+Kv*iev;
S1(i)=S(1);              
S2(i)=S(2);
S3(i)=S(3);

%signS=[sign(S1(1));sign(S2(1));sign(S3(1))];
satS=[sign(S1(i)/Teta1);sign(S2(i)/Teta2);sign(S3(i)/Teta3)];

%bom nhieu ngau nhien vao
f1d(i)=-fM1*sin(phi(i))-fM2*sin(pi/3-phi(i))+fM3*sin(pi/3+phi(i))+fA1*cos(phi(i))+fA2*cos(2*pi/3+phi(i))+fA3*cos(4*pi/3+phi(i));
f2d(i)=fM1*cos(phi(i))-fM2*cos(pi/3-phi(i))-fM3*cos(pi/3+phi(i))+fA1*sin(phi(i))+fA2*sin(2*pi/3+phi(i))+fA3*sin(4*pi/3+phi(i));
f3d(i)=L*(fM1+fM2+fM3);
fd=[f1d(i);f2d(i);f3d(i)];

%tin hieu dieu khien
M_bar=(1/anpha)*H'*M;
V_bar=(1/anpha)*H'*V;
ud=(1/anpha)*H'*fd;
u=M_bar*H;

%u=u*(-(H_inv_dot+C*H_inv*V_bar)*dq+(H_inv_dot*K+C*H_inv*K)*e+H_inv*K*de+(H_inv_dot+C*H_inv*K)*dqr+H_inv*d2qr+H_inv*ud)-Q*S-P*satS;  
u=u*(-(H_inv_dot+Kv*H_inv+H_inv*V_bar)*dq+(H_inv_dot*K+Kv*H_inv*K)*e+H_inv*K*de+(H_inv_dot+Kv*H_inv)*dqr+H_inv*d2qr+H_inv*ud)-Q*S-P*satS;  
u1(1)=u(1); 
u2(1)=u(2);
u3(1)=u(3);

%z 
w1(i)=z(1);   
w2(i)=z(2); 
w3(i)=z(3);

%V 
v(i)=sqrt(dx(i)*dx(i)+dy(i)*dy(i));

end;

%ket qua mo phong
%**************************
%qu? ??o 1 
%     plot(x1,y1,'k*',x2,y2,'k*',x3,y3,'k*',x4,y4,'k*',x5,y5,'k*',x6,y6,'k*',[xr1,xr2,xr3,xr4,xr5],[yr1,yr2,yr3,yr4,yr5],'r--',xr,yr,'b:','LineWidth', 3);
%     title('quy dao tham chieu va quy dao thuc');
%     xlabel('truc X (m)');ylabel('Truc Y (m)');
%     legend('vitri 1','vitri 2','vitri 3','vitri 4','vitri 5','vitri 6','thucte','tham chieu',1);
%     axis equal;
%*********************
%qu? ??o s? 8 
plot(x1,y1,'k*',x2,y2,'k*',x3,y3,'k*',x4,y4,'k*',x5,y5,'k*',x6,y6,'k*',x7,y7,'k*',x8,y8,'k*',x9,y9,'k*',[xr1,xr2,xr3,xr4,xr5,xr6,xr7,xr8],[yr1,yr2,yr3,yr4,yr5,yr6,yr7,yr8],'r--',xr,yr,'b:','LineWidth', 3);
title('quy dao tham chieu va quy dao thuc');
xlabel('truc X (m)');ylabel('Truc Y (m)');
legend('vitri 1','vitri 2','vitri 3','vitri 4','vitri 5','vitri 6','vitri 7','vitri 8','vitri 9','thucte','tham chieu',1);
axis equal;
%*********************
%qu? ??o hình tròn
%     plot(x1,y1,'k*',x2,y2,'k*',x3,y3,'k*',x4,y4,'k*',[xr1,xr2,xr3,xr4],[yr1,yr2,yr3,yr4],'r--',xr,yr,'b:','LineWidth', 2);
%     title('quy dao tham chieu va quy dao thuc');
%     xlabel('truc X (m)');ylabel('Truc Y (m)');
%     legend('vitri 1','vitri 2','vitri 3','vitri 4','vitri 5','thucte','tham chieu',1);
%     axis equal;
%*********************


figure;plot(t,e1*1000,'r-',t,e2*1000,'k-',t,e3*180/pi,'b-');xlabel('Thoi gian (s)');ylabel('Vec to sai so vi tri ep');legend('e1','e2','e3',1);
figure;plot(t,ev1*180/pi,'r-',t,ev2*180/pi,'k-',t,ev3*180/pi,'b-');xlabel('Thoi gian (s)');ylabel('Vec to sai so van toc ev');legend('ev1','ev2','ev3',1);
%figure;plot(t,e1*1000,'k-',t,e2*1000,'k-.',t,e3*180/pi,'k--');xlabel('Time(s)');ylabel('Tracking error vector e1(mm), e2(mm), e3 (deg)');legend('e1','e2','e3',1);
figure;	plot(t,w1,'r-',t,w2,'k-',t,w3,'b-');xlabel('Thoi gian (s)');ylabel('Van toc goc cua cac banh xe [rad/s]');%legend('w1','w2','w3',1);	
%figure;	plot(t,w1,'r-',t,w2,'k-',t,w3,'b-');xlabel('Time(s)');ylabel('Angular velocities of three wheels z(deg/s)');
figure;	plot(xr,yr,'k --',x,y,':','LineWidth', 1);xlabel('Truc X [m]');ylabel('Truc Y [m]');axis equal; legend('SP','CV',3);
%figure;	plot(xr,yr,'k-',x,y,'k--');xlabel('X coordinate (m)');ylabel('Y coordinate (m)');%title('Movement of omnidirectional mobile platform');
figure;	plot(t,v,'r-');xlabel('Thoi gian (s)');ylabel('Van toc dai cua OMR [m/s]');%title('Linear velocity of OMP (m/s)');
figure;	plot(t,w,'b-');xlabel('Thoi gian (s)');ylabel('Van toc goc cua OMR [rad/s]');%title('Angular velocity of OMP (rad/s)');	
figure;	plot(t,S1,'r-',t,S2,'k-',t,S3,'b-');xlabel('Thoi gian(s)');ylabel('vector mat truot S [rad/s]');%legend('S1','S2','S3',1);
figure;	plot(t,S1,'r-');xlabel('Time (s)');ylabel('Sliding surface S1 [m/s]');	
figure;	plot(t,S2,'b-');xlabel('Time (s)');ylabel('Sliding surface S2 [m/s]');
figure;	plot(t,S3,'k-');xlabel('Time (s)');ylabel('Sliding surface S3 [rad/s]');	
figure;	plot(t,S3*180/pi,'k-');xlabel('Time (s)');ylabel('Sliding surface S3 (deg/s)');
%figure;    plot(t,tt1,'r-',t,tt2,'k-',t,tt3,'b-');xlabel('Time(s)');ylabel('Torque input control vector [N.m]');%legend('tt1','tt2','tt3',1);
figure;    plot(t,f1d,'r-',t,f2d,'k-',t,f3d,'b-');xlabel('Time(s)');ylabel('Friction force vector [N.m]');%legend('f1d','f2d','f3d',1);
