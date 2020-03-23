%%%%%%%%%%%%%%non-smooth control on structural vibration%%%%%%%%%%%%%%%%%%%
% clear;
close all;
clc;
t=0.02;
%1、initial value
M=3.456*10^5*eye(8);
% M=zeros(8,8);
% for iM=1:1:8
%     for jM=1:1:8
%         if iM<=jM
%             M(iM,jM)=3.456*10^5;
%         else
%             M(iM,jM)=0;
%         end
%     end
% end
% C=10^3*[490 467 410 386 348 298 243 196]*[];

C=[490+467 -467 0 0 0 0 0 0;
    -467 467+410 -410 0 0 0 0 0;
    0 -410 410+386 -386 0 0 0 0;
    0 0 -386 386+348 -348 0 0 0;
    0 0 0 -348 348+298 -298 0 0;
    0 0 0 0 -298 298+243 -243 0;
    0 0 0 0 0 -243 243+196 -196;
    0 0 0 0 0 0 -196 196]*10^3;
alpha1=0.1;
k1=3.4;k2=3.26;k3=2.85;k4=2.69;k5=2.43;k6=2.07;k7=1.69;k8=1.37;
Ke=10^5*alpha1*[k1+k2 -k2 0 0 0 0 0 0;
                -k2 k2+k2 -k3 0 0 0 0 0;
                0 -k3 k3+k4 -k4 0 0 0 0;
                0 0 -k4 k4+k5 -k5 0 0 0;
                0 0 0 -k5 k5+k6 -k6 0 0;
                0 0 0 0 -k6 k6+k7 -k7 0;
                0 0 0 0 0 -k7 k7+k8 -k8;
                0 0 0 0 0 0 -k8 k8];
Dy1=2.4;Dy2=2.3;Dy3=2.2;Dy4=2.1;Dy5=2.0;Dy6=1.9;Dy7=1.7;Dy8=1.5;
KI=10^5*10^(-2)*(1-alpha1)*[k1*Dy1+k2*Dy2 -k2*Dy2 0 0 0 0 0 0;
                    -k2*Dy2 k2*Dy2+k3*Dy3 -k3*Dy3 0 0 0 0 0;
                    0 -k3*Dy3 k3*Dy3+k4*Dy4 -k4*Dy4 0 0 0 0;
                    0 0 -k4*Dy4 k4*Dy4+k5*Dy5 -k5*Dy5 0 0 0;
                    0 0 0 0-k5*Dy5 k5*Dy5+k6*Dy6 -k6*Dy6 0 0;
                    0 0 0 0 -k6*Dy6 k6*Dy6+k7*Dy7 -k7*Dy7 0;
                    0 0 0 0 0 -k7*Dy7 k7*Dy7+k8*Dy8 -k8*Dy8;
                    0 0 0 0 0 0 -k8*Dy8 k8*Dy8];
buchanglong=30; 

H=eye(8);%[0 0 0 0 0 0 0 1]';
I0=[1 1 1 1 1 1 1 1]';%[1 0 0 0 0 0 0 0]';
u=zeros(8,buchanglong/t);

z1=zeros(8,buchanglong/t);%x
z2=zeros(8,buchanglong/t);%v
z3=zeros(8,buchanglong/t);%dx
dv=[0 0 0 0 0 0 0 0]';
%%%%%%%%%%%%%%%%%%%%LQR算法状态%%%%%%%%

u1=zeros(8,buchanglong/t);
z11=zeros(8,buchanglong/t);%x
z21=zeros(8,buchanglong/t);%v
z31=zeros(8,buchanglong/t);%dx
accLQR=zeros(8,buchanglong/t);
dv1=[0 0 0 0 0 0 0 0]';
% z=[z1,z2,z3]';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%no control status%%%%%%%%

u2=zeros(8,buchanglong/t);
z12=zeros(8,buchanglong/t);%x
z22=zeros(8,buchanglong/t);%v
z32=zeros(8,buchanglong/t);%dx
accWukong=zeros(8,buchanglong/t);
dv2=[0 0 0 0 0 0 0 0]';
% z=[z1,z2,z3]';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A1=1.0;beta1=0.5;n1=95;gama1=0.5;%%%%%%%%%%%% non-linear parameter



Dy=[Dy1 Dy2 Dy3 Dy4 Dy5 Dy6 Dy7 Dy8]';


% for j=1:1:8
    for i=2:1:buchanglong/t
%%%%%%%%%%%%%%math model%%%%%%%%%%%%%%%%%%%%%%
ddXg(i)=acc2(i);%sin(t);%%%seismic wave from file a%%%%%

z1(:,i)=z1(:,i-1)+z3(:,i-1)*t;
% z2(:,i)=;

for j=1:1:8
    dx(j)=z3(j,i);
    v(j)=z2(j,i);
    dv(j)=Dy(j)^(-1)*(A1*dx(j)-beta1*abs(dx(j))*abs(v(j))^(n1-1)*v(j)-gama1*dx(j)*abs(v(j))^n1);
    %inline('Dy(j)^(-1)*(A1*dx(j)-beta1*abs(dx(j))*abs(v(j))^(n1-1)*v(j)-gama1*dx(j)*abs(v(j))^n1)');
%     dx(j)=z3(j,i);
%     v(j)=z2(j,i);
% v=[v1,v2,v3,v4,v5,v6,v7,v8];
%     kk1=feval(dv(j),dx(j),v(j));
%     kk2=feval(dv(j),dx(j)+t/2,v(j)+(t*kk1)/2);
%     kk3=feval(dv(j),dx(j)+t/2,v(j)+(t*kk2)/2);
%     kk4=feval(dv1,dx(j)+t,v(j)+t*kk3);
%     dv(j)=(t/6)*(kk1+2*kk2+2*kk3+kk4);
end
z2(:,i)=z2(:,i-1)+dv*t;


z=[z1(:,i-1)' z2(:,i-1)' z3(:,i-1)']';
z3(:,i)=z3(:,i-1)+([-inv(M)*Ke -inv(M)*KI -inv(M)*C]*z-I0*ddXg(i)+inv(M)*H*u(:,i-1))*t;
%%%%%%%%%%%%non smooth control algorithm%%%%%%%%%%%%%%%%%
kNonsmooth1=3;kNonsmooth2=2;a1=0.43;a2=0.6;
fai(:,i)=-kNonsmooth1.*sign(z1(:,i)).*(abs(z1(:,i))).^a1-kNonsmooth2.*sign(z3(:,i)).*(abs(z3(:,i))).^a2;


%%%%%%%%%%%%inverse transform%%%%%%%%%%%%%%%%%%%%
S_mao=[-inv(M)*Ke -inv(M)*KI -inv(M)*C];
u(:,i)=inv(H)*M*(fai(:,i)-S_mao*z+I0*ddXg(i));



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%LQR mathematical model%%%%%%%%%%%%%%%%%%%%%%%%%
z11(:,i)=z11(:,i-1)+z31(:,i-1)*t;
% z2(:,i)=;

for j=1:1:8
    dx1(j)=z31(j,i);
    v1(j)=z21(j,i);
    dv1(j)=Dy(j)^(-1)*(A1*dx1(j)-beta1*abs(dx1(j))*abs(v1(j))^(n1-1)*v1(j)-gama1*dx1(j)*abs(v1(j))^n1);
   end
z21(:,i)=z21(:,i-1)+dv1*t;

zLQR=[z11(:,i-1)' z21(:,i-1)' z31(:,i-1)']';
z31(:,i)=z31(:,i-1)+([-inv(M)*Ke -inv(M)*KI -inv(M)*C]*zLQR-I0*ddXg(i)+inv(M)*H*u1(:,i-1))*t;

accLQR(:,i)=[-inv(M)*Ke -inv(M)*KI -inv(M)*C]*zLQR-I0*ddXg(i)+inv(M)*H*u1(:,i-1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%LQR%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A=[zeros(8,8) eye(8);
    -inv(M)*Ke -inv(M)*C];
B=[zeros(8,8);
    inv(M)*H];
Q=80*[Ke zeros(8,8);zeros(8,8) M];
% Q=(Q+Q')/2;
R=10^(-6)*eye(8);
%  Q=eye(16);
% R=10^(-7)*eye(8);
G=lqr(A,B,Q,R);
u1(:,i)=-G*[z11(:,i)' z31(:,i)']';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%no control mathematical model%%%%%%%%%%%%%%%%%%%%%%%%%
z12(:,i)=z12(:,i-1)+z32(:,i-1)*t;
% v2=zeros(8,1);
for j=1:1:8
    dx2(j)=z32(j,i);
    v2(j)=z22(j,i);
    dv2(j)=Dy(j)^(-1)*(A1*dx2(j)-beta1*abs(dx2(j))*abs(v2(j))^(n1-1)*v2(j)-gama1*dx2(j)*abs(v2(j))^n1);
end
z22(:,i)=z22(:,i-1)+dv2*t;

zWukong=[z12(:,i-1)' z22(:,i-1)' z32(:,i-1)']';
z32(:,i)=z32(:,i-1)+([-inv(M)*Ke -inv(M)*KI -inv(M)*C]*zWukong-I0*ddXg(i))*t;

accWukong(:,i)=[-inv(M)*Ke -inv(M)*KI -inv(M)*C]*zWukong-I0*ddXg(i);

% for j=1:1:8
% if u(j,i) >3.456*10^5*8*0.1
%     u(j,i)=3.456*10^5*8*0.1;
% elseif u(j,i) <-3.456*10^5*8*0.1
%     u(j,i)=-3.456*10^5*8*0.1;
% else
% u(j,i)=u(j,i);
% end 
% end
% 
% for j=1:1:8
% if u1(j,i) >3.456*10^5*8*0.1
%     u1(j,i)=3.456*10^5*8*0.1;
% elseif u1(j,i) <-3.456*10^5*8*0.1
%     u1(j,i)=-3.456*10^5*8*0.1;
% else
% u1(j,i)=u1(j,i);
% end 
% end


% for j=1:1:8
% if j==1
%     u(j,i)=u(j,i);
% else
%     u(j,i)=0;
% end 
% end

    end
% end

%%%%%%%%%%simulation%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1)
    plot([0.02:0.02:buchanglong],u(1,:),'k','linewidth',1.5)
    hold on
    plot([0.02:0.02:buchanglong],u1(1,:),'b--','linewidth',1.5)
    hold off
    xlabel ('t(s)','fontsize',20)
    ylabel({'First Floor';'Controlling Force(N)'},'fontsize',20)
h=legend('nonsmooth control','LQR control');
set(h,'Fontsize',14);
set(gca,'FontSize',14); 


figure(2)
    subplot(2,1,1)
    plot([0.02:0.02:buchanglong],z1(1,:),'k','linewidth',1.5)
    hold on
    plot([0.02:0.02:buchanglong],z11(1,:),'b--','linewidth',1.5)
    hold off
    xlabel ('t(s)','fontsize',20)
    ylabel({'First Floor';'Displentment(m)'},'fontsize',20)
h=legend('nonsmooth control','LQR control');
set(h,'Fontsize',14);
set(gca,'FontSize',14); 
    subplot(2,1,2)
    plot([0.02:0.02:buchanglong],z12(1,:),'r','linewidth',1.5)
    xlabel ('t(s)','fontsize',20)
    ylabel({'First Floor';'Displentment(m)'},'fontsize',20)
h=legend('no control');
set(h,'Fontsize',14);
set(gca,'FontSize',14); 


    figure(3)
    plot([0.02:0.02:buchanglong],z3(1,:),'k','linewidth',1.5)
    hold on
    plot([0.02:0.02:buchanglong],z31(1,:),'b--','linewidth',1.5)
     plot([0.02:0.02:buchanglong],z32(1,:),'r-.','linewidth',1.5)
    hold off
    xlabel ('t(s)','fontsize',20)
    ylabel({'First Floor';'Velocity(m/s)'},'fontsize',20)
h=legend('nonsmooth control','LQR control','no control');
set(h,'Fontsize',14);
set(gca,'FontSize',14); 
    
    
    figure(4)
    plot([0.02:0.02:buchanglong],fai(1,:),'k','linewidth',1.5)
    hold on
    plot([0.02:0.02:buchanglong],accLQR(1,:),'r','linewidth',1.5)
     plot([0.02:0.02:buchanglong],accWukong(1,:),'g--','linewidth',1.5)
    hold off
    xlabel ('t(s)','fontsize',20)
    ylabel({'First Floor';'Acceleration(m/s^2)'},'fontsize',20)
h=legend('nonsmooth control','LQR control','no control');
set(h,'Fontsize',14);
set(gca,'FontSize',14); 


    figure(5)
    plot(acc1,acc2)
     xlabel ('t(s)','fontsize',20)
    ylabel({'EI Centro';'Acceleration(m/s^2)'},'fontsize',20)
% h=legend('nonsmooth control','LQR control','no control');
% set(h,'Fontsize',14);
set(gca,'FontSize',14); 



table1=zeros(8,10);
table1(:,1)=[1,2,3,4,5,6,7,8]';
table1(:,2)=Dy*0.01;

for i=1:1:8
%     for j=1:1:8
     table1(i,3)=max(abs(z12(i,:)));
     table1(i,4)=max(abs(accWukong(i,:)));
     table1(i,5)=max(abs(z11(i,:)));
     table1(i,6)=max(abs(accLQR(i,:)));
     table1(i,7)=max(abs(z1(i,:)));
     table1(i,8)=max(abs(fai(i,:)));
     
     table1(i,9)=(1-max(abs(z1(i,:)))/max(abs(z11(i,:))))*100;
     table1(i,10)=(1-max(abs(fai(i,:)))/max(abs(accLQR(i,:))))*100;
end
fprintf('result，1st col floor, 2nd col Dy, 3rd col max displacement under no control, 4th col max acc under no control, 5th max displacement LQR, 6th col max acc under no control, 7th max displacement non smooth, 8th col max acc under no control')
table1
fprintf('max control force under non smooth')
round(max(max(u))/10^3)
fprintf('max control force under LQR')
round(max(max(u1))/10^3)

% 
fprintf('percentage of displacement')
table1(:,9)

fprintf('percentage of acc')
table1(:,10)



%   figure(6)
% plot(table1(:,7),table1(:,1),'b')      
% hold on 
% plot(table1(:,5),table1(:,1),'k--')   
% % plot(table1(:,3),table1(:,1),'k-.')   
% hold off
%  xlabel ('Max Displentment(m)','fontsize',20)
%     ylabel({'story number'},'fontsize',20)
% h=legend('nonsmooth control','LQR control','no control');
% set(h,'Fontsize',14);
% set(gca,'FontSize',14); 
