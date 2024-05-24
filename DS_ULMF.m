clear;close all;clc;
N=5;
L=4000;
C=5;
UU=randn(N,C*L);
w0=randn(N,1);

Pup1=0.01;
Pup2=0.02;
Pup3=0.05;
Pup4=0.1;
Pup5=0.2;
Pup6=0.4;
%% Update rates calculation

Pup_inverse1=1-Pup1;
Pup_inverse2=1-Pup2;
Pup_inverse3=1-Pup3;
Pup_inverse4=1-Pup4;
Pup_inverse5=1-Pup5;
Pup_inverse6=1-Pup6;
syms t x

ft=1/2/sqrt(3);
gx=int(ft,t,-x,x);
eqn1 = gx==Pup_inverse1;
eqn2 = gx==Pup_inverse2;
eqn3 = gx==Pup_inverse3;
eqn4 = gx==Pup_inverse4;
eqn5 = gx==Pup_inverse5;
eqn6 = gx==Pup_inverse6;

x1=vpasolve(eqn1,x);
x2=vpasolve(eqn2,x);
x3=vpasolve(eqn3,x);
x4=vpasolve(eqn4,x);
x5=vpasolve(eqn5,x);
x6=vpasolve(eqn6,x);
Threshold1=double(x1);
Threshold2=double(x2);
Threshold3=double(x3);
Threshold4=double(x4);
Threshold5=double(x5);
Threshold6=double(x6);

%% Step Sizes
mu_LMF=0.002;
mu_LMF1=0.0025;
mu_LMF2=0.0025;
mu_LMF3=0.002;
mu_LMF4=0.002;
mu_LMF5=0.002;
mu_LMF6=0.002;
t=2;
for n=1:C*L
    VV2(n)=t*rand-0.5*t;
end
sigma=sqrt(t^2/12);
VV2 = VV2.^6;
s = sum(VV2);
j = s/(C*L);
%% Misadjustments
rol_LMF1= mu_LMF1*2*j*Pup1*N/3/sigma^4;
rol_LMF2= mu_LMF2*2*j*Pup2*N/3/sigma^4;
rol_LMF3= mu_LMF3*2*j*Pup3*N/3/sigma^4;
rol_LMF4= mu_LMF4*2*j*Pup4*N/3/sigma^4;
rol_LMF5= mu_LMF5*2*j*Pup5*N/3/sigma^4;
rol_LMF6= mu_LMF6*2*j*Pup6*N/3/sigma^4;


%% Operating


tic
KK=100;

t=2;
sigma=sqrt(t^2/12);

ek_LMS1=zeros(1,C*L);
o1=zeros(1,KK);
o2=zeros(1,KK);
o3=zeros(1,KK);
o4=zeros(1,KK);
o5=zeros(1,KK);
o6=zeros(1,KK);

o11=zeros(1,KK);
o22=zeros(1,KK);
o33=zeros(1,KK);
o44=zeros(1,KK);
o55=zeros(1,KK);
o66=zeros(1,KK);
for kk=1:KK
    delta1=zeros(C*L,3);
    delta2=zeros(1,C*L);
    delta3=zeros(1,C*L);
    delta4=zeros(1,C*L);
    delta5=zeros(1,C*L);
    delta6=zeros(1,C*L);

    %% Noise

    for n=1:C*L

        VV(n)=t*rand-0.5*t;
    end
    %
    sigma1=sqrt(((t/2)^3-Threshold1^3*sigma^3)/(3*t/2));
    sigma2=sqrt(((t/2)^3-Threshold2^3*sigma^3)/(3*t/2));
    sigma3=sqrt(((t/2)^3-Threshold3^3*sigma^3)/(3*t/2));
    sigma4=sqrt(((t/2)^3-Threshold4^3*sigma^3)/(3*t/2));
    sigma5=sqrt(((t/2)^3-Threshold5^3*sigma^3)/(3*t/2));
    %         sigma6=sqrt(((t/2)^3-Threshold6^3*sigma^3)/(3*t/2));
    %%  Desired Output
    DD=w0'*UU+VV;

    %% Initial system parameter
    w_initial=zeros(N,1);
    w_LMF=w_initial;
    w_LMF1=w_initial;
    w_LMF2=w_initial;
    w_LMF3=w_initial;
    w_LMF4=w_initial;
    w_LMF5=w_initial;
    w_LMF6=w_initial;
    w_LMP=w_initial;

    for i=1:C*L
        dk=DD(i);
        uk=UU(:,i);
        %% LMF
        Err_LMF(kk,i)=(w_LMF-w0)'*(w_LMF-w0);
        ek_LMF = dk- w_LMF'*uk;
        w_LMF=w_LMF+mu_LMF*(ek_LMF)^3*uk;
        %% DS-ULMF1
        Err_LMF1(kk,i)=(w_LMF1-w0)'*(w_LMF1-w0);
        ek_LMF1 = dk- w_LMF1'*uk;
        if  abs( ek_LMF1/sigma)>(Threshold1*(1+rol_LMF1))

            delta=1;
            o11(kk)=o11(kk)+1;
        else
            delta=0;
            o11(kk)=o11(kk);
        end

        w_LMF1=w_LMF1+mu_LMF1*delta*(ek_LMF1)^3*uk;
        %% DS-ULMF2
        Err_LMF2(kk,i)=(w_LMF2-w0)'*(w_LMF2-w0);
        ek_LMF2 = dk- w_LMF2'*uk;
        if  abs( ek_LMF2/sigma)>(Threshold2*(1+rol_LMF2))

            delta=1;
            o22(kk)=o22(kk)+1;
        else
            delta=0;
            o22(kk)=o22(kk);
        end

        w_LMF2=w_LMF2+mu_LMF2*delta*(ek_LMF2)^3*uk;
        %% DS-ULMF3
        Err_LMF3(kk,i)=(w_LMF3-w0)'*(w_LMF3-w0);
        ek_LMF3 = dk- w_LMF3'*uk;
        if  abs( ek_LMF3/sigma)>(Threshold3*(1+rol_LMF3))

            delta=1;
            o33(kk)=o33(kk)+1;
        else
            delta=0;
            o33(kk)=o33(kk);
        end
        w_LMF3=w_LMF3+mu_LMF3*delta*(ek_LMF3)^3*uk;
        %% DS-ULMF4
        Err_LMF4(kk,i)=(w_LMF4-w0)'*(w_LMF4-w0);
        ek_LMF4 = dk- w_LMF4'*uk;
        if  abs( ek_LMF4/sigma)>(Threshold4*(1+rol_LMF4))

            delta=1;
            o44(kk)=o44(kk)+1;
        else
            delta=0;
            o44(kk)=o44(kk);
        end

        w_LMF4=w_LMF4+mu_LMF4*delta*(ek_LMF4)^3*uk;
        %% DS-ULMF5
        Err_LMF5(kk,i)=(w_LMF5-w0)'*(w_LMF5-w0);
        ek_LMF5 = dk- w_LMF5'*uk;
        if  abs( ek_LMF5/sigma)>(Threshold5*(1+rol_LMF5))

            delta=1;
            o55(kk)=o55(kk)+1;
        else
            delta=0;
            o55(kk)=o55(kk);
        end

        w_LMF5=w_LMF5+mu_LMF5*delta*(ek_LMF5)^3*uk;
        %% DS-ULMF6
        Err_LMF6(kk,i)=(w_LMF6-w0)'*(w_LMF6-w0);
        ek_LMF6 = dk- w_LMF6'*uk;
        if  abs( ek_LMF6/sigma)>(Threshold6*(1+rol_LMF6))

            delta=1;
            o66(kk)=o66(kk)+1;
        else
            delta=0;
            o66(kk)=o66(kk);
        end

        w_LMF6=w_LMF6+mu_LMF6*delta*(ek_LMF6)^3*uk;
    end

    disp(kk);
end
toc

gap=1000;
x= 1:gap:C*L;
figure(1),hold on

Err8=mean(Err_LMF1);
Err9=mean(Err_LMF2);
Err10=mean(Err_LMF3);
Err11=mean(Err_LMF4);
Err12=mean(Err_LMF5);
Err13=mean(Err_LMF6);
Err14=mean(Err_LMF);

plot(x,10* log10(Err14(1:gap:C*L)),'o','MarkerSize',8,'LineWidth',1,'Color','g');
plot(x,10* log10(Err13(1:gap:C*L)),'p','MarkerSize',8,'LineWidth',1,'Color','m');
plot(x,10* log10(Err12(1:gap:C*L)),'p','MarkerSize',8,'LineWidth',1,'Color','k');
plot(x,10* log10(Err11(1:gap:C*L)),'s','MarkerSize',8,'LineWidth',1,'Color','y');
plot(x,10* log10(Err10(1:gap:C*L)),'+','MarkerSize',8,'LineWidth',1,'Color','b');
plot(x,10* log10(Err9(1:gap:C*L)),'v','MarkerSize',8,'LineWidth',1,'Color','c');
plot(x,10* log10(Err8(1:gap:C*L)),'h','MarkerSize',8,'LineWidth',1,'Color','r');



plot(10*log10(mean(Err_LMF)),'g','linewidth',1);
plot(10*log10(mean(Err_LMF6)),'Color','m','linewidth',1);
plot(10*log10(mean(Err_LMF5)),'Color','k','linewidth',1);
plot(10*log10(mean(Err_LMF4)),'Color','y','linewidth',1);
plot(10*log10(mean(Err_LMF3)),'Color','b','linewidth',1);
plot(10*log10(mean(Err_LMF2)),'Color','c','linewidth',1);
plot(10*log10(mean(Err_LMF1)),'Color','r','linewidth',1);

legend('LMF (Pup=100%,\eta=0.002)','DS-ULMF (Pup=40%,\eta=0.002)','DS-ULMF (Pup=20%,\eta=0.002)','DS-ULMF(Pup=10%,\eta=0.002)','DS-ULMF(Pup=5%,\eta=0.002)','DS-ULMF(Pup=2%,\eta=0.0025)','DS-ULMF(Pup=1%,\eta=0.0025)');

ylim([-43,8]);
xlabel('Iteration');
ylabel('MSD(dB)');
box on;
grid on;


figure(2)
[f,xi]=ksdensity(VV);
plot(xi,f,'-r');
grid on;
xlabel('Magnitude');
ylabel('Possibility');

op1=mean(o1);
op2=mean(o2);
op3=mean(o3);
op4=mean(o4);
op5=mean(o5);


op11=mean(o11);
op22=mean(o22);
op33=mean(o33);
op44=mean(o44);
op55=mean(o55);
op66=mean(o66);

sigma*Threshold1*(1+rol_LMF1)
sigma*Threshold2*(1+rol_LMF2)
sigma*Threshold3*(1+rol_LMF3)
sigma*Threshold4*(1+rol_LMF4)
sigma*Threshold5*(1+rol_LMF5)

op1/C/L
op2/C/L
op3/C/L
op4/C/L
op5/C/L

op11/C/L
op22/C/L
op33/C/L
op44/C/L
op55/C/L
op66/C/L