clear;close all;clc;
N=5;
L=4000;
C=5;

UU=randn(N,C*L);
w0=randn(N,1);
alpha=2;
P=alpha;
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
mu_LMS=0.003;
mu_LMS1=0.003;
mu_LMS2=0.003;
mu_LMS3=0.003;
mu_LMS4=0.003;
mu_LMS5=0.003;
mu_LMS6=0.003;
mu_LMS7=0.003;
mu_LMS8=0.003;


fi1=1-(1-Pup1)^3;
fi2=1-(1-Pup2)^3;
fi3=1-(1-Pup3)^3;
fi4=1-(1-Pup4)^3;
fi5=1-(1-Pup5)^3;
fi6=1-(1-Pup6)^3;
%% Misadjustments
rol_LMS1= mu_LMS1*Pup1*N/(1-mu_LMS1*Pup1*N);
rol_LMS2= mu_LMS2*Pup2*N/(1-mu_LMS2*Pup2*N);
rol_LMS3= mu_LMS3*Pup3*N/(1-mu_LMS3*Pup3*N);
rol_LMS4= mu_LMS4*Pup4*N/(1-mu_LMS4*Pup4*N);
rol_LMS5= mu_LMS5*Pup5*N/(1-mu_LMS5*Pup5*N);
rol_LMS6= mu_LMS6*Pup6*N/(1-mu_LMS6*Pup6*N);

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
for kk=1:KK
    %% Noise
    delta1=zeros(C*L,3);
    delta2=zeros(1,C*L);
    delta3=zeros(1,C*L);
    delta4=zeros(1,C*L);
    delta5=zeros(1,C*L);
    delta6=zeros(1,C*L);


    for n=1:C*L

        VV(n)=t*rand-0.5*t;
    end

    sigma1=sqrt(((t/2)^3-Threshold1^3*sigma^3)/(3*t/2));
    sigma2=sqrt(((t/2)^3-Threshold2^3*sigma^3)/(3*t/2));
    sigma3=sqrt(((t/2)^3-Threshold3^3*sigma^3)/(3*t/2));
    sigma4=sqrt(((t/2)^3-Threshold4^3*sigma^3)/(3*t/2));
    sigma5=sqrt(((t/2)^3-Threshold5^3*sigma^3)/(3*t/2));
    sigma6=sqrt(((t/2)^3-Threshold6^3*sigma^3)/(3*t/2));
    %%  Desired Output
    DD=w0'*UU+VV;

    %% Initial system parameter

    w_initial=randn(N,1);
    w_LMS1=w_initial;
    w_LMS2=w_initial;
    w_LMS3=w_initial;
    w_LMS4=w_initial;
    w_LMS5=w_initial;
    w_LMS6=w_initial;
    w_LMS=w_initial;

    for i=1:C*L
        dk=DD(i);
        uk=UU(:,i);

        %% LMS
        Err_LMS(kk,i)=(w_LMS-w0)'*(w_LMS-w0);
        ek_LMS = dk- w_LMS'*uk;
        w_LMS=w_LMS+mu_LMS*ek_LMS*uk;
        %% DS-ULMS1
        Err_LMS1(kk,i)=(w_LMS1-w0)'*(w_LMS1-w0);
        ek_LMS1 = dk- w_LMS1'*uk;
        if  abs( ek_LMS1/sigma)>(Threshold1*(1+rol_LMS1))

            delta1(i,1)=1;
            % delta1(i,2)=ek_LMS1(i);
            delta1(i,3)=VV(i);
            o1(kk)=o1(kk)+1;
        else
            delta1(i,1)=0;
            o1(kk)=o1(kk);
        end
        w_LMS1=w_LMS1+mu_LMS1*delta1(i,1)*uk*ek_LMS1;
        %%  DS-ULMS2
        Err_LMS2(kk,i)=(w_LMS2-w0)'*(w_LMS2-w0);
        ek_LMS2 = dk- w_LMS2'*uk;
        if  abs( ek_LMS2/sigma)>(Threshold2*(1+rol_LMS2))

            delta=1;
            o2(kk)=o2(kk)+1;
        else
            delta=0;
            o2(kk)=o2(kk);
        end

        w_LMS2=w_LMS2+mu_LMS2*delta*uk*ek_LMS2;
        %%  DS-ULMS3
        Err_LMS3(kk,i)=(w_LMS3-w0)'*(w_LMS3-w0);
        ek_LMS3 = dk- w_LMS3'*uk;
        if  abs( ek_LMS3/sigma)>(Threshold3*(1+rol_LMS3))

            delta=1;
            o3(kk)=o3(kk)+1;
        else
            delta=0;
            o3(kk)=o3(kk);
        end

        w_LMS3=w_LMS3+mu_LMS3*delta*uk*ek_LMS3;
        %%  DS-ULMS4
        Err_LMS4(kk,i)=(w_LMS4-w0)'*(w_LMS4-w0);
        ek_LMS4 = dk- w_LMS4'*uk;
        if  abs( ek_LMS4/sigma)>(Threshold4*(1+rol_LMS4))

            delta=1;
            o4(kk)=o4(kk)+1;
        else
            delta=0;
            o4(kk)=o4(kk);
        end

        w_LMS4=w_LMS4+mu_LMS4*delta*uk*ek_LMS4;
        %%  DS-ULMS5
        Err_LMS5(kk,i)=(w_LMS5-w0)'*(w_LMS5-w0);
        ek_LMS5 = dk- w_LMS5'*uk;
        if  abs( ek_LMS5/sigma)>(Threshold5*(1+rol_LMS5))

            delta=1;
            o5(kk)=o5(kk)+1;
        else
            delta=0;
            o5(kk)=o5(kk);
        end

        w_LMS5=w_LMS5+mu_LMS5*delta*uk*ek_LMS5;
        %%  DS-ULMS6
        Err_LMS6(kk,i)=(w_LMS6-w0)'*(w_LMS6-w0);
        ek_LMS6 = dk- w_LMS6'*uk;
        if  abs( ek_LMS6/sigma)>(Threshold6*(1+rol_LMS6))

            delta=1;
            o6(kk)=o6(kk)+1;
        else
            delta=0;
            o6(kk)=o6(kk);
        end

        w_LMS6=w_LMS6+mu_LMS6*delta*uk*ek_LMS6;

    end
    %% Theoretical post-selection SS-MSD

    TH(kk)=mu_LMS^2*N*sigma^2/(1-(1-mu_LMS)^2);
    TH1(kk)=mu_LMS1^2*N*sigma1^2/(1-(1-mu_LMS1)^2);
    TH2(kk)=mu_LMS2^2*N*sigma2^2/(1-(1-mu_LMS2)^2);
    TH3(kk)=mu_LMS3^2*N*sigma3^2/(1-(1-mu_LMS3)^2);
    TH4(kk)=mu_LMS4^2*N*sigma4^2/(1-(1-mu_LMS4)^2);
    TH5(kk)=mu_LMS5^2*N*sigma5^2/(1-(1-mu_LMS5)^2);
    TH6(kk)=mu_LMS6^2*N*sigma6^2/(1-(1-mu_LMS6)^2);
    disp(kk);
end
toc

gap=1000;
x= 1:gap:C*L;
figure(1),hold on

Err1=mean(Err_LMS1);
Err2=mean(Err_LMS2);
Err3=mean(Err_LMS3);
Err4=mean(Err_LMS4);
Err5=mean(Err_LMS5);
Err6=mean(Err_LMS6);

Err7=mean(Err_LMS);

plot(x,10* log10(Err7(1:gap:C*L)),'x','MarkerSize',10,'LineWidth',1,'Color','g');
plot(x,10* log10(Err6(1:gap:C*L)),'o','MarkerSize',10,'LineWidth',1,'Color','m');
plot(x,10* log10(Err5(1:gap:C*L)),'*','MarkerSize',10,'LineWidth',1,'Color','k');
plot(x,10* log10(Err4(1:gap:C*L)),'d','MarkerSize',10,'LineWidth',1,'Color','y');
plot(x,10* log10(Err3(1:gap:C*L)),'s','MarkerSize',10,'LineWidth',1,'Color','b');
plot(x,10* log10(Err2(1:gap:C*L)),'p','MarkerSize',10,'LineWidth',1,'Color','c');
plot(x,10* log10(Err1(1:gap:C*L)),'h','MarkerSize',10,'LineWidth',1,'Color','r');

plot(10*log10(mean(TH)*ones(1,C*L)),'--g','linewidth',2);
plot(10*log10(mean(TH6)*ones(1,C*L)),'--m','linewidth',2);
plot(10*log10(mean(TH5)*ones(1,C*L)),'--k','linewidth',2);
plot(10*log10(mean(TH4)*ones(1,C*L)),'--y','linewidth',2);
plot(10*log10(mean(TH3)*ones(1,C*L)),'--b','linewidth',2);
plot(10*log10(mean(TH2)*ones(1,C*L)),'--c','linewidth',2);
plot(10*log10(mean(TH1)*ones(1,C*L)),'--r','linewidth',2);



plot(10*log10(mean(Err_LMS)),'g','linewidth',1);
plot(10*log10(mean(Err_LMS1)),'r','linewidth',1);
plot(10*log10(mean(Err_LMS2)),'c','linewidth',1);
plot(10*log10(mean(Err_LMS3)),'b','linewidth',1);
plot(10*log10(mean(Err_LMS4)),'y','linewidth',1);
plot(10*log10(mean(Err_LMS5)),'k','linewidth',1);
plot(10*log10(mean(Err_LMS6)),'m','linewidth',1);

legend('NumColumns', 2);
legend('LMS (Pup=100%)','DS-ULMS(Pup=40%)','DS-ULMS (Pup=20%)','DS-ULMS (Pup=10%)','DS-ULMS (Pup=5%)','DS-ULMS (Pup=2%)', ...
    'DS-ULMS (Pup=1%)','TH(Pup=100%)','TH(Pup=40%)','TH(Pup=20%)','TH(Pup=10%)','TH(Pup=5%)','TH(Pup=2%)','TH(Pup=1%)');

ylim([-43,13]);

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
op6=mean(o6);

op1/C/L
op2/C/L
op3/C/L
op4/C/L
op5/C/L
op6/C/L

mean(Err_LMS1(:,end))
mean(Err_LMS2(:,end))
mean(Err_LMS3(:,end))
mean(Err_LMS4(:,end))
mean(Err_LMS5(:,end))
mean(Err_LMS6(:,end))