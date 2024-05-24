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
Pup7=0.2;
Pup8=0.4;
%% Update rates calculation

Pup_inverse1=1-Pup1;
Pup_inverse2=1-Pup2;
Pup_inverse3=1-Pup3;
Pup_inverse4=1-Pup4;
Pup_inverse5=1-Pup5;
Pup_inverse6=1-Pup6;
Pup_inverse7=1-Pup7;
Pup_inverse8=1-Pup8;
syms t x
ft=1/2/sqrt(3);
gx=int(ft,t,-x,x);
eqn1 = gx==Pup_inverse1;
eqn2 = gx==Pup_inverse2;
eqn3 = gx==Pup_inverse3;
eqn4 = gx==Pup_inverse4;
eqn5 = gx==Pup_inverse5;
eqn6 = gx==Pup_inverse6;
eqn7 = gx==Pup_inverse7;
eqn8 = gx==Pup_inverse8;
x1=vpasolve(eqn1,x);
x2=vpasolve(eqn2,x);
x3=vpasolve(eqn3,x);
x4=vpasolve(eqn4,x);
x5=vpasolve(eqn5,x);
x6=vpasolve(eqn6,x);
x7=vpasolve(eqn7,x);
x8=vpasolve(eqn8,x);
Threshold1=double(x1);
Threshold2=double(x2);
Threshold3=double(x3);
Threshold4=double(x4);
Threshold5=double(x5);
Threshold6=double(x6);
Threshold7=double(x7);
Threshold8=double(x8);
%% Forgetting Factors
lamda_RLS=0.997;
lamda_RLS1=0.997;
lamda_RLS2=0.997;
lamda_RLS3=0.997;
lamda_RLS4=0.997;
lamda_RLS5=0.997;
lamda_RLS6=0.997;
lamda_RLS7=0.997;
lamda_RLS8=0.997;
%% Misadjustments
rol_RLS1= Pup1*N*(1-lamda_RLS1)/(2-(1-lamda_RLS1)*Pup1);
rol_RLS2= Pup2*N*(1-lamda_RLS2)/(2-(1-lamda_RLS2)*Pup2);
rol_RLS3= Pup3*N*(1-lamda_RLS3)/(2-(1-lamda_RLS3)*Pup3);
rol_RLS4= Pup4*N*(1-lamda_RLS4)/(2-(1-lamda_RLS4)*Pup4);
rol_RLS5= Pup5*N*(1-lamda_RLS5)/(2-(1-lamda_RLS5)*Pup5);
rol_RLS6= Pup6*N*(1-lamda_RLS6)/(2-(1-lamda_RLS6)*Pup6);
rol_RLS7= Pup6*N*(1-lamda_RLS7)/(2-(1-lamda_RLS7)*Pup6);
rol_RLS8= Pup6*N*(1-lamda_RLS8)/(2-(1-lamda_RLS8)*Pup6);
%% Operating

tic
KK=50;

t=2;
sigma=sqrt(t^2/12);

ek_LMF1=zeros(1,C*L);
o1=zeros(1,KK);
o2=zeros(1,KK);
o3=zeros(1,KK);
o4=zeros(1,KK);
o5=zeros(1,KK);
o6=zeros(1,KK);


for kk=1:KK
    %% DS-UPDF Noise
    delta1=zeros(C*L,3);
    delta2=zeros(1,C*L);
    delta3=zeros(1,C*L);
    delta4=zeros(1,C*L);
    delta5=zeros(1,C*L);
    delta6=zeros(1,C*L);




    for n=1:C*L

        VV(n)=t*rand-0.5*t;
    end
    V1=(abs(VV)>((1-Pup1)*t/2)).*VV;
    V2=(abs(VV)>((1-Pup2)*t/2)).*VV;
    V3=(abs(VV)>((1-Pup3)*t/2)).*VV;
    V4=(abs(VV)>((1-Pup4)*t/2)).*VV;
    V5=(abs(VV)>((1-Pup5)*t/2)).*VV;
    V6=(abs(VV)>((1-Pup6)*t/2)).*VV;
    V7=(abs(VV)>((1-Pup7)*t/2)).*VV;
    V8=(abs(VV)>((1-Pup8)*t/2)).*VV;

    sigma1=sqrt(((t/2)^3-(Threshold1*sigma)^3)/3/(t/2));
    sigma2=sqrt(((t/2)^3-(Threshold2*sigma)^3)/3/(t/2));
    sigma3=sqrt(((t/2)^3-(Threshold3*sigma)^3)/3/(t/2));
    sigma4=sqrt(((t/2)^3-(Threshold4*sigma)^3)/3/(t/2));
    sigma5=sqrt(((t/2)^3-(Threshold5*sigma)^3)/3/(t/2));
    sigma6=sqrt(((t/2)^3-(Threshold6*sigma)^3)/3/(t/2));
    sigma7=sqrt(((t/2)^3-(Threshold7*sigma)^3)/3/(t/2));
    sigma8=sqrt(((t/2)^3-(Threshold8*sigma)^3)/3/(t/2));

    %%  Desired Output
    DD=w0'*UU+VV;
    D1=w0'*UU+V1;
    D2=w0'*UU+V2;
    D3=w0'*UU+V3;
    D4=w0'*UU+V4;
    D5=w0'*UU+V5;
    D6=w0'*UU+V6;
    D7=w0'*UU+V7;
    D8=w0'*UU+V8;
    %% Initial system parameter

    w_initial=randn(N,1);
    w_RLS=w_initial;
    w_RLS1=w_initial;
    w_RLS2=w_initial;
    w_RLS3=w_initial;
    w_RLS4=w_initial;
    w_RLS5=w_initial;
    w_RLS6=w_initial;
    w_RLS7=w_initial;
    w_RLS8=w_initial;

    Pn = eye(N)*1;
    Pn1 = eye(N)*1;
    Pn2 = eye(N)*1;
    Pn3 = eye(N)*1;
    Pn4 = eye(N)*1;
    Pn5 = eye(N)*1;
    Pn6 = eye(N)*1;
    Pn7 = eye(N)*1;
    Pn8 = eye(N)*1;
    for i=1:C*L
        dk=DD(i);
        dk1=D1(i);
        dk2=D2(i);
        dk3=D3(i);
        dk4=D4(i);
        dk5=D5(i);
        dk6=D6(i);
        dk7=D7(i);
        dk8=D8(i);
        uk=UU(:,i);

        %% RLS
        Err_RLS(kk,i) = (w_RLS-w0)' * (w_RLS-w0);
        ek_RLS = dk - w_RLS' * uk;
        kn = Pn * uk / ( lamda_RLS+ uk' * Pn * uk );
        Pn = 1/lamda_RLS * ( Pn - kn * uk' * Pn);
        w_RLS = w_RLS +kn * ek_RLS;
        %% RLS1
        Err_RLS1(kk,i) = (w_RLS1-w0)' * (w_RLS1-w0);
        ek_RLS1 = dk1 - w_RLS1' * uk;
        kn1 = Pn1 * uk / ( lamda_RLS1+ uk' * Pn1 * uk );
        Pn1= 1/lamda_RLS1 * ( Pn1 - kn1 * uk' * Pn1);
        w_RLS1 = w_RLS1 +kn1 * ek_RLS1;
        %% RLS2
        Err_RLS2(kk,i) = (w_RLS2-w0)' * (w_RLS2-w0);
        ek_RLS2 = dk2 - w_RLS2' * uk;
        kn2 = Pn2 * uk / ( lamda_RLS2+ uk' * Pn2 * uk );
        Pn2= 1/lamda_RLS2 * ( Pn2 - kn2 * uk' * Pn2);
        w_RLS2 = w_RLS2 +kn2 * ek_RLS2;
        %% RLS3
        Err_RLS3(kk,i) = (w_RLS3-w0)' * (w_RLS3-w0);
        ek_RLS3 = dk3 - w_RLS3' * uk;
        kn3 = Pn3 * uk / ( lamda_RLS3+ uk' * Pn3 * uk );
        Pn3= 1/lamda_RLS3 * ( Pn3 - kn3 * uk' * Pn3);
        w_RLS3 = w_RLS3+kn3 * ek_RLS3;
        %% RLS4
        Err_RLS4(kk,i) = (w_RLS4-w0)' * (w_RLS4-w0);
        ek_RLS4 = dk4 - w_RLS4' * uk;
        kn4 = Pn4 * uk / ( lamda_RLS4+ uk' * Pn4 * uk );
        Pn4= 1/lamda_RLS4 * ( Pn4 - kn4 * uk' * Pn4);
        w_RLS4 = w_RLS4 +kn1 * ek_RLS4;
        %% RLS5
        Err_RLS5(kk,i) = (w_RLS5-w0)' * (w_RLS5-w0);
        ek_RLS5 = dk5 - w_RLS5' * uk;
        kn5 = Pn5* uk / ( lamda_RLS5+ uk' * Pn5 * uk );
        Pn5= 1/lamda_RLS5 * ( Pn5 - kn5 * uk' * Pn5);
        w_RLS5 = w_RLS5 +kn5 * ek_RLS5;
        %% RLS6
        Err_RLS6(kk,i) = (w_RLS6-w0)' * (w_RLS6-w0);
        ek_RLS6 = dk6 - w_RLS6' * uk;
        kn6 = Pn6 * uk / ( lamda_RLS6+ uk' * Pn6 * uk );
        Pn6= 1/lamda_RLS6 * ( Pn6 - kn6 * uk' * Pn6);
        w_RLS6 = w_RLS6 +kn6 * ek_RLS6;
        %% RLS7
        Err_RLS7(kk,i) = (w_RLS7-w0)' * (w_RLS7-w0);
        ek_RLS7 = dk7 - w_RLS7' * uk;
        kn7 = Pn7 * uk / ( lamda_RLS7+ uk' * Pn7 * uk );
        Pn7= 1/lamda_RLS7 * ( Pn7 - kn6 * uk' * Pn7);
        w_RLS7 = w_RLS7 +kn7 * ek_RLS7;
        %% RLS8
        Err_RLS8(kk,i) = (w_RLS8-w0)' * (w_RLS8-w0);
        ek_RLS8 = dk8 - w_RLS8' * uk;
        kn8 = Pn8 * uk / ( lamda_RLS8+ uk' * Pn8 * uk );
        Pn8= 1/lamda_RLS8 * ( Pn8 - kn8 * uk' * Pn8);
        w_RLS8 = w_RLS8 +kn8 * ek_RLS8;
    end

%% Theoretical Prior-selection SS-MSD
    TH(kk)=(1-lamda_RLS)/(1+lamda_RLS)*N*sigma^2;
    TH1(kk)=(1-lamda_RLS1)/(1+lamda_RLS1)*N*sigma1^2;
    TH2(kk)=(1-lamda_RLS2)/(1+lamda_RLS2)*N*sigma2^2;
    TH3(kk)=(1-lamda_RLS3)/(1+lamda_RLS3)*N*sigma3^2;
    TH4(kk)=(1-lamda_RLS4)/(1+lamda_RLS4)*N*sigma4^2;
    TH5(kk)=(1-lamda_RLS5)/(1+lamda_RLS5)*N*sigma5^2;
    TH6(kk)=(1-lamda_RLS6)/(1+lamda_RLS6)*N*sigma6^2;
    TH7(kk)=(1-lamda_RLS7)/(1+lamda_RLS7)*N*sigma7^2;
    TH8(kk)=(1-lamda_RLS8)/(1+lamda_RLS8)*N*sigma8^2;
    disp(kk);
end
toc

gap=1000;
x= 1:gap:C*L;
figure(1),hold on

Err1=mean(Err_RLS1);
Err2=mean(Err_RLS2);
Err3=mean(Err_RLS3);
Err4=mean(Err_RLS4);
Err5=mean(Err_RLS5);
Err6=mean(Err_RLS6);
% Err7=mean(Err_RLS7);
% Err8=mean(Err_RLS8);

Err9=mean(Err_RLS);

plot(x,10* log10(Err9(1:gap:C*L)),'x','MarkerSize',10,'LineWidth',1,'Color','g');
%  plot(x,10* log10(Err8(1:gap:C*L)),'<','MarkerSize',10,'LineWidth',1,'Color','r');
%  plot(x,10* log10(Err7(1:gap:C*L)),'^','MarkerSize',10,'LineWidth',1,'Color','[0.5 0.5 0.1]');
plot(x,10* log10(Err6(1:gap:C*L)),'o','MarkerSize',10,'LineWidth',1,'Color','m');
plot(x,10* log10(Err5(1:gap:C*L)),'*','MarkerSize',10,'LineWidth',1,'Color','k');
plot(x,10* log10(Err4(1:gap:C*L)),'d','MarkerSize',10,'LineWidth',1,'Color','y');
plot(x,10* log10(Err3(1:gap:C*L)),'s','MarkerSize',10,'LineWidth',1,'Color','b');
plot(x,10* log10(Err2(1:gap:C*L)),'p','MarkerSize',10,'LineWidth',1,'Color','c');
plot(x,10* log10(Err1(1:gap:C*L)),'h','MarkerSize',10,'LineWidth',1,'Color','r');




plot(10*log10(mean(TH)*ones(1,C*L)),'--g','linewidth',2);
%  plot(10*log10(mean(TH8)*ones(1,C*L)),'--','Color','r','linewidth',2);
%  plot(10*log10(mean(TH7)*ones(1,C*L)),'--','Color','[0.5 0.5 0.1]','linewidth',2);
plot(10*log10(mean(TH6)*ones(1,C*L)),'--m','linewidth',2);
plot(10*log10(mean(TH5)*ones(1,C*L)),'--k','linewidth',2);
plot(10*log10(mean(TH4)*ones(1,C*L)),'--y','linewidth',2);
plot(10*log10(mean(TH3)*ones(1,C*L)),'--b','linewidth',2);
plot(10*log10(mean(TH2)*ones(1,C*L)),'--c','linewidth',2);
plot(10*log10(mean(TH1)*ones(1,C*L)),'--r','linewidth',2);

plot(10*log10(mean(Err_RLS)),'g','linewidth',1);
plot(10*log10(mean(Err_RLS1)),'r','linewidth',1);
plot(10*log10(mean(Err_RLS2)),'c','linewidth',1);
plot(10*log10(mean(Err_RLS3)),'b','linewidth',1);
plot(10*log10(mean(Err_RLS4)),'y','linewidth',1);
plot(10*log10(mean(Err_RLS5)),'k','linewidth',1);
plot(10*log10(mean(Err_RLS6)),'m','linewidth',1);
% plot(10*log10(mean(Err_RLS7)),'Color','[0.5 0.5 0.1]','linewidth',1);
% plot(10*log10(mean(Err_RLS8)),'Color','r','linewidth',1);
legend('NumColumns', 2);
legend('RLS (DS-UPDF:Pup=100%)','RLS (DS-UPDF:Pup=40%)','RLS (DS-UPDF:Pup=20%)','RLS (DS-UPDF:Pup=10%)','RLS (DS-UPDF:Pup=5%)','RLS (DS-UPDF:Pup=2%)', ...
    'RLS (DS-UPDF:Pup=1%)','TH(DS-UPDF:Pup=100%)','TH(DS-UPDF:Pup=40%)','TH(DS-UPDF:Pup=20%)','TH(DS-UPDF:Pup=10%)','TH(DS-UPDF:Pup=5%)','TH(DS-UPDF:Pup=2%)','TH(DS-UPDF:Pup=1%)');

ylim([-43,12]);

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


mean(Err_RLS1(:,end))
mean(Err_RLS2(:,end))
mean(Err_RLS3(:,end))
mean(Err_RLS4(:,end))
mean(Err_RLS5(:,end))
mean(Err_RLS6(:,end))