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
Pup7=0.6;
Pup8=0.8;
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


    sigma1=sqrt(((t/2)^3-Threshold1^3*sigma^3)/(3*t/2));
    sigma2=sqrt(((t/2)^3-Threshold2^3*sigma^3)/(3*t/2));
    sigma3=sqrt(((t/2)^3-Threshold3^3*sigma^3)/(3*t/2));
    sigma4=sqrt(((t/2)^3-Threshold4^3*sigma^3)/(3*t/2));
    sigma5=sqrt(((t/2)^3-Threshold5^3*sigma^3)/(3*t/2));
    sigma6=sqrt(((t/2)^3-Threshold6^3*sigma^3)/(3*t/2));
    sigma7=sqrt(((t/2)^3-Threshold7^3*sigma^3)/(3*t/2));
    sigma8=sqrt(((t/2)^3-Threshold8^3*sigma^3)/(3*t/2));
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

    w_LMS1=w_initial;
    w_LMS2=w_initial;
    w_LMS3=w_initial;
    w_LMS4=w_initial;
    w_LMS5=w_initial;
    w_LMS6=w_initial;
    w_LMS7=w_initial;
    w_LMS8=w_initial;
    w_LMS=w_initial;

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

        %% LMS
        Err_LMS(kk,i)=(w_LMS-w0)'*(w_LMS-w0);
        ek_LMS = dk- w_LMS'*uk;
        w_LMS=w_LMS+mu_LMS*abs(ek_LMS)^(P-1)*sign(ek_LMS)*uk;
        %% LMS1
        Err_LMS1(kk,i)=(w_LMS1-w0)'*(w_LMS1-w0);
        ek_LMS1(i) = dk1- w_LMS1'*uk;
        w_LMS1=w_LMS1+mu_LMS1*abs(ek_LMS1(i))^(alpha-2)*uk*ek_LMS1(i);
        %% LMS2
        Err_LMS2(kk,i)=(w_LMS2-w0)'*(w_LMS2-w0);
        ek_LMS2 = dk2- w_LMS2'*uk;
        w_LMS2=w_LMS2+mu_LMS2*abs(ek_LMS2)^(alpha-2)*uk*ek_LMS2;
        %% LMS3
        Err_LMS3(kk,i)=(w_LMS3-w0)'*(w_LMS3-w0);
        ek_LMS3 = dk3- w_LMS3'*uk;
        w_LMS3=w_LMS3+mu_LMS3*abs(ek_LMS3)^(alpha-2)*uk*ek_LMS3;
        %% LMS4
        Err_LMS4(kk,i)=(w_LMS4-w0)'*(w_LMS4-w0);
        ek_LMS4 = dk4- w_LMS4'*uk;
        w_LMS4=w_LMS4+mu_LMS4*abs(ek_LMS4)^(alpha-2)*uk*ek_LMS4;
        %% LMS5
        Err_LMS5(kk,i)=(w_LMS5-w0)'*(w_LMS5-w0);
        ek_LMS5 = dk5- w_LMS5'*uk;
        w_LMS5=w_LMS5+mu_LMS5*abs(ek_LMS5)^(alpha-2)*uk*ek_LMS5;
        %% LMS6
        Err_LMS6(kk,i)=(w_LMS6-w0)'*(w_LMS6-w0);
        ek_LMS6 = dk6- w_LMS6'*uk;
        w_LMS6=w_LMS6+mu_LMS6*abs(ek_LMS6)^(alpha-2)*uk*ek_LMS6;

        %    %% LMS7
        %         Err_LMS7(kk,i)=(w_LMS7-w0)'*(w_LMS7-w0);
        %         ek_LMS7 = dk7- w_LMS7'*uk;
        %     w_LMS7=w_LMS7+mu_LMS7*abs(ek_LMS7)^(alpha-2)*uk*ek_LMS7;

        %        %% LMS8
        %         Err_LMS8(kk,i)=(w_LMS8-w0)'*(w_LMS8-w0);
        %         ek_LMS8 = dk8- w_LMS8'*uk;
        %     w_LMS8=w_LMS8+mu_LMS8*abs(ek_LMS8)^(alpha-2)*uk*ek_LMS8;
    end

%% Theoretical Prior-selection SS-MSD

    TH(kk)=mu_LMS^2*N*sigma^2/(1-(1-mu_LMS)^2);
    TH1(kk)=mu_LMS1^2*N*sigma1^2/(1-(1-mu_LMS1)^2);
    TH2(kk)=mu_LMS2^2*N*sigma2^2/(1-(1-mu_LMS2)^2);
    TH3(kk)=mu_LMS3^2*N*sigma3^2/(1-(1-mu_LMS3)^2);
    TH4(kk)=mu_LMS4^2*N*sigma4^2/(1-(1-mu_LMS4)^2);
    TH5(kk)=mu_LMS5^2*N*sigma5^2/(1-(1-mu_LMS5)^2);
    TH6(kk)=mu_LMS6^2*N*sigma6^2/(1-(1-mu_LMS6)^2);

    %    TH7(kk)=mu_LMS7^2*N*sigma7^2*y^2/(1-(1-mu_LMS7*y^2)^2);
    %    TH8(kk)=mu_LMS8^2*N*sigma8^2*y^2/(1-(1-mu_LMS8*y^2)^2);
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
% Err7=mean(Err_LMS7);
% Err8=mean(Err_LMS8);
Err9=mean(Err_LMS);

plot(x,10* log10(Err9(1:gap:C*L)),'x','MarkerSize',10,'LineWidth',1,'Color','g');
% plot(x,10* log10(Err8(1:gap:C*L)),'<','MarkerSize',10,'LineWidth',1,'Color','[0.6 0.2 0.1]');
% plot(x,10* log10(Err7(1:gap:C*L)),'^','MarkerSize',10,'LineWidth',1,'Color','[0.5 0.5 0.1]');
plot(x,10* log10(Err6(1:gap:C*L)),'o','MarkerSize',10,'LineWidth',1,'Color','m');
plot(x,10* log10(Err5(1:gap:C*L)),'*','MarkerSize',10,'LineWidth',1,'Color','k');
plot(x,10* log10(Err4(1:gap:C*L)),'d','MarkerSize',10,'LineWidth',1,'Color','y');
plot(x,10* log10(Err3(1:gap:C*L)),'s','MarkerSize',10,'LineWidth',1,'Color','b');
plot(x,10* log10(Err2(1:gap:C*L)),'p','MarkerSize',10,'LineWidth',1,'Color','c');
plot(x,10* log10(Err1(1:gap:C*L)),'h','MarkerSize',10,'LineWidth',1,'Color','r');

plot(10*log10(mean(TH)*ones(1,C*L)),'--g','linewidth',2);

% plot(10*log10(mean(TH8)*ones(1,C*L)),'--','linewidth',2,'Color','[0.6 0.2 0.1]');
% plot(10*log10(mean(TH7)*ones(1,C*L)),'--','linewidth',2,'Color','[0.5 0.5 0.1]');
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
% plot(10*log10(mean(Err_LMS7)),'Color','[0.5 0.5 0.1]','linewidth',1);
% plot(10*log10(mean(Err_LMS8)),'Color','[0.6 0.2 0.1]','linewidth',1);

legend('NumColumns', 2);
legend('LMS (DS-UPDF:Pup=100%)','LMS(DS-UPDF:Pup=40%)','LMS (DS-UPDF:Pup=20%)','LMS (DS-UPDF:Pup=10%)','LMS (DS-UPDF:Pup=5%)','LMS (DS-UPDF:Pup=2%)', ...
    'LMS (DS-UPDF:Pup=1%)','TH(DS-UPDF:Pup=100%)','TH(DS-UPDF:Pup=40%)','TH(DS-UPDF:Pup=20%)','TH(DS-UPDF:Pup=10%)','TH(DS-UPDF:Pup=5%)','TH(DS-UPDF:Pup=2%)','TH(DS-UPDF:Pup=1%)');

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


mean(Err_LMS1(:,end))
mean(Err_LMS2(:,end))
mean(Err_LMS3(:,end))
mean(Err_LMS4(:,end))
mean(Err_LMS5(:,end))
mean(Err_LMS6(:,end))