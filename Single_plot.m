% System Chloroform
% This program is to get PV line for T = 500 K only
% to analyze cubic behavior of Van der Waals equation in liquid-vapor coexistence region
R=0.0831451; % in  l·bar·mol^−1·K^−1
a=15.34; % in bar·l^2·mol^-2
b=0.1019; % in l·mol^-1
T=500;
Vl=0.1985; % refer Pdata for 500 K from main program
Vg= 0.5927;% refer Pdata for 500 K from main program
Psat= 41.0360;% refer Pdata for 500 K from main program
V=zeros(1200,1);
P=zeros(1200,1);
for i=1:1:1200
    V(i,1)=0.120+i/1000; %making array for volume at step of 0.001 l
    P(i,1)= (R*T/(V(i,1)-b))-(a/(V(i,1)*V(i,1)));% Finding Pressure for corresponding volume at given temperature using EOS
end
plot(V,P); %Plotting PV Line
ylim([0 90])
title("P-V Diagram for Chloroform Using Van der Waals Equation");
xlabel("V in litres");
ylabel("P in bar");
hold on;
V=zeros(0,0);
P=zeros(0,0);
jjj=1;
for i=0.120:0.001:Vl
    V=[V i];
    P= [P (R*T/(V(1,jjj)-b))-(a/(V(1,jjj)*V(1,jjj)))];
    jjj=jjj+1;
end
for i=Vl:0.001:Vg
    V=[V i];
    P=[P Psat];
    jjj=jjj+1;
end
plot(V,P);
ylim([0 90])