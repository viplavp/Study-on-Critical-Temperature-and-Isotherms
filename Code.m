% System Chloroform using Van der Waals equation of State
R=0.0831451; % in  l·bar·mol^−1·K^−1
a=15.34; % in bar·l^2·mol^-2
b=0.1019; % in l·mol^-1
Pc=53.2868; % in bar
Vc=0.241;   % in l
Tc=8*a/(27*R*b); %in K using formular derived from Van der Waal EOS
jj=1;
%printing Critical point
Pc
Vc
Tc
for T= 455:5:525  % loop from temperature of 455 K to 525 K with step of 5 K
    V=zeros(1200,1);
    P=zeros(1200,1);
    for i=1:1:1200
        V(i,1)=0.120+i/1000; %making array for volume at step of 0.001 litres
        P(i,1)= (R*T/(V(i,1)-b))-(a/(V(i,1)*V(i,1)));% Finding Pressure for corresponding volume at given temperature using EOS
    end
    plot(V,P); %Plotting P-V Diagram
    ylim([0 90])
    title("P-V Diagram for Chloroform Using Van der Waals Equation");
    xlabel("V in litres");
    ylabel("P in bar");
    hold on;
end
Parr1=zeros(1,18); %Parr1, Parr2, Varr1 and Varr2  would store saturation data used for dome plotting
Parr2=zeros(1,18);
Varr1=zeros(1,18);
Varr2=zeros(1,18);
Pdata=zeros(17,3);
for T= 445:5:525  % loop from temperature of 445 K to 525 K with step of 5 K
    % Psat calculation Part
    % Refer Assignment Report to know more about function F
    F = @(x) [(R*T*log(x(1)/(x(1)-b))+b*R*T/(x(1)-b)-R*T*log(x(1)/(x(1)-b)-a/(x(1)*R*T))-2*a/x(1)-R*T*log(x(2)/(x(2)-b))-b*R*T/(x(2)-b)+R*T*log(x(2)/(x(2)-b)-a/(x(2)*R*T))+2*a/x(2));
        (x(1)^3-(b+R*T/x(3))*x(1)^2+a*x(1)/x(3)-a*b/x(3));
        (x(2)^3-(b+R*T/x(3))*x(2)^2+a*x(2)/x(3)-a*b/x(3))];
    %Taking initial guess of Psat using zero derivative points as discussed
    %in the Report of Assignment
    roo = roots([R*T -2*a 4*a*b -2*a*b*b]);
    pres1= R*T/(roo(1)-b)-a/(roo(1)^2);  %pres1 corresponds to maxima point pressure
    pres2= R*T/(roo(2)-b)-a/(roo(2)^2);  %pres2 corresponds to minima point pressure
    pres=(pres1+pres2)/2; %Taking Average of minima pressure and maxima pressure to get initial guess for Psat
    vr = roots([pres -pres*b-R*T a -a*b]);% Finding Volume corresponding to pres
    x0 = [min(vr) max(vr) pres]; %Initial Guess Vector in form of [Vl Vg Psat]
    [x,fval] = fsolve(F,x0) % Solving three simutaneous Non linear equations defined by the function F
    Parr1(1,jj)=x(3);
    Parr2(1,19-jj)=x(3);
    Varr1(1,jj)=x(1);
    Varr2(1,19-jj)=x(2);
    Pdata(jj,:)=x;   % store saturation data
    jj=jj+1;
    
end
Parr1(1,jj)=Pc;
Varr1(1,jj)=Vc;
Parr2(1,1)=Pc;
Varr2(1,1)=Vc;
for T= 531:6:567 %Plotting near-critical and supercritical region
    V=zeros(1200,1);
    P=zeros(1200,1);
    for i=1:1:1200
        V(i,1)=0.120+i/1000;
        P(i,1)= (R*T/(V(i,1)-b))-(a/(V(i,1)*V(i,1)));
    end
    plot(V,P);
    ylim([0 90])
    hold on;
end
%Plotting the splines for the DOME
xx1 = Varr1(1,1)-0.005:0.01:Varr1(1,18);
yy1 = spline(Varr1,Parr1,xx1);  % Make spline
p=plot(xx1,yy1); % Plot spline
p.Color = 'k'; %set colour of dome to black
p.LineStyle = '--'; %make dome curve as broken line
xx2 = Varr2(1,1):0.01:Varr2(1,18);
yy2 = spline(Varr2,Parr2,xx2);
p1=plot(xx2,yy2);  % Plot Spline
p1.Color = 'k'; %set colour of dome to black
p1.LineStyle = '--'; %make dome curve as broken line
hold off;
