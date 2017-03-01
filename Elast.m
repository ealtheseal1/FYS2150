%Data til Ovelsen Elastisitet

load kjoring1

m_teoretisk = [0.0 0.5 1.0 1.5 2.0 2.5 3.0 3.5]; %kg
utsving = [9.99 9.27 8.55 7.84 7.12 6.42 5.71 5.01]; %mm
h2 = [9.99-0.04 9.27-0.04 8.55-0.03 7.84-0.05 7.12-0.05 6.42-0.07 5.71-0.09 5.01-0.1];
h3 = [9.95 9.24 8.51 7.79 7.10 6.35 5.62 4.90]; %mm
% Forskyvning i h(0) til 9.95 mm pga staven vris seg fordi ulik masse
% distribusjon på vektskål

deltah = [0.04 0.04 0.03 0.03 0.05 0.05 0.07 0.09 1.0];
deltah3 = 0.2; %mm



h = fittype('alpha*x+beta');
[fit1,gof,fitinfo] = fit(m_teoretisk(:),h3(:),h,'StartPoint',[0 0]);
residuals = fitinfo.residuals;
I = abs( residuals) > 1.5 * std( residuals );
outliers = excludedata(m_teoretisk(:),h3(:),'indices',I);

fit2 = fit(m_teoretisk(:),h3(:),h,'StartPoint',[0 0],...
           'Exclude',outliers)

x = m_teoretisk(:);
y = h3(:);
       
p=polyfit(x,y,1);
%Tilpasningsparametre m og c i y=m*x+c:
m=p(1)
c=p(2)
yline = polyval(p,x);
       
D = sum(x.^2) - 1/length(x) *(sum(x)).^2;
d = y - m*x - c;

deltaMkvadr = 1/(length(y)-2) * sum(d.^2)./D;

deltaM = sqrt(deltaMkvadr)

deltaCkvadr = 1/(length(x)-2) * (D/length(x) + mean(x))* sum(d.^2)/D;

deltaC = sqrt(deltaCkvadr)

figure(1)
plot(m_teoretisk,utsving ,'o')
hold on
plot(m_teoretisk, h3, '+')
hold on
plot(x,yline)
xlabel('masse [kg]')
ylabel('h(m) , [mm]')
legend()

massestav = [2.485 2.487]; %kg

%Øvelse 2

diameter = [15.97 16.06 16.00 15.98 16.02]*1e-3; %m
lengde = 1.335; %m
deltaL = (0.5*1e-3)+sqrt(6)*0.5+(1.4*1e-3); %lengde mellom knivbladene meter
deltaSTOREL = (0.1*1e-3)+sqrt(7)*0.5 +(1.4*1e-3); %hele lengden [m]

Eteoretisk = (4*(lengde.^3)*9.81)./(3.0*pi*abs(m*1.0e-3)*((mean(diameter)).^4))

sE = Eteoretisk*sqrt(((2*std(diameter))/mean(diameter))^2 + ((2*std(fw))/mean(fw))^2 + (deltaL/lengde)^2 + ((std(massestav))/mean(massestav))^2)

%Øvelse 4



grunntonefrekvens = sqrt((Eteoretisk*pi*(mean(diameter))^2)./(16*(mean(massestav))*lengde))
maaltfrekvens = mean(fw);

deltaF = sqrt( maaltfrekvens^2 * ((sE./Eteoretisk)^2 - ((2*std(diameter))/(mean(diameter)))^2 - (deltaL/lengde)^2 - ((std(massestav))/mean(massestav))^2)) *0.5


