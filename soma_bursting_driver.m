clear all

%initial condition 
    if exist('initial.dat', 'file')% File exists.  
       load initial.dat
       Y=initial;
else
Y =[0*ones(1,1) 1e-4*ones(1,1) 0.01*ones(1,11)]';
end
%parameters
global gerg gcalbar gl glca  Icap_max fsca gkcabar 
global glna girk gnabar_s gkhhbar_s gkabar_s g_H Istim fr
gcalbar = 139e-6; 
gl=280e-6;
gnabar_s = 6000e-6;
fr=0;
gkabar_s = 1680e-6;
fsca=0.018;%0.036
Icap_max=0.011;%0.011
glca = 2.45e-6;  
glna = 0;girk = 0;
gerg=130e-6; 
g_H=78e-6; 
gkhhbar_s = 1117e-6;
gkcabar=0*70e-6;
Istim=0;

% %================ running ODE
options=odeset('RelTol',1e-9,'Stats','on');
%[t,y]=ode45(@soma_bursting_ODE,[0 30000],Y);
[t,y]=ode15s(@soma_bursting_ODE,[0 20000],Y);
figure(5); hold on
%plot(t(1:end)*0.001,y(1:end,1),'k')
    
plot(0.001*t,y(:,1),'k')


data=[0.001*t,y(:,1)];
%data=[0.001*t,y(:,10)+y(:,11)];
endvalue=y(end,:)';
save vs.dat data -ascii
save state.dat endvalue -ascii
