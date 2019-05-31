% currents in mA/cm^2; conductance in S/cm^2, time in msec; V in mV; ion
% concentration in mM=(0.001 mol/liter); diameter, length in mu meter;

function dy = soma_bursting_ODE(t,y)

global gerg gcalbar gl glca Icap_max
global glna girk gnabar_s gkhhbar_s gkabar_s g_H Istim fr
global gkcabar fsca steadystate
%-----------------------------------------------------------initialization
dy=zeros(13,1);
Vs=y(1); Cai_s=y(2);
M_s=y(3); H_s=y(4); N_s=y(5); P_s=y(6); Q1_s=y(7);
Dl_s=y(8); Hs_s=y(9); 
o_erg_s=y(10); 
i_erg_s=y(11);
m_H_s=y(12); Q2_s=y(13);

%-------------------------------------------------------------- parameters 
blockacurrent = 1;
Ek = -90.0;

%/* other parameters */
ALPHA=0; PBAR=0.0225; gclamp=0*0.05; Vset=40; 
ggabaa_s=0.0e-6;ggabaa_d=0.0e-6; ggabaa_p=0.0e-6; Pnmda=0.0e-6*10.503/0.8965; ratio = 1.0000;
gampa = 0.0e-6*3.626/0.8965; gampak=0.0e-6*3.626/0.8965; 
cmd=1.0; cmp=1.0; cms=1.0; R=8.3140; F=96520.0; T=308.15; lamda=0.75; lamdaca=0.3; 
Nao=145.0; Ko=2.5; Ki=140.0; Cao=2.0; Mgo=1.2; Ra=40.0;
Kmca=0.0005; Kmna=10.0; Kmmg=50.7; Kmn=0.0001; Kml=0.00045; hill=4.0;
dd=1.5; dp=3.0; ds=15.0; Ld=350.0; Lp=150.0; Ls=25.0; fp=2.0; fs=2.0; fd=2.0; 
q = 9.0; barinap_d=0.0092; barinap_p=0.0092; barinap_s=0.0092; baricap=0.00191; %// Maximum Ca pump current
dom=100.00; taudom=1.0; d1=0.84; d2=1.0; k1=0.18; k2=0.011; bbar=0.28; abar=0.48;
	
% ----------- Electrode current 
%Istim = 0; 
istim = 0.1*Istim/(pi*ds*Ls);

%----------- Leak currents: 
ilca_s = glca*(Vs - 50);
ilk_s_irk = girk*(Vs +90)./(1+exp((Vs+45)/20));
ilna_s = glna*(Vs - 60);
il_s=gl*(Vs+65);
il_s = ilna_s + ilk_s_irk + il_s;

%--------- Gating variables of sodium (m, h, hs) + DR (n) + KA (p, q)
hsinf_s = 1/(1 + exp((Vs+54.8)/1.57)); %Kun's model
ninf_s = 1/(1 + exp(-(Vs + 25.0)/12.0));
minf_s = 1/(1 + exp(-(Vs + 30.09)/13.2));

hinf_s = 1/(1 + exp(((Vs + 54)/12.8)));
pinf_s = 1/(1 + exp( -(Vs+35.1)/13.4));
qinf_s = 1/(1 + exp( (Vs+80)/6));

    am=-1.9565e1;bm=-5.0542e-1;cm=7.9992e-1;dm=3.0212;em=-7.4630e-3;taum_shift=0.01;
mtau_s = taum_shift + 1/(cm*(am+bm*Vs)/(exp(am+bm*Vs)-1) + dm*exp(em*Vs));
    ah=5.0754e-4;bh=6.3213e-2;ch=9.7529;dh=1.3442e-1;tauh_shift=0.4;
htau_s = 1/(ah*exp(-bh*Vs)+ch*exp(dh*Vs)) + tauh_shift;
hstau_s = 20+580/(1+exp(Vs/1)); 
%   sc=20;loc=-50;sh=-0.02;sl=350;
%   x=1.0-(sh*(Vs-loc));if sh~=0 y=(-1.0/sh)*(log(x)); else y=Vs-loc; end
%if x>=0 ntau_s=1+sc*exp(-y*y/sl); else ntau_s=1; end
x=[-61.1253    4.4429  -36.8869   -9.7083    0.0052   27.2598    0.8876];
ntau_s=x(6)*(1./(1 + exp(-(Vs-x(1))/x(2)))).*(1./(1 + exp(-(Vs-x(3))/x(4)))+x(5))+x(7);
%ptau_s=1.77*exp(-0.051*Vs);  if Vs<=-65 ptau_s=63; end
x=[ -71.5402   26.0594  -62.5026   -6.5199   -0.5108   95.5813   48.2438];
ptau_s=x(6)*(1/(1 + exp(-(Vs-x(1))/x(2))))*(1/(1 + exp(-(Vs-x(3))/x(4)))+x(5)) +x(7);
q1tau_s=6.1*exp(0.015*Vs);
  %if Vs>=-40 q2tau_s=26.4515*exp(0.0227*Vs); else q2tau_s=60./(1 + exp( (Vs+50)/6)); end
   x=[294.0087   55.8321  -52.5933   -4.9104   -5.2348   84.8594   35.3239];
q2tau_s= x(1)+x(2)*(1./(1 + exp(-(Vs-x(3))/x(4)))+x(5)).*(1./(1 + exp((Vs-x(6))/x(7))));

%---------L type calcium current + cacilum pump + SK current
if Cai_s > 0.0
	cinf_s = 1.0/(1.0 + (0.00019/Cai_s)^4);
	icap_s = Icap_max/(1 + (0.0005/Cai_s)^1);
else
	cinf_s = 0;
	icap_s = 0;
end

dlinf_s = 1/(1 + exp(-(Vs + 45.0)/7.5));
%dltau_s = 18.0*exp(-((Vs + 70.0)*(Vs + 70.0))/625.0) + 0.30;
% time constant of calcium from Putzier 2009
a1=-0.020876;b1=39.726;c1=4.711;a2=0.19444;b2=15.338;c2=224.21;
linoid=a1*(Vs+b1)/(exp(-(Vs+b1)/c1)-1);exponential=a2*exp(-(Vs+b2)/c2);
dltau_s = 1/(linoid+exponential);
ikca_s = gkcabar*cinf_s*(Vs +90);

%----------H current
m_Hinf=1/(1+exp((Vs+77.6)/17.317));
tau_H_s = 26.21+3136/(1+exp(-(Vs+22.686)/29.597));

%---------ERG potasssium 
 alpha_o_s = 0.0036*exp(0.0759*Vs);     	beta_o_s  = 1.2523e-5*exp(-0.0671*Vs);
 alpha_i_s = 91.11*exp(0.1189*Vs);      	beta_i_s  = 12.6*exp(0.0733*Vs);
      
    ical_s = gcalbar*Dl_s*(Vs - 50);
    ikhh_s = gkhhbar_s*N_s^3*(Vs +90);
    ika_s = gkabar_s*P_s*(Q1_s/2+Q2_s/2)*(Vs +90);
    ina_s = gnabar_s*M_s*M_s*M_s*H_s*(fr+(1-fr)*Hs_s)*(Vs - 60); 
    %ina_s = gnabar_s*M_s*M_s*M_s*H_s*(fr+(1-fr)*Hs)*(Vs - 60); %bifurcation Hs constant
    I_H_s = g_H*m_H_s^2*(Vs+29);
    %o_erg=beta_i_s*OIerg/(alpha_i_s+beta_i_s); %bifurcation OIerg constant
    %Ierg_s = gerg*o_erg*(Vs - Ek);
    Ierg_s = gerg*o_erg_s*(Vs - Ek); 

%-------------summary of currents
icalcium_s = ilca_s + ical_s + icap_s;
isoma = ika_s + ikhh_s + ina_s + il_s + ilca_s + ical_s + ikca_s + Ierg_s + I_H_s - istim ;
		
%-------------------------------------------------------------- ODE governing dynamics of model*/
dy(1) =  -1000*isoma/cms; %Ith(ydot, V_s)
dy(2) =  -2*fsca*icalcium_s/(ds*0.0001*F); %Ith(ydot, Cai_s)
dy(3:7) = [ (minf_s - M_s)/mtau_s;  %Ith(ydot, M_s) =
              (hinf_s - H_s)/htau_s;  %Ith(ydot, H_s) = 
              (ninf_s - N_s)/ntau_s*1.2;  %Ith(ydot, N_s) =
              (pinf_s - P_s)/ptau_s;  %Ith(ydot, P_s) =
              (qinf_s - Q1_s)/q1tau_s];  %Ith(ydot, Q_s) =
dy(8) = (dlinf_s - Dl_s)/dltau_s; %Ith(ydot, Dl_s) =
dy(9) = (hsinf_s - Hs_s)/hstau_s; %Ith(ydot, Hs_s)  %STEP2
dy(10) =  alpha_o_s*(1-o_erg_s-i_erg_s) + beta_i_s*i_erg_s - (alpha_i_s+beta_o_s)*o_erg_s;% d o_erg /dt
dy(11) =  alpha_i_s*o_erg_s - beta_i_s*i_erg_s; 
dy(12) = (1/(1+exp((Vs+77.6)/17.317))-m_H_s)/tau_H_s;  % I_H unspecfic hyperpolarizing current
dy(13)=  (qinf_s - Q2_s)/q2tau_s;
        
