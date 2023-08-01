function simulate_global_EO_Vcmax25(Year,MERRA_2_file,leaf_absorptance_file,save_file)

%% Reading climate data
T=ncread(MERRA_2_file,'T2M');
T(T<=0)=NaN;
Tg=nanmean(T,3);
Ps=ncread(MERRA_2_file,'Ps');
I=ncread(MERRA_2_file,'I');
Dg=ncread(MERRA_2_file,'Dg');


CO2_NASA=[340.120000000000;341.480000000000;343.150000000000;344.870000000000;346.350000000000;347.610000000000;349.310000000000;351.690000000000;353.200000000000;354.450000000000;355.700000000000;356.540000000000;357.210000000000;358.960000000000;360.970000000000;362.740000000000;363.880000000000;366.840000000000;368.540000000000;369.710000000000;371.320000000000;373.450000000000;375.980000000000;377.700000000000;379.980000000000;382.090000000000;384.020000000000;385.830000000000;387.640000000000;390.100000000000;391.850000000000;394.060000000000;396.740000000000;398.810000000000;401.010000000000;404.410000000000;406.760000000000;408.720000000000;411.660000000000;414.240000000000];
CO2_con=CO2_NASA(Year-1980);


%% Reading leaf absorptance data
alpha=ncread(leaf_absorptance_file,'fleaf');

%% atmospheric CO2 partial pressure
Ca=CO2_con.*Ps./1000000;
%% intercellular O2 concentration
Oi=0.21.*Ps;

%% Michaelis¨CMenten coefficients at 25¡ãC
Kc=404.9.*Ps./1000000;
Ko=278.4.*Ps./1000;
K=Kc.*(1+Oi./Ko);

%% CO2 compensation point in the absence of mitochondrial respiration at 25 degree
Gamma=42.75.*Ps./1000000;


%% the viscosity of water at 25¡ãC relative to its value at 25¡ãC
eta=1;

%% the ratio of the cost factor for photosynthesis to the cost factor for transpiration
beta=146;

%% simulate the optimal ratio of intercellular to ambient CO2
sigma=sqrt(beta.*(K+Gamma)./(1.6.*eta));
chi=Gamma./Ca+(1-Gamma./Ca).*sigma./(sigma+sqrt(Dg));

%% intercellular CO2 concentration,
Ci=Ca.*chi;


%% the fraction of incident quanta utilized in electron transport for Schemes 4 and 5
phi_Scheme4=alpha.*(0.352+0.022.*25-0.00034.*25.*25).*0.5;
phi_Scheme5=0.7455.*(0.352+0.022.*25-0.00034.*25.*25).*0.5;

%% absorbed photosynthetic photon flux density for Schemes 4 and 5
APAR_Scheme4=phi_Scheme4.*I;
APAR_Scheme5=phi_Scheme5.*I;

%% estimate rjv25 from Tg using the empirical relationship
rjv25=exp(-0.0193.*Tg+1.02);

%% simulate Vcmax25 for Schemes 4 and 5
m=(Ci-Gamma)./(Ci+2.*Gamma);
mc=(Ci-Gamma)./(Ci+K);
M=4.*mc./m;
Vcmax25_Scheme4=APAR_Scheme4.*sqrt((1./M).^2-(1./rjv25).^2);
Vcmax25_Scheme5=APAR_Scheme5.*sqrt((1./M).^2-(1./rjv25).^2);

%% save results
ncwrite(save_file,'Vcmax25_Scheme4',Vcmax25_Scheme4);
ncwrite(save_file,'Vcmax25_Scheme5',Vcmax25_Scheme5);

end