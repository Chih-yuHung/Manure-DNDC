#This file calculate heat transfer and it iterate everyday

#Environmental input
ReL<-wind*ri/Vair                                   #Reynold's number, unitless,F42
Nu<-ifelse(ReL<5*10^5,0.453*(ReL^(0.5))*(Pr^(1/3))  #Nusselt Number,F43
           ,(0.037*(ReL^(4/5))-871)*(Pr^(1/3)))    
hcv.ms<-(Nu*ka)/ri                                  #Heat transfer coefficient,F44       

#Radiative heat transfer
declination.s<-23.45*sin((2*pi*(284+T.day)/365))             # seasonal declination(degree),F101:KG101
sin.alpha<-pmax((cos(deg2rad(L))*cos(deg2rad(declination.s))
          *cos(deg2rad(H))+sin(deg2rad(L))
          *sin(deg2rad(declination.s))),0)                   # sunlight degree, F102:KG102

#This's a part to calculate shadow area due to the tank wall, it's not in Rennie, 2017
source("6.1 shade effect.R",echo=F)
###End for shadow calculation

m<-ifelse(sin.alpha>0,Pa/(101325*sin.alpha),0)       # Optical air mass number, #F103-KG103
Sb<-ifelse(sin.alpha>0, Eb*(tau^m)*sin.alpha,0)      # solar bean radiation (W/m2),F104-KG104
Sd<-ifelse(sin.alpha>0,0.3*(1-tau^m)*Eb*sin.alpha,0) # Diffusive radiation (w/m2),F105-KG105
Sr.total<-sum(Sb,Sd)                                 # F322, Total solar radiation
q.net.rad<-alpha.s*((Sb+Sd)/Sr.total)*((SR*1000*1000)/T.delta) #Net solar radiation, F106:KG106
q.net.rad<-q.net.rad*light.d                         #apply shade coefficient  

#Relative humidity from measured data
#Rh estimated based on RH6 and RH15 with T.hour,F110:KG110
Rh<-c(1:288)  
Rh[1:71]<--((RH6-RH15)/2)*cos((-9-T.hour[1:71])*pi/15)+((RH6+RH15)/2)
Rh[72:180]<-((RH6-RH15)/2)*cos((6-T.hour[72:180])*pi/9)+((RH6+RH15)/2)
Rh[181:288]<-((RH6-RH15)/2)*cos((6-T.hour[181:288])*pi/9)+((RH6+RH15)/2)


#Estimate air temp.
sunrise<-T.hour[which(sin.alpha>0,arr.ind=TRUE)[1]] #determine sunrise time,F119:KG119

#hr, sunrise reference, sunrise= 0, an hour before sunrise =23,F120:KG120
sunrise.ref<-ifelse(T.hour<sunrise,24+T.hour-sunrise
       ,T.hour-sunrise)
#x, in Schaub,1991,F121:KG120
x<-ifelse(sunrise.ref>=0&sunrise.ref<=14-sunrise
         ,x<-cos(sunrise.ref*pi/(14-sunrise))
         ,x<-cos((sunrise.ref-((14-sunrise)+1))*(pi/(23-(14-sunrise)))))

#phase, to determine the phase of sunrise, 1.before sunrise, 2, 
#after sunrise before sunset, 3, after sunset
T.air<-ifelse(T.hour<sunrise,(AirTmax0-AirTmin1)/2*x+((AirTmax0+AirTmin1)/2)
       ,ifelse(sunrise.ref>=0 & sunrise.ref<=(14-sunrise),
               (-((AirTmax1-AirTmin1)/2)*x+(AirTmax1+AirTmin1)/2)
               ,(((AirTmax1-AirTmin2)/2)*x+(AirTmax1+AirTmin2)/2)))

T.air.K<-T.air+273.15

#WVPD,F111:KG111
WVPD<-Teten.H2Oa*exp((Teten.H2Ob*T.air)/(Teten.H2Oc+T.air))*(1-Rh/100)
#Evaporation per second (kg/s), F112:KG112
E<-rho.w*(WVPD)*wind.f/(24*3600*1000)*Au
Evap.depth.d<-sum(E*T.delta)/rho.w/Au #Incorporate daily evaporation, depth together
#air emissivity, F116:KG116
#in excel, the last part was (F95+273.15)^(1/7), F95 is T.hour
#I replace F95+273.15 to T.air.K
e.ac<-1.72*(((Teten.H2Oa*exp((Teten.H2Ob*(T.air))/(Teten.H2Oc+(T.air)))*(Rh/100))/T.air.K)^(1/7))

#cloud-corrected air emissivity, F117:KG117 
e.a<-(1-0.84*cc)*e.ac+0.84*cc

#Soil temperature, 300 cells, 2.995m
S.Temp[300,]<-annualT.K   #The deepest soil temperature was assumed to be annual T
#Need delta.z[30] source manure volume,F9100
delta.depth<-delta.z[30]/2+dep.s/2#F9103
soil.c<-T.delta/(den.s*(Au*dep.s))#constant of soil, F9100

#Thermal conductivity/specific heat correction, F167:KG196
#Manure temperature calculation,F132:KG162 
#soil temperature cacultation, E9107:KH9107
#The process is to calculate 5 mins thermal conductivity 
#and then 5 mins Manure temperature from soil temp. and pre. manure temp
#soil temp was from pre. soil temp and manure temp.
#use the Manure temp in previous 5 mins and calculate thermal conductivity
Cp<-c(1:288)#Specific heat of manure, two values, frozen or liquid F108:KG108
Ts<-c(1:288)#manure temp at degree C, F109:KG109, not sure the purpose here
Ts[1]<-ini.M.Temp[1]-273.15
T.conductivity<-matrix(ncol=288,nrow=30) # Conductivity, F166:KG196
delta.T.evap<-c(1:288) #delta T-evap, F114:KG114
delta.T.radevap<-c(1:288)   #delta T-rad+evap, F115:KG115


for (j in 1:288) {
  if (j ==1) {
    T.conductivity[,j]<-ifelse(ini.M.Temp>=273.15|ini.M.Temp<272.15,k.m/C.pm,k.m/C.pm.fusion)
    Cp[j]<-ifelse(ini.M.Temp[j]>=273.15|ini.M.Temp[j]<272.15,C.pm,C.pm.fusion)
    delta.T.evap[j]<-(-(E[j]*Lambda*T.delta))/(Cp[j]*(rho.m*Au*delta.z[1]))
    delta.T.radevap[j]<-(q.net.rad[j]*Au*T.delta-(e.sigma*Au*T.delta*(((epsilon*(ini.M.Temp[j])^4))-(e.a[j]*T.air.K[j]^4))))/(rho.m*Cp[j]*Au*delta.z[1])+delta.T.evap[j]
    M.Temp[1,j]<-(ini.M.Temp[j]+time.weight[1]*T.conductivity[1,j]*(Au*hcv.ms*T.air.K[j]+Au/delta.zd[1]*ini.M.Temp[j+1]))/(1+time.weight[1]*T.conductivity[1,j]*(Au*hcv.ms+Au/delta.zd[1]))+delta.T.radevap[j]
    M.Temp[2:29,j]<-(ini.M.Temp[2:29]+time.weight[2:29]*T.conductivity[2:29,j]*(Au/delta.zu[2:29]*ini.M.Temp[1:28]+Au/delta.zd[2:29]*ini.M.Temp[3:30]))/(1+time.weight[2:29]*T.conductivity[2:29,j]*(Au/delta.zu[2:29]+Au/delta.zd[2:29]))
    M.Temp[30,j]<-(ini.M.Temp[30]+time.weight[30]*T.conductivity[30,j]*(Au/delta.zu[30]*ini.M.Temp[29]+Au/delta.depth*ini.S.Temp[1]))/(1+time.weight[30]*T.conductivity[30,j]*(Au/delta.zu[30]+Au/delta.depth))
    S.Temp[1,j]<-(ini.S.Temp[1]+soil.c*ks.cp*(Au/delta.depth*ini.M.Temp[30]+Au/dep.s*ini.S.Temp[2]))/(1+soil.c*ks.cp*(Au/delta.depth+Au/dep.s))
    S.Temp[2:299,j]<-(ini.S.Temp[2:299]+soil.c*ks.cp*(Au/dep.s*ini.S.Temp[1:298]+Au/dep.s*ini.S.Temp[3:300]))/(1+soil.c*ks.cp*(Au/dep.s+Au/dep.s))
  } else {
    T.conductivity[,j]<-ifelse(M.Temp[,j-1]>=273.15|M.Temp[,j-1]<272.15,k.m/C.pm,k.m/C.pm.fusion)  
    Cp[j]<-ifelse(M.Temp[1,j-1]>=273.15|M.Temp[1,j-1]<272.15,C.pm,C.pm.fusion)
    delta.T.evap[j]<-(-(E[j]*Lambda*T.delta)/(Cp[j]*(rho.m*Au*delta.z[1])))
    delta.T.radevap[j]<-(q.net.rad[j]*Au*T.delta-(e.sigma*Au*T.delta*((epsilon*(M.Temp[1,j-1])^4)-(e.a[j]*T.air.K[j]^4))))/(rho.m*Cp[j]*Au*delta.z[1])+delta.T.evap[j]
    M.Temp[1,j]<-(M.Temp[1,j-1]+time.weight[1]*T.conductivity[1,j]*(Au*hcv.ms*T.air.K[j]+Au/delta.zd[1]*M.Temp[2,j-1]))/(1+time.weight[1]*T.conductivity[1,j]*(Au*hcv.ms+Au/delta.zd[1]))+delta.T.radevap[j]
    M.Temp[2:29,j]<-(M.Temp[2:29,j-1]+time.weight[2:29]*T.conductivity[2:29,j]*(Au/delta.zu[2:29]*M.Temp[1:28,j-1]+Au/delta.zd[2:29]*M.Temp[3:30,j-1]))/(1+time.weight[2:29]*T.conductivity[2:29,j]*(Au/delta.zu[2:29]+Au/delta.zd[2:29]))
    M.Temp[30,j]<-(M.Temp[30,j-1]+time.weight[30]*T.conductivity[30,j]*(Au/delta.zu[30]*M.Temp[29,j-1]+Au/delta.depth*S.Temp[1,j-1]))/(1+time.weight[30]*T.conductivity[30,j]*(Au/delta.zu[30]+Au/delta.depth))
    S.Temp[1,j]<-(S.Temp[1,j-1]+soil.c*ks.cp*(Au/delta.depth*M.Temp[30,j-1]+Au/dep.s*S.Temp[2,j-1]))/(1+soil.c*ks.cp*(Au/delta.depth+Au/dep.s))
    S.Temp[2:299,j]<-(S.Temp[2:299,j-1]+soil.c*ks.cp*(Au/dep.s*S.Temp[1:298,j-1]+Au/dep.s*S.Temp[3:300,j-1]))/(1+soil.c*ks.cp*(Au/dep.s+Au/dep.s))
  }
}
Ts[2:288]<-M.Temp[2:288] # don't know the purpose

#Temp and depth adjustment, F200:Q238
#Current enthalpy, J209:J238
Enthalpy.c<-ifelse(M.Temp[,288]<272.15,M.Temp[,288]*rho.m*M.volume*C.pm/10^6
                   ,ifelse(M.Temp[,288]>=273.15,(272.15*rho.m*M.volume*C.pm+rho.m*M.volume*C.pm.fusion+(M.Temp[,288]-273.15)*rho.m*M.volume*C.pm)/10^6
                           ,(272.15*rho.m*M.volume*C.pm+(M.Temp[,288]-272.15)*rho.m*M.volume*C.pm.fusion)/10^6))

#New enthalpy
depthchange.d<-M.storage/Au/365-Evap.depth.d          #L34
In.M.volume<-Au*depthchange.d                         #L41
In.M.volume.p<-Au*precip.d                            #precip volume
depth.factor<-depthchange.d/M.depth                   #N204
depth.factor.p<-precip.d/M.depth                      #change results from precipitation
delta.z.new<-delta.z*(1+depth.factor)                 #L209:238
M.volume.new<-delta.z.new*Au                          #new manure volume,M209:M238
delta.z.new.p<-delta.z*depth.factor.p                   #depth factor for precipitations
M.volume.new.p<-delta.z.new.p*Au                      #new precipitation volume
#incoming Manure temp
In.M.temp<-Avg.Barn.temp+Barn.temp.amp*sin(2*pi/365*T.day+Temp.cost) #Incoming manure temp, L49,L39
#Enthalpy after manure added, N209:N238
Enthalpy.c.new<-Enthalpy.c+
  (M.volume.new-M.volume)*rho.m*((In.M.temp+273.15)*C.pm+C.pm.fusion)/1000000+
  (M.volume.new.p)*rho.m*((Tmean+273.15)*C.pm+C.pm.fusion)/1000000

Enthalpy.V<-Enthalpy.c.new/(M.volume.new+M.volume.new.p)  #Enthalpy/V, O209:O238

#Final temp after depth adjustment,Q209:Q238
#This is actually the manure temperature after manure addition and we used this 
#to be the new initial manure temp for the next day
#not the manure temp at the end of the day!
Final.M.Temp<-ifelse(Enthalpy.V<E.272,272.15*Enthalpy.V/E.272,
                     ifelse(Enthalpy.V>=E.273,273.15+(Enthalpy.V-E.273)*10^6/(C.pm*rho.m)
                            ,272.15+(Enthalpy.V-E.272)/fusion))

if (mean(Final.M.Temp)>=(50+273.15)) {
  cat("Manure temperature too high to be true")
  break
}



