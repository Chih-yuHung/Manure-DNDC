#Manure volume change, J60:O89
zd[1]<-(Zmmax/2)/sum((1+grid.c)^c(0:14))
delta.z[1:15]<-((1+grid.c)^(0:14))*zd[1]
delta.z[16:30]<-((1+grid.c)^(14:0))*zd[1]
for (j in 2:30){
  zd[j]<-zd[j-1]+delta.z[j]
}
for (j in 2:30){
  zu[j]<-zd[j-1]
}
zp<-(zu+zd)/2
for (j in 1:29){
  delta.zd[j]<-abs(zp[j]-zp[j+1])
  delta.zd[30]<-0
}
delta.zu<-c(0,delta.zd[1:29])
M.volume<-Au*delta.z                   # calculate manure volume, m3, T60-T89, 
M.volume.current<-sum(M.volume)        # current volume,L42
time.weight<-T.delta/(rho.m*M.volume)  # V60:V89  

#Temp
VT<-ini.M.Temp*M.volume                    # Initial temp x volume
Avg.M.temp<-sum(VT)/sum(M.volume)-273.15   # avg. manure temp, F86 , degree C
suf.M.temp<-ini.M.Temp[1]-273.15           # surface manure tmep, F87,  degree C
mid.M.temp<-mean(ini.M.Temp[15:16])-273.15 # middle manure tmepe, F88 , degree C
bot.M.temp<-ini.M.Temp[30]-273.15          # Bottom manure temp, F89 , degree C  

