library(REdaS); library(xlsx); library(beepr) ;library(dplyr); library(imputeTS)
#Set location, initial date and end time; date origin in R, 1970-1-1
Location<-"Arlanda"
start.date<-"2020-6-18"  
end.date<-"2023-6-17" 
#removal dates in multiple years, can be different dates
removal.dates<-as.numeric(as.Date(c("2020-10-12","2021-4-25"
                                    ,"2020-10-12","2021-4-25"
                                    ,"2021-10-12","2022-4-25"
                                    ,"2022-10-12","2023-4-25"),by="days"))
removal.end<-as.numeric(as.Date(c("2020-10-22","2021-5-27"
                                  ,"2020-10-22","2021-5-27"
                                  ,"2021-10-22","2022-5-27"
                                  ,"2022-10-22","2023-5-27"),by="days"))
#The duration of removals
removal.day<-unique(removal.end-removal.dates)
removal.duration<-list()
    for (i in 3:8){
        removal.duration[[i-2]]<-c(removal.dates[i]:removal.end[i])}

#Input environmental data
Envir.daily<-read.csv("input/daily env input_Arlanda_June18.csv",header=T)
#To produce an extra year for balance soil temperature
Envir.daily<-Envir.daily[c(1:365,1:1095),]
#to know how many days we have for the loop
d.length<-nrow(Envir.daily)
#initial manure temp
ini.M.Temp<-read.csv("input/Initial M temp.csv",header=T) 
ini.M.Temp<-ini.M.Temp[,"Initial.Temp"] #change to vector

#Read parameters
source("3.Parameters.R",echo=F)           #Parameters we can change
source("4.Constants_sweden.R",echo = F)   #Constants no need to change

#Store daily manure.temp in the 30 layers
manure.temp<-c()
manure.depth<-c()

#Start the simulation here. 
starttime<-Sys.time()
for (i in 1:d.length){
#To read daily data from environment input.
print(paste("Date ID =",Output$`Date ID`[i],"DOY=",Output$DOY[i])) 
T.day<-Output$DOY[i]
AirTmax0<-Envir.daily$AirTmax0[i]
AirTmax1<-Envir.daily$AirTmax1[i]
AirTmin1<-Envir.daily$AirTmin1[i]
AirTmin2<-Envir.daily$AirTmin2[i]
SR<-Envir.daily$SR[i]
wind<-Envir.daily$wind[i]              #m/s
wind.v<-Envir.daily$wind[i]            #daily wind speed at 2m, m/s
wind.f<-(2.36+1.67*wind.v)*Au^(-0.05)  #wind function (for evaporation), mm/d/kPa, G16
cc<-min(Envir.daily$cloud[i],1)    #cloud clover
precip.d<-Envir.daily$precip[i]/1000
source("3.1.Alpha.s_adjustment.R",echo=F) #adjust the albedo in winter
RH6<-Envir.daily$RH.6[i]
RH15<-Envir.daily$RH.15[i]
#Manure depth adjustment
#Reset M.depth and ini.M.temp after soil temperature stabilization
if (sum(i==365|i==730|i==1095|i==1460)==1){
M.depth<-0.599032
Zmmax<-M.depth
}

# If current date = removal dates then calculate the depth change and cat removal date
if(sum(removal.dates[c(3,5,7)] == Output$`Date ID`[i]) == 1) {
    removal.depth.d<-(M.depth-removal.depth[1])/removal.day[1] #the difference between the peak and the minimum 
    cat(paste("manure removal date =",i))
    }

if(sum(removal.dates[c(4,6,8)] == Output$`Date ID`[i]) == 1) {
   removal.depth.d<-(M.depth-removal.depth[2])/removal.day[2]    
    cat(paste("manure removal date =",i))
}


#To calculate manure volume
source("5.Manure volume_sweden.R",echo=F)
#To calculate solar radiation, soil temp, and manure temp at 300 sec steps.
source("6.Solar radiation and soil temp_sweden_shade.R",echo=F)
#To calculate final hourly temp
source("7.hourly temp_sweden.R",echo=F)

#To obtain temp and depth at the end of the day
#I need manure temperature in different depth every day for the last year
#It would be a table with DOY, manure.depth, and temp
if(is.element(i,tail(1:d.length,n=365))){
    manure.temp<-c(manure.temp,M.temp.d)
    aa<-seq(0,M.depth,length.out=30)
    manure.depth<-c(manure.depth,aa)
}

#Write the results, only write after first year
Output[i,6]<-Avg.M.temp.d     #Daily manure temperature (C), before depth adjustment
print(paste("Sequence",i,"And Manure temp",Avg.M.temp.d))
print(paste("alpha.s=",alpha.s))
Output[i,7]<-M.depth*100      #Daily manure depth (cm)
Output[i,8]<-M.volume.current #Daily manure volume(m3)
Output[i,9]<-Evap.depth.d*100 #Daily Evaporation (cm)
Output[i,10]<-precip.d*100    #Daily Precipitation (cm)
Output[i,11]<-sum(q.net.rad)  #Net solar radiation

#write the result for manure temperature at 0.5m
source("8.Manure temperature at 0.5.R",echo=F)
        
#daily changing depth of manure for next day
if (Output$`Date ID`[i] %in% removal.duration[[1]]|
    Output$`Date ID`[i] %in% removal.duration[[3]]|
    Output$`Date ID`[i] %in% removal.duration[[5]]) {
    M.depth<-M.depth-removal.depth.d
    Zmmax<-M.depth
} else if(Output$`Date ID`[i] %in% removal.duration[[2]]|
          Output$`Date ID`[i] %in% removal.duration[[4]]|
          Output$`Date ID`[i] %in% removal.duration[[6]]) {
    M.depth<-M.depth-removal.depth.d
    Zmmax<-M.depth
} else {
M.depth<-M.depth+depthchange.d
Zmmax<-M.depth    
}
#Save the new temperatures
ini.M.Temp<-Final.M.Temp
#This the soil temperature(to depth 2.995m)
ini.S.Temp<-S.Temp[,288]
}
endtime<-Sys.time()
#Output<-Output[366:d.length,] # first year is for soil temp stabilization
totaltime<-endtime-starttime
totaltime

# Summary for the results, part 2. Only output after all simulation is done. 
#The data of manure tank, i.e output L1:M18
Output.tank<-data.frame(matrix(ncol=2,nrow = 18))
Output.tank[1:18,1]<-c("Location","SurfaceArea(m2)","Starting.Depth(m)"
                       ,"Starting.Volume.m3","Total.Solids(%)",""
                       ,"Tank.Storage","Total.Tank.Volume(m3)"
                       ,"Yearly Maximum Storage Volume.m3","Yearly.Rain.Volume.m3"
                       ,"Yearly.Manure.Storage.Volume.m3","Tank.Diameter.m",""
                       ,"Average.Tm.C","Max.Tm.C","Min.Tm.C","","Max.d.cm")
Output.tank[1,2]<-Location
Output.tank[2:3,2]<-c(Au,M.depth)        #area and initial depth, m2 and m
Output.tank[4,2]<-as.numeric(Output.tank[2,2])*as.numeric(Output.tank[3,2]) #starting volume.
Output.tank[5,2]<-Total.solid                 #%
Output.tank[c(6,7,13,17),2]<-""                      #blank
Output.tank[8,2]<-Tank.v                      #Total Tank volume, m3
Output.tank[9,2]<-Max.storage                 #yearly manure storage,m3
Output.tank[10,2]<-Rain.v                     #yearly rain volume in storage,m3
Output.tank[11,2]<-M.storage                  #yearly Manure storage,m3
Output.tank[12,2]<-ri*2                       #Diameter of tank,m
Output.tank[14,2]<-mean(Output$Temperature.C) #Avg. manure Temperature for the estimate years
Output.tank[15,2]<-max(Output$Temperature.C)  #Max manure Temperature
Output.tank[16,2]<-min(Output$Temperature.C)  #Min manure Temperature
Output.tank[18,2]<-max(Output$Depth.cm)       #Maximum Manure Depth

#tabulate a daily manure temp and depth 
DOY<-rep(c(167:365,1:166),each=30)
daily.data<-as.data.frame(cbind(DOY,manure.depth,manure.temp))

#output to an excel file
write.csv(Output,paste("output/",Location,Sys.Date(),".csv",sep=""),row.names = FALSE)

