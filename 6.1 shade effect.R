#This's a part to calculate shadow area due to the tank wall, it's not in Rennie, 2017
wall.h<-Htank-M.depth                              # the wall height above manure surface, m, double
cot.alpha<-(1-sin.alpha^2)^(1/2)/sin.alpha         # angle, vector of numbers(288)
cos.theta<-(wall.h*cot.alpha/2)/ri                 # the angle in the circle-circle intersection,vector of numbers (288)
deg.theta<-acos(cos.theta)                         # vector of numbers (288) 
Intersection.h<-ri*(1-cos.theta^2)^(1/2)           # the height (m) of triangle in the circle-circle intersection ,vector of numbers (288) 
shadow<-pi*ri^2-(4*pi*ri^2*deg.theta/(2*pi)
                 -4*(wall.h*cot.alpha)/2*Intersection.h/2)  # shadow area, m2 ,vector of numbers (288) 
light.d<-1-(shadow/Au)                             # the percentage that sunlight on the surface, between 0-1
light.d[is.nan(light.d)]<-1
###End for shadow calculation
