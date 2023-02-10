# Ho inserito alcuni punti in cui il programma deve essere integrato per funzionare
# correttamente e fare tutti i calcoli necessari, non ci sono errori ma solo parti evidenziate da integrare 



#Program to read data from particle counter
rm(list = ls())

library("openair")
library("maps")
library("lattice")



#set working directory, change as appropriate
setwd("C:/home/laura/Scrivania/File_prove_Lab")

file_in <- "21-11-10_989C.CSV"

#read data table, skip first 2 lines after the header
datos <- read.table(file=file_in, sep=",", na.string="null", as.is=TRUE, header=TRUE, skip = 2)

#remove first row
datos <- datos[-1,]

#in the "Time" column, substitute T  with a space and remove Z 
datos$Time <- sub("T", " ", datos$Time)
datos$Time <- sub("Z", "", datos$Time)

#read the "Time" as date, specifying the format
datos$date <- as.POSIXct(datos$Time, format = "%Y-%m-%d %H:%M:%S", "GMT")


#convert particle counts to PM [microgram/m-3]; attention to the units: particle density given in g cm-3; diameter in micron; particle number are divided by 0.1 L flow rate
bin.opc<-matrix(c(0.30, 0.50, 1.0, 2.5, 5.0, 10.0, 15.0),1,7)
D.opc<-matrix(NA,1,6)

m.opc<-matrix(NA,nrow(datos),6)
colnames(m.opc)<-c(paste0("bin",seq(1:6),"_opc"))
opc.s<-datos[,13:18]*10^-3  #convert counts/0.1 L to counts/cm-3 

for(i in 1:6){
  LB<-bin.opc[,i]
  UB<-bin.opc[,i+1]
  D.opc[,i]<-LB*(1/4*(1+(UB/LB)^2)*(1+(UB/LB)))^(1/3)
  m.opc[,i]<-1.65*opc.s[,i]*(D.opc[,i]*10^-4)^3*pi/6
}
C.opc<-m.opc*10^13






########  Parti da integrare 

#  Le prossime righe servono a calcolare il PM1, PM2.5 e PM10 e a mettere il risultato 
#  delle operazioni all'interno del dataframe datos come colonne aggiuntive


pm10 <- cbind(apply(C.opc[,1:5],1,sum,na.rm=FALSE))   # calcolo del pm10

# servono due righe analoghe ma con i giusti parametri per il calcolo del pm1 e pm2.5

datos <- cbind(datos, pm1, pm2.5, pm10)  # questa riga mette il risultato delle operazioni nel dataframe


#########



#########

# Le righe seguenti servono a creare i grafici per il confronto tra PM calcolati da voi e PM dello strumento, 
# ho lasciato solo quelle per il confronto tra PM1, vanno aggiunti PM2.5 e PM10


#compare pm retrieved from the sensor with that calculated manually, using scatterplot
scatterPlot(datos, x ="PM.1.0", y = "pm1", method="scatter", xlab = "measured PM1 (ug/m3)", ylab ="calculated PM1 (ug/m3)", linear = TRUE)

#compare pm retrieved from the sensor with that c,alculated manually, using time series plot
timePlot(datos, pollutant = c("PM.1.0", "pm1"), group = TRUE, name.pol=c("measured PM1", "calculated PM1"), ylab="PM1 (ug/m3)")


#########






#calculate indexes for bias and correlation and save as txt file

stat_PM1<-modStats(datos, mod="PM.1.0", obs ="pm1")
stat_PM2.5<-modStats(datos, mod="PM.2.5", obs ="pm2.5")
stat_PM10<-modStats(datos, mod="PM.10.0", obs ="pm10")
stat_PM <- cbind(stat_PM1, stat_PM2.5, stat_PM10)
write.table(stat_PM, file="evaluation statistics PM.txt",  dec=".", sep ="\t", col.names=T, row.names=F, append=F)

data_manual <- datos[,20:22]
data_autom <- datos[,10:12]
corr<-cor(data_manual, data_autom, method="spearman")
write.table(corr, file="spearman correlations PM.txt",  dec=".", sep ="\t", col.names=T, row.names=T, append=F)





######## Parti da integrare 

# In questa parte del programma vengono creati dei nuovi dataframe (subdata1, subdata2, subdata3),
# che contengono i vari periodi con condizioni abientali diverse da analizzare e su cui fare le medie
# Ho lasciato solo la parte di programma che crea il subdata1, fa le medie su questo periodo e calcola le distribuzioni
# dovete rifare lo stesso aggiustando i parametri per gli altri due periodi


#select first subset of data (windows close)
start.date <- as.POSIXct("2020-11-30 13:16", format = "%Y-%m-%d %H:%M", "GMT" )
end.date <- as.POSIXct("2020-11-30 13:46", format = "%Y-%m-%d %H:%M", "GMT")
subdata1 <- subset(datos, date >= start.date & date <= end.date)



# average the data for specified time period, save to new data .txt file, decimal point, tab separated
datos_first <-timeAverage(subdata1, avg.time = "30 min", statistic = "mean", vector.ws=FALSE)

write.table(datos_first, file="first period average.txt",  dec=".", sep ="\t", col.names=T, row.names=F, append=F)


#extract just particle number from the two subsets
PN1 <- datos_first[,11:16]


#convert particle number distribution into particle surface distribution, and then plot
S.opc <- pi*D.opc^2
PS1 <- S.opc*PN1


#convert particle number distribution into particle volume distribution, and then plot
V.opc <- pi/6*D.opc^3
PV1 <- V.opc*PN1


############









############

# Le righe seguenti servono a creare un grafico contenente le size distribution medie dei primi due periodi insieme, 
# se volete aggiungere il terzo periodo dovete modificarli
# anche i momenti vengono calcolati per i primi due periodi e vanno integrati per il terzo periodo





#plot average particle size distribution for the two periods                 
plot(D.opc, PN1, xlab = expression(paste("dD (", mu, "m)")), ylab = expression(paste("dN/dln(D) (cm"^"-3",")")), log ="xy", yaxt="n", type ="o", col ="red")
par(new=TRUE)
plot(D.opc, PN2, type="o", yaxt="n",xaxt="n", ann =FALSE, col ="blue", log="xy")

y1 <- floor(log10(range(PN1)))
pow_y <- seq(y1[1], y1[2]+1)
ticksat_y <- as.vector(sapply(pow_y, function(p) (1:10)*10^p))
axis(2, 10^pow_y)
axis(2, ticksat_y, labels=NA, tcl=-0.25, lwd=0, lwd.ticks=1)

x1 <- floor(log10(range(D.opc)))
pow_x <- seq(x1[1], x1[2]+1)
ticksat_x <- as.vector(sapply(pow_x, function(p) (1:10)*10^p))
axis(1, 10^pow_x)
axis(1, ticksat_x, labels=NA, tcl=-0.25, lwd=0, lwd.ticks=1)

legend(3.5, 800, legend =c("first subset", "second subset"), col =c("red", "blue"), lty=c(1,1), cex=0.8)






#plot average particle surface distribution for the two periods                 
plot(D.opc, PS1, xlab = expression(paste("dD (", mu, "m)")), ylab = expression(paste("dS/dln(D) ( ", mu, "m"^"2","cm"^"-3",")")), log ="xy",  yaxt="n", type ="o", col ="red")
par(new=TRUE)
plot(D.opc, PS2, type="o", ann =FALSE, yaxt="n",xaxt="n",  col ="blue", log="xy")

y1 <- floor(log10(range(PS1)))
pow_y <- seq(y1[1], y1[2]+1)
ticksat_y <- as.vector(sapply(pow_y, function(p) (1:10)*10^p))
axis(2, 10^pow_y)
axis(2, ticksat_y, labels=NA, tcl=-0.25, lwd=0, lwd.ticks=1)

x1 <- floor(log10(range(D.opc)))
pow_x <- seq(x1[1], x1[2]+1)
ticksat_x <- as.vector(sapply(pow_x, function(p) (1:10)*10^p))
axis(1, 10^pow_x)
axis(1, ticksat_x, labels=NA, tcl=-0.25, lwd=0, lwd.ticks=1)

legend(3.5, 800, legend =c("first subset", "second subset"), col =c("red", "blue"), lty=c(1,1), cex=0.8)






#plot average particle volume distribution for the two periods                 
plot(D.opc, PV1, xlab = expression(paste("dD (", mu, "m)")), ylab = expression(paste("dV/dln(D) ( ", mu, "m"^"3","cm"^"-3",")")), log ="xy",  yaxt="n",  type ="o", col ="red")
par(new=TRUE)
plot(D.opc, PV2, type="o", yaxt="n",xaxt="n", ann =FALSE, col ="blue", log="xy")

y1 <- floor(log10(range(PV1)))
pow_y <- seq(y1[1], y1[2]+1)
ticksat_y <- as.vector(sapply(pow_y, function(p) (1:10)*10^p))
axis(2, 10^pow_y)
axis(2, ticksat_y, labels=NA, tcl=-0.25, lwd=0, lwd.ticks=1)

x1 <- floor(log10(range(D.opc)))
pow_x <- seq(x1[1], x1[2]+1)
ticksat_x <- as.vector(sapply(pow_x, function(p) (1:10)*10^p))
axis(1, 10^pow_x)
axis(1, ticksat_x, labels=NA, tcl=-0.25, lwd=0, lwd.ticks=1)

legend(3.5, 800, legend =c("first subset", "second subset"), col =c("red", "blue"), lty=c(1,1), cex=0.8)







#find moments of distribution
library(e1071)
MP01 <-moment(as.numeric(PN1[1,]), order = 0, center = TRUE, absolute = FALSE, na.rm =TRUE)
MP02 <-moment(as.numeric(PN2[1,]), order = 0, center = FALSE, absolute = FALSE, na.rm =TRUE)
MP11 <-moment(as.numeric(PN1[1,]), order = 1, center = FALSE, absolute = FALSE, na.rm =TRUE)
MP12 <-moment(as.numeric(PN2[1,]), order = 1, center = FALSE, absolute = FALSE, na.rm =TRUE)
MP21 <-moment(as.numeric(PN1[1,]), order = 2, center = FALSE, absolute = FALSE, na.rm =TRUE)
MP22 <-moment(as.numeric(PN2[1,]), order = 2, center = FALSE, absolute = FALSE, na.rm =TRUE)
MP31 <-moment(as.numeric(PN1[1,]), order = 3, center = FALSE, absolute = FALSE, na.rm =TRUE)
MP32 <-moment(as.numeric(PN2[1,]), order = 3, center = FALSE, absolute = FALSE, na.rm =TRUE)

