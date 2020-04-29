library(dplyr)
library(readxl)
require(reshape)
library(RCurl)
here::here() %>% setwd()

x <- getURL('https://raw.githubusercontent.com/MinCiencia/Datos-COVID19/master/output/producto3/CasosTotalesCumulativo.csv')
infectados <- read.csv(text = x)

y <- getURL('https://raw.githubusercontent.com/MinCiencia/Datos-COVID19/master/output/producto14/FallecidosCumulativo.csv')
muertos <- read.csv(text=y)

muertos = cbind(muertos[,1],data.frame(matrix(0,nrow=dim(muertos)[1],ncol=19)),muertos[-1])

datos = muertos
num = datos[,-1]
nm= num[,1]
i=2

while (i<dim(num)[2]+1){
    nm=cbind(nm,num[,i]-num[,i-1])
    i=i+1
}

casi=cbind(muertos[1],nm)
i=1
dat=c()

while(i<17+1){
    aux=cbind(matrix(rep(casi[i,1]),length(casi[1,][-1])),t(casi[i,][-1]))
    dat=rbind(dat,aux)
    i=i+1
}

datos=data.frame(dat)
b = rename(datos, c(X1='new_deaths',V1='country'))

datos = infectados
num = datos[,-1]
nm= num[,1]
i=2

while (i<dim(num)[2]+1){
    nm=cbind(nm,num[,i]-num[,i-1])
    i=i+1
}

casi=cbind(muertos[1],nm)
i=1
dat=c()
while(i<17+1){
    aux=cbind(t(casi[i,][-1]))
    dat=rbind(dat,aux)
    i=i+1
}
inf =dat
c = data.frame(inf)
d = rename(c, c(X1='new_cases'))
df = cbind(b,d)

write.csv(df,file = "res.csv",row.names = F)

muTransform <- function(zMedian){
  mu <- log(zMedian) 
}
sigmaTransform <- function(zMean, mu){
  sigma <- sqrt(2*(log(zMean) - mu))
}
# Hospitalisation to death distribution
hospitalisation_to_death_truncated <- function(x, mu, sigma) {
  plnorm(x + 1, mu, sigma) - plnorm(x, mu, sigma)
}
hospitalisation_to_death_truncated_low <- function(x){
  hospitalisation_to_death_truncated(x, muLow, sigmaLow)
}
hospitalisation_to_death_truncated_mid <- function(x){
  hospitalisation_to_death_truncated(x, muMid, sigmaMid)
}
hospitalisation_to_death_truncated_high <- function(x){
  hospitalisation_to_death_truncated(x, muHigh, sigmaHigh)
}
# Función del CFR
scale_cfr <- function(data_1_in, delay_fun){
  case_incidence <- data_1_in$new_cases
  death_incidence <- data_1_in$new_deaths
  cumulative_known_t <- 0 # cumulative cases with known outcome at time tt
  # Sum over cases up to time tt
  for(ii in 1:nrow(data_1_in)){
    known_i <- 0 # number of cases with known outcome at time ii
    for(jj in 0:(ii - 1)){
      known_jj <- (case_incidence[ii - jj]*delay_fun(jj))
      known_i <- known_i + known_jj
    }
    cumulative_known_t <- cumulative_known_t + known_i # Tally cumulative known
  }
  # naive CFR value
  b_tt <- sum(death_incidence)/sum(case_incidence) 
  # corrected CFR estimator
  p_tt <- sum(death_incidence)/cumulative_known_t
  data.frame(nCFR = b_tt, cCFR = p_tt, total_deaths = sum(death_incidence), 
             cum_known_t = round(cumulative_known_t), total_cases = sum(case_incidence))
    }
underReportingEstimates <- function(data, delay_fun){ 
  dplyr::group_by(data, country) %>%
    dplyr::do(scale_cfr(., delay_fun)) %>%
    dplyr::filter(cum_known_t > 0 & cum_known_t >= total_deaths)  %>%
    dplyr::mutate(nCFR_UQ = binom.test(total_deaths, total_cases)$conf.int[2],
                  nCFR_LQ = binom.test(total_deaths, total_cases)$conf.int[1],
                  cCFR_UQ = binom.test(total_deaths, cum_known_t)$conf.int[2],
                  cCFR_LQ = binom.test(total_deaths, cum_known_t)$conf.int[1],
                  underreporting_estimate = cCFRBaseline / (100*cCFR),
                  lower = cCFREstimateRange[1] / (100 * cCFR_UQ),
                  upper = cCFREstimateRange[2] / (100 * cCFR_LQ),
                  quantile25 = binom.test(total_deaths, cum_known_t, conf.level = 0.5)$conf.int[1],
                  quantile75 = binom.test(total_deaths, cum_known_t, conf.level = 0.5)$conf.int[2]) %>% 
    dplyr::filter(total_deaths > 0)}

cCFRBaseline <- 1.4
cCFREstimateRange <- c(1.2, 1.7)
zmeanLow <- 8.7
zmedianLow <- 6.7
muLow <- muTransform(zmedianLow)
sigmaLow <- sigmaTransform(zmeanLow, muLow)
zmeanMid <- 13
zmedianMid <- 9.1
muMid <- muTransform(zmedianMid)
sigmaMid <- sigmaTransform(zmeanMid, muMid)
zmeanHigh <- 20.9
zmedianHigh <- 13.7
muHigh <- muTransform(zmedianHigh)
sigmaHigh <- sigmaTransform(zmeanHigh, muHigh)

allTogetherClean=read.csv('res.csv')

allTogetherLow <- underReportingEstimates(allTogetherClean, hospitalisation_to_death_truncated_low) 
allTogetherMid <- underReportingEstimates(allTogetherClean, hospitalisation_to_death_truncated_mid) 
allTogetherHigh <- underReportingEstimates(allTogetherClean, hospitalisation_to_death_truncated_high)

finalRes <- dplyr::tibble(
      country = allTogetherMid$country,
      total_cases = allTogetherMid$total_cases,
      total_deaths = allTogetherMid$total_deaths,
      underreporting_estimate  = pmin(allTogetherLow$underreporting_estimate, allTogetherMid$underreporting_estimate, allTogetherHigh$underreporting_estimate),
      lower = pmin(allTogetherLow$lower, allTogetherMid$lower, allTogetherHigh$lower),
      upper = pmax(allTogetherLow$upper, allTogetherMid$upper, allTogetherHigh$upper))

names(finalRes)<-c('Localidad','Casos_totales','Muertes_Totales','Tasa_de_no_reportados','Cota_inferior_IC','Cota_superior_IC')

nfecha= length(subset(allTogetherClean,country=='Arica y Parinacota')[,1])

fecha=data.frame(seq(as.Date('2020-03-03'), by='day', length=nfecha))

x <- getURL('https://raw.githubusercontent.com/MinCiencia/Datos-COVID19/master/output/producto3/CasosTotalesCumulativo.csv')
infectados <- read.csv(text = x)

y <- getURL('https://raw.githubusercontent.com/MinCiencia/Datos-COVID19/master/output/producto14/FallecidosCumulativo.csv')
muertos <- read.csv(text=y)

muertos = cbind(muertos[,1],data.frame(matrix(0,nrow=dim(muertos)[1],ncol=19)),muertos[-1])

names(muertos)<-c('s')

k=1
while(k<17+1){
nombre = muertos[k,1] #desde 1 hasta 17
i=1
datos_por_fecha = c()
while (i<nfecha+1){
    df2=subset(allTogetherClean,country == nombre)[1:i,]
    write.csv(df2,file = "gato.csv",row.names = F)
    base=read.csv('gato.csv')
    d1<-try(underReportingEstimates(base, hospitalisation_to_death_truncated_low),silent=TRUE)
    d2<-try(underReportingEstimates(base, hospitalisation_to_death_truncated_mid),silent=TRUE)
    d3<-try(underReportingEstimates(base, hospitalisation_to_death_truncated_high),silent=TRUE)  
    if (length(d1) == 1 | length(d2) == 1 | length(d3) == 1){
        res <- data.frame(NA,NA,NA,NA,NA,NA,NA,NA,NA)
        names(res)<-c( 'nCFR' ,'nCFR_LQ' ,'nCFR_UQ' ,'cCFR' ,'cCFR_LQ' ,'cCFR_UQ','Tasa_de_no_reportados','Cota_inferior_IC','Cota_superior_IC')
        } else {
            allTogetherLow <- underReportingEstimates(base, hospitalisation_to_death_truncated_low)
            allTogetherMid <- underReportingEstimates(base, hospitalisation_to_death_truncated_mid) 
            allTogetherHigh <- underReportingEstimates(base, hospitalisation_to_death_truncated_high)
            finalRes <- dplyr::tibble(
              country = allTogetherMid$country,
              total_cases = allTogetherMid$total_cases,
              total_deaths = allTogetherMid$total_deaths,
              underreporting_estimate  = pmin(allTogetherLow$underreporting_estimate, allTogetherMid$underreporting_estimate, allTogetherHigh$underreporting_estimate),
              lower = pmin(allTogetherLow$lower, allTogetherMid$lower, allTogetherHigh$lower),
              upper = pmax(allTogetherLow$upper, allTogetherMid$upper, allTogetherHigh$upper))
              names(finalRes)<-c('Localidad','Casos_totales','Muertes_Totales','Tasa_de_no_reportados','Cota_inferior_IC','Cota_superior_IC')
              extra = cbind(allTogetherHigh[1,2],allTogetherHigh[1,8],allTogetherHigh[1,7],allTogetherHigh[1,3],allTogetherHigh[1,10],allTogetherHigh[1,9])
              if(dim(finalRes)[1]==0){
                  res <- data.frame(NA,NA,NA,NA,NA,NA,NA,NA,NA)
                  names(res)<-c( 'nCFR' ,'nCFR_LQ' ,'nCFR_UQ' ,'cCFR' ,'cCFR_LQ' ,'cCFR_UQ','Tasa_de_no_reportados','Cota_inferior_IC','Cota_superior_IC')
              } else {
                  resu <- cbind(finalRes[1,4:6])
                  if (resu[1,3]>1){
                        resu[1,3]=1
                    }
                  res<-cbind(extra,resu)
              }
        }
        datos_por_fecha=rbind(datos_por_fecha,res)
    i=i+1
    }
q=cbind(fecha,data.frame(cbind(t(infectados[k,][-1]),t(muertos[k,][-1]))))
names(q)<-c('Fecha','Casos acumulados a la fecha','Muertes acumuladas a la fecha')
write.csv(cbind(q,datos_por_fecha),file = paste(nombre,'.csv',sep=''),row.names = F)
k=k+1
}
