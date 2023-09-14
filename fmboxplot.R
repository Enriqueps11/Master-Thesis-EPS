###############################################################
#This function creates the multivariate boxplot taking into account the indexes obtained in indexes.R
source("C:/Users/enriq/Desktop/UNIVERSIDAD/Máster en Big Data/TFM/Code/indexes.R")

param <- function(X, t, rangeval, nbasis, norder) {
  results<-indM(X, t, rangeval, nbasis, norder)
  indexes<-results$ind
  sdata<-results$sdata
  n<-nrow(indexes)
  hypo<-sort(indexes$dtaMHI)
  epi<-sort(indexes$dtaMEI)
  #75% con mayor Hypograph
  q_h<-hypo[round(n/4)+1]
  sorted_h<-list()
  for (i in 1:dim(X)[3]) {
    sorted_h[[i]] <- as.data.frame(subset(sdata[,,i],indexes$dtaMHI>=q_h))
  }
  #75% con menor Epigraph
  q_e<-epi[round(n*3/4)+1]
  sorted_e<-list()
  for (i in 1:dim(X)[3]) {
    sorted_e[[i]] <- as.data.frame(subset(sdata[,,i],indexes$dtaMEI<=q_e))
  }
  #Bucle para crear q1 y q3
  q1<-array(dim = c(dim(X)[3],length(t)))
  for (j in 1:dim(X)[3]) {
    data<-as.data.frame(sorted_h[j])
    for (i in 1:ncol(data)) {
      q1[j,i]<-min(data[,i])
    }
  }
  q3<-array(dim = c(dim(X)[3],length(t)))
  for (j in 1:dim(X)[3]) {
    data<-as.data.frame(sorted_e[j])
    for (i in 1:ncol(data)) {
      q3[j,i]<-max(data[,i])
    }
  }
  
  #IQR: q3-q1
  IQR<-array(dim = c(dim(X)[3],length(t)))
  IQR <- q3-q1
  
  
  #Lower whisker 
  w1 <- q1 - 0.8*IQR
  #Upper whisker 
  w3 <- q3 + 0.8*IQR
  
  #Encontrar todos los outliers y concatenarlos
  outliers<-c()
  for (i in 1:dim(X)[3]) {
    d<-sdata[,,i]
    w1_<-w1[i,]
    w3_<-w3[i,]
    for (k in 1:nrow(d)) {
      if((sum(d[k,]>w3_|d[k,]<w1_)>0)&(!is.element(k, outliers))){outliers<-c(outliers,k)}
    }
  }
  outliers<-sort(outliers)
  
  res<-list("q1"=q1,"q3"=q3,"w1"=w1,"w3"=w3,"sdata"=sdata,"outliers"=outliers)
  return(res)
}
fmboxplot <- function(i,q1,q3,w1,w3,sdata,outliers) {
  d<-sdata[,,i]
  w1<-w1[i,]
  w3<-w3[i,]
  p <- ggplot()
  p <- p + geom_line(aes_string(x=t,y = q1[i,]), color="#2dd9ac")
  p <- p + geom_line(aes_string(x=t,y = q3[i,]), color="#213c97")
  p <- p + geom_ribbon(aes(x= t, ymin = q1[i,], ymax = q3[i,]), fill = "lightskyblue1", alpha = 0.5)
  p <- p + geom_line(aes_string(x=t,y = w1), color="#59bec3")
  p <- p + geom_line(aes_string(x=t,y = w3), color="#060552")
  if(!is.null(outliers)){
    for (k in outliers) {
      p <- p + geom_line(aes_string(x=t,y = d[k,]), color="red")
    }
  }
  for (j in 1:10) {
    if(!is.element(j,outliers)){p <- p + geom_line(aes_string(x=t,y = d[j,]), color="#FF00FF",linetype = "solid")}
  }
  return(p)
}

fmboxplot_0 <- function(i,q1,q3,w1,w3,sdata,outliers) {
  d<-sdata[,,i]
  w1<-w1[i,]
  w3<-w3[i,]
  p <- ggplot()
  p <- p + geom_line(aes_string(x=t,y = q1[i,]), color="#2dd9ac")
  p <- p + geom_line(aes_string(x=t,y = q3[i,]), color="#213c97")
  p <- p + geom_ribbon(aes(x= t, ymin = q1[i,], ymax = q3[i,]), fill = "lightskyblue1", alpha = 0.5)
  p <- p + geom_line(aes_string(x=t,y = w1), color="#59bec3")
  p <- p + geom_line(aes_string(x=t,y = w3), color="#060552")
  if(!is.null(outliers)){
    for (k in outliers) {
      p <- p + geom_line(aes_string(x=t,y = d[k,]), color="red")
    }
  }
  return(p)
}

###############################################################################
#Función para datos individuales
###############################################################################

fmboxplot_i <- function(X, t, rangeval, nbasis, norder) {
  # we obtain the indexes for the observations
  results<-ind(X, t, rangeval, nbasis, norder)
  indexes<-results$ind
  sdata<-results$sdata
  n<-nrow(indexes)
  hypo<-sort(indexes$dtaMHI)
  epi<-sort(indexes$dtaMEI)
  #75% con mayor Hypograph
  q_h<-hypo[round(n/4)+1]
  sorted_h <- as.data.frame(subset(sdata,indexes$dtaMHI>=q_h))
  #75% con menor Epigraph
  q_e<-epi[round(n*3/4)+1]
  sorted_e <- as.data.frame(subset(sdata,indexes$dtaMEI<=q_e))
  #Bucle para crear q1 y q3
  q1<-rep(0,length(t))
  for (i in 1:ncol(sorted_h)) {
    q1[i]<-min(sorted_h[,i])
  }
  q3<-rep(0,length(t))
  for (i in 1:ncol(sorted_e)) {
    q3[i]<-max(sorted_e[,i])
  }
  
  #IQR: q3-q1
  IQR <- q3-q1
  #Lower whisker 
  w1 <- q1 - 0.8*IQR
  #Upper whisker 
  w3 <- q3 + 0.8*IQR
  
  #Dataframe con todos los parámetros
  values<-data.frame(q1,q3,w1,w3)
  
  #Calculo outliers
  outliers<-c()
  for (k in 1:nrow(sdata)) {
    if((sum(sdata[k,]>w3|sdata[k,]<w1)>0)&(!is.element(k, outliers))){outliers<-c(outliers,k)}
  }
  
  #Plot
  p <- ggplot()
  p <- p + geom_line(aes_string(x=t,y = values$q1), color="#2dd9ac")
  p <- p + geom_line(aes_string(x=t,y = values$q3), color="#213c97")
  p <- p + geom_ribbon(aes(x= t, ymin = values$q1, ymax = values$q3), fill = "lightskyblue1", alpha = 0.5)
  p <- p + geom_line(aes_string(x=t,y = values$w1), color="#59bec3")
  p <- p + geom_line(aes_string(x=t,y = values$w3), color="#060552")
  if(!is.null(outliers)){
    for (i in outliers) {
    p <- p + geom_line(aes_string(x=t,y = sdata[i,]), color="red")
    }
  }
  for (j in 1:10) {
    if(!is.element(j,outliers)){p <- p + geom_line(aes_string(x=t,y = sdata[j,]), color="#FF00FF",linetype = "solid")}
  }
  res<-list("plot"=p,"outliers"=outliers)
  return(res)
}

fmboxplot_i0 <- function(X, t, rangeval, nbasis, norder) {
  # we obtain the indexes for the observations
  results<-ind(X, t, rangeval, nbasis, norder)
  indexes<-results$ind
  sdata<-results$sdata
  n<-nrow(indexes)
  hypo<-sort(indexes$dtaMHI)
  epi<-sort(indexes$dtaMEI)
  #75% con mayor Hypograph
  q_h<-hypo[round(n/4)+1]
  sorted_h <- as.data.frame(subset(sdata,indexes$dtaMHI>=q_h))
  #75% con menor Epigraph
  q_e<-epi[round(n*3/4)+1]
  sorted_e <- as.data.frame(subset(sdata,indexes$dtaMEI<=q_e))
  #Bucle para crear q1 y q3
  q1<-rep(0,length(t))
  for (i in 1:ncol(sorted_h)) {
    q1[i]<-min(sorted_h[,i])
  }
  q3<-rep(0,length(t))
  for (i in 1:ncol(sorted_e)) {
    q3[i]<-max(sorted_e[,i])
  }
  
  #IQR: q3-q1
  IQR <- q3-q1
  #Lower whisker 
  w1 <- q1 - 0.8*IQR
  #Upper whisker 
  w3 <- q3 + 0.8*IQR
  
  #Dataframe con todos los parámetros
  values<-data.frame(q1,q3,w1,w3)
  
  #Calculo outliers
  outliers<-c()
  for (k in 1:nrow(sdata)) {
    if((sum(sdata[k,]>w3|sdata[k,]<w1)>0)&(!is.element(k, outliers))){outliers<-c(outliers,k)}
  }
  
  #Plot
  p <- ggplot()
  p <- p + geom_line(aes_string(x=t,y = values$q1), color="#2dd9ac")
  p <- p + geom_line(aes_string(x=t,y = values$q3), color="#213c97")
  p <- p + geom_ribbon(aes(x= t, ymin = values$q1, ymax = values$q3), fill = "lightskyblue1", alpha = 0.5)
  p <- p + geom_line(aes_string(x=t,y = values$w1), color="#59bec3")
  p <- p + geom_line(aes_string(x=t,y = values$w3), color="#060552")
  if(!is.null(outliers)){
    for (i in outliers) {
      p <- p + geom_line(aes_string(x=t,y = sdata[i,]), color="red")
    }
  }
  res<-list("plot"=p,"outliers"=outliers)
  return(res)
}

###############################################################################
#Función para TPR y FPR
###############################################################################

rates <- function(outliers,outliers_i) {
  real_outliers<-c(1:10)
  if(is.null(outliers)){
    TPR_multi<-0
    FPR_multi<-0
  }else{
    correct<-0
    for (i in outliers) {
      if(i %in% real_outliers){correct<-correct+1}
    }
    TPR_multi<-correct/length(real_outliers) # Porcentaje de acierto
    FPR_multi<-(length(outliers)-correct)/length(outliers) #Porcentaje false positives
  }
  if(is.null(outliers_i)){
    TPR_ind<-0
    FPR_ind<-0
  }else{
    correct<-0
    for (i in outliers_i) {
      if(i %in% real_outliers){correct<-correct+1}
    }
    TPR_ind<-correct/length(real_outliers) # Porcentaje de acierto
    FPR_ind<-(length(outliers_i)-correct)/length(outliers_i) #Porcentaje false positives
  }
  TPR<-c(TPR_multi,TPR_ind)
  FPR<-c(FPR_multi,FPR_ind)
  df <- data.frame("True positive rate" = TPR, "False positive rate" = FPR)
  nombres_columnas <- colnames(df)
  nombres_columnas_ajustados <- gsub("\\.", " ", nombres_columnas)
  colnames(df) <- nombres_columnas_ajustados
  nombres_filas <- c("Multivariate","Individual")
  rownames(df) <- nombres_filas
  return(df)
}

rates0 <- function(outliers,outliers_i) {
  TPR_multi<-0
  FPR_multi<-length(outliers)/100
  TPR_ind<-0
  FPR_ind<-length(outliers_i)/100
  TPR<-c(TPR_multi,TPR_ind)
  FPR<-c(FPR_multi,FPR_ind)
  df <- data.frame("True positive rate" = TPR, "False positive rate" = FPR)
  nombres_columnas <- colnames(df)
  nombres_columnas_ajustados <- gsub("\\.", " ", nombres_columnas)
  colnames(df) <- nombres_columnas_ajustados
  nombres_filas <- c("Multivariate","Individual")
  rownames(df) <- nombres_filas
  return(df)
}