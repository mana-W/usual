#cell 19 2d plot

find_ctrl <- function(x){
  n=1
  i=0
  while(i==0){
    if(length(which(geneslist[[n]]==x))==0){
      n=n+1
    }else{
      ctrlset = sample(geneslist[[n]][geneslist[[n]]!=x],100)
      i=1
    }
  }
  return(c(ctrlset))
}

cellscore_fun <- function(testmoud,ctrlmoudcells){
   mean_Gt <- mean(Erexp[testmoud,cells])
   mean_Gc <- mean(Erexp[ctrlmoud,cells])
   SC <- mean_Gt - mean_Gc
   return(SC)
}
cellscore_fun_all <- function(x){
  #scmes2 <- cellscore_fun(moudle$MES2[moudle$MES2!="None"],MES2ctrl,x)
  scmes <- cellscore_fun(moudle$MES[moudle$MES!="None"],MESctrl,x)
  scac <- cellscore_fun(moudle$AC[moudle$AC!="None"],ACctrl,x)
  scopc <- cellscore_fun(moudle$OPC[moudle$OPC!="None"],OPCctrl,x)
  scnpc <- cellscore_fun(moudle$NPC[moudle$NPC!="None"],NPCctrl,x)
  #scnpc2 <- cellscore_fun(moudle$NPC2[moudle$NPC2!="None"],NPC2ctrl,x)
  scg1s <- cellscore_fun(moudle$G1.S[moudle$G1.S!="None"],G1Sctrl,x)
  scg2m <- cellscore_fun(moudle$G2.M[moudle$G2.M!="None"],G2Mctrl,x)
  return(data.frame(
            MES=scmes,
            #MES1=scmes1,
            AC=scac,
            OPC=scopc,
            NPC=scnpc,
            #NPC2=scnpc2,
            G1.S=scg1s,
            G2.M=scg2m
        ))
}
x_pos_fun <- function(x){
  y <- as.numeric(unlist(x[c(1:4)]))
  if(as.character(x[8])=="MES"){
    return(log2((y[1]-y[2])+1))
  }else if(as.character(x[8])=="OPC"){
    return(log2((y[3]-y[4])+1))
  }else if(as.character(x[8])=="AC"){
    return(-log2((y[2]-y[1])+1))
  }else if(as.character(x[8])=="NPC"){
    return(-log2((y[4]-y[3])+1))
  }
}


Eexp <- read.table("GSM3828672_Smartseq2_GBM_IDHwt_processed_TPM.tsv",row.names = 1,header=T)
moudle <- read.csv("moudle_4.csv",header=T)
exp <- t(apply(Eexp,1,function(x){(2^x-1)*10}))
#Eexp <- t(apply(exp,1,function(x){log2((x/10)+1)}))
Erexp <- t(apply(Eexp,1,function(x){x-mean(x)}))

allgenes <- apply(exp,1,function(x){
      x = as.numeric(unlist(x))
      x=x[x!=0]
      return(log2(mean(x)+1))
      })
#allgenes <- allgenes[allgenes>0]
allgenes <- sort(allgenes,decreasing = T)
allgenes <- allgenes[allgenes > 4]
geneslist <- list()
n <- round(length(allgenes)/30)
i=1
index=1
while(i<31){
  if(i<30){
    genes <- names(allgenes)[c(index:(index+n))]
    geneslist <- c(geneslist,list(genes))
    index <- index+n+1
    i=i+1
  }else{
    genes <- names(allgenes)[c(index:length(allgenes))]
    geneslist <- c(geneslist,list(genes))
    i=i+1
  }
}

MESctrl <- unlist(lapply(moudle[["MES"]][moudle[["MES"]]!="None"],find_ctrl))
#MES1ctrl <- unlist(lapply(moudle[["MES1"]][moudle[["MES1"]]!="None"],find_ctrl))
ACctrl <- unlist(lapply(moudle[["AC"]][moudle[["AC"]]!="None"],find_ctrl))
OPCctrl <- unlist(lapply(moudle[["OPC"]][moudle[["OPC"]]!="None"],find_ctrl))
NPCctrl <- unlist(lapply(moudle[["NPC"]][moudle[["NPC"]]!="None"],find_ctrl))
#NPC2ctrl <- unlist(lapply(moudle[["NPC2"]][moudle[["NPC2"]]!="None"],find_ctrl))
G1Sctrl <- unlist(lapply(moudle[["G1.S"]][moudle[["G1.S"]]!="None"],find_ctrl))
G2Mctrl <- unlist(lapply(moudle[["G2.M"]][moudle[["G2.M"]]!="None"],find_ctrl))

SCall <- do.call("rbind",lapply(colnames(Erexp),cellscore_fun_all))
row.names(SCall)<-colnames(Erexp)
write.csv(SCall,"SCall.csv",row.names = F,quote = F)

SCall$D_y <- apply(SCall,1,function(x){max(x[3],x[4])-max(x[1],x[2])})
SCall$type <- apply(SCall,1,function(x){
                      i=which(x[c(1:4)]==max(x[c(1:4)]))
                      return(names(SCall)[i])
                    })

SCall$x <- apply(SCall,1,x_pos_fun)
ggplot(SCall,aes(x,D_y))+geom_point()+theme_bw()+xlab("AC------MES")+ylab("AC------OPC")
write.csv(SCall,"SCall.csv",row.names = F,quote = F)


