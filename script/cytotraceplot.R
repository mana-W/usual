

#CytoTRACE result 
library(cowplot)
library(ggplot2)

cytoTRACEplotgene <- function(results,phenot,gene,...){
df <- data.frame(results$CytoTRACE)
df["Celltype"] <- phenot
df["Gene"] <- results$exprMatrix[gene,]
if(is.null(cols)){
fig1 <- ggplot(df,aes(x=results.CytoTRACE,y=Gene,color=Celltype))+
  geom_point()+geom_smooth(method= loess,se = F)+
  ggtitle(gene)+theme_bw()+scale_colour_brewer(palette = "Set1")+
  theme(legend.position = c(0.82, 0.76),legend.background = element_rect(color = "black", linetype = "solid", size = .35))
}else{
  fig1 <- ggplot(df,aes(x=results.CytoTRACE,y=Gene,color=Celltype))+
    geom_point()+geom_smooth(method= loess,se = F)+
    ggtitle(gene)+theme_bw()+scale_colour_manual(values = cols)+
    theme(legend.position = c(0.82, 0.76),legend.background = element_rect(color = "black", linetype = "solid", size = .35))
  
}
fig2 <- ggplot(df,aes(x=results.CytoTRACE,y=Gene))+
  geom_point()+geom_smooth(method= loess,se = F,aes(color="red4"))+
  ggtitle(gene)+theme_bw()+NoLegend()
fig <- plot_grid(fig1,fig2)
return(fig)
}
test <- cytoTRACEplotgene(results,phenot,"Ocm")
pdf("Ocm_curve.pdf",width = 7.5,height = 3.9)
test
dev.off()
