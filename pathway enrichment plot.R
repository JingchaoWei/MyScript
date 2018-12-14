library(stringr)
library(magrittr)
library(ggplot2)
#柱状图
tmp <- kk@result[kk@result$qvalue<0.01,]#提取数据框
#将generatio由字符转化为数值
str(tmp)
Char2Num <- function(i) {
  tmp1 <- tmp[i,3]%>%str_split(pattern = '/')%>%unlist()%>%as.numeric()
  tmp2 <- divide_by(tmp1[1],tmp1[2])
  return(tmp2)
}
tmp3 <- vector()
for (i in 1:nrow(tmp)) {
  tmp4 <- Char2Num(i)
  tmp3 <- c(tmp3,tmp4)
}
tmp3
tmp$GeneRatio <- tmp3
str(tmp)
#plot
library(ggplot2)
ggplot(tmp,aes(x=GeneRatio,y=Description))+
  geom_point(aes(size=Count,color=qvalue))+
  scale_colour_gradient(low="blue",high="red")+ 
  scale_x_continuous()+
  labs( color=expression(qvalue), x="Gene Ratio",
        y="Pathway name", title="Pathway enrichment")+ 
  theme_bw()+
  theme(axis.text.y = element_text(size = rel(1)), 
        axis.title.x = element_text(size=rel(1)), 
        axis.title.y = element_blank() ) 
#气泡图
rm(list=ls()) 
library(Cairo) 
library(stringr) 
pathway = read.table("./enh_statistics/A549_KEGG.tsv",header=T,sep="\t") 
pathway$Term<-str_split_fixed(pathway$Term,":",2)[,2] 
png_path="./figure/KEGG.png" 
CairoPNG(png_path, width = 5.9, height = 3, units='in', dpi=600) 
ggplot(pathway,aes(x=Fold.Enrichment,y=Term))+
  geom_point(aes(size=Count,color=-1*log10(PValue)))+ 
  scale_colour_gradient(low="green",high="red")+
  labs(color=expression(-log[10](P.value)),
       x="Fold enrichment",y="Pathway name",title="Pathway enrichment")+
  theme_bw()+ theme( axis.text.y = element_text(size = rel(1.3)),
                     axis.title.x = element_text(size=rel(1.3)), 
                     axis.title.y = element_blank()) 
                                                      