library(IMvigor210CoreBiologies)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(ggsci)
library(agricolae)

data(cds)

expMatrix <- counts(cds)
eff_length2 <- fData(cds)[,c("entrez_id","length","symbol")]
rownames(eff_length2) <- eff_length2$entrez_id
head(eff_length2)
feature_ids <- rownames(expMatrix)
expMatrix <- expMatrix[feature_ids %in% rownames(eff_length2),]
mm <- match(rownames(expMatrix),rownames(eff_length2))
eff_length2 <- eff_length2[mm,]

x <- expMatrix/eff_length2$length
matrix <- t(t(x)/colSums(x))*1e6
summary(duplicated(rownames(matrix)))

matrix <- IOBR::anno_eset(eset = matrix,
                        annotation = eff_length2,
                        symbol = "symbol",
                        probe = "entrez_id",
                        method = "mean")
tumor_type <- "blca"
if(max(matrix)>100) matrix <- log2(matrix+1)

meta <- pData(cds)


meta2 <- meta[grep("NE",meta$Overall_Response,invert = T),]

column <- table(meta2$Overall_Response,meta2$Lund) %>% data.frame()
colnames(column) <- c("response","Lund","Freq")

p1 <- ggplot(column,aes(Lund,weight=Freq,fill=response))+geom_bar(position="fill")+#coord_flip()+
  theme(panel.grid.major =element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black")) + #scale_fill_npg()
  scale_fill_manual(values = brewer.pal(5,"Set3")[-2]) #[c(1,5,17,14,21)]
p1 + theme_bw() + labs(y="Proportion") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))



all(colnames(matrix)==rownames(meta))
meta$CD39 <- t(matrix)[,"ENTPD1"]
meta$PDL1 <- t(matrix)[,"CD274"]

ggplot(meta,aes(x=Lund,y=CD39,fill=Lund))+geom_boxplot()+
  geom_jitter(size=0.5)+scale_fill_npg()+theme_bw()
fit <- aov(CD39~Lund,data=meta)
summary(fit)
out <- LSD.test(result,"Lund")
plot(out)


ggplot(meta,aes(x=Lund,y=PDL1,fill=Lund))+geom_boxplot()+
  geom_jitter(size=0.5)+scale_fill_npg()+theme_bw()
fit <- aov(PDL1~Lund,data=meta)
summary(fit)
out <- LSD.test(result,"Lund")
plot(out)