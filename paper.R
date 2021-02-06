


#---fig 1 #------------ 
library(ggplot2)
library(EasyStat)
library(tidyverse)
data = read.csv("./env/env5_NEW_change_nine.csv",row.names = 1)
head(data)

axis_order = c("Control","A","O","S","AO","AS","OS","AOS")


data$group = gsub("CK","Control",data$group)
data$group = gsub("Y","Raw soil",data$group)



data = data %>% filter(group %in% axis_order)
data$group = as.factor(data$group)
data = data[,c(1:2,5,7,10)]
head(data)
colnames(data)[3:5] = c("AFe(ug/g)","ASi(ug/g)","DOC(mg/kg)")


# result = MuiKwWlx(data = data,num = c(3:5))
# result
result = MuiaovMcomper(data = data,num = c(3:5))
result

result <- ord_sig(data = result,ID = NULL)
result$AFe = gsub(" ","",result$AFe)
result$ASi = gsub(" ","",result$ASi)
result$DOC = gsub(" ","",result$DOC)


result1 = FacetMuiPlotresultBox(data = data,num = c(3:5),result = result,sig_show ="abc",ncol = 1 )
p = result1[[1]] + theme_bw()+ scale_x_discrete(limits = axis_order)
p  


A = result1[[2]]


aa = c()
i = 1
A$group

for (i in 1:length(A$group)) {
  if (A$group[i] == "Control" ) {
    aa[i] = "Control"
  }
  if (A$group[i] %in% c("A","B","O","S") ) {
    aa[i] = "single"
  }
  if (A$group[i] %in% c("AB","AO","AS","BO","BS","OS") ) {
    aa[i] = "two"
  }
  if (A$group[i] %in% c("ABO","ABS","AOS","BOS") ) {
    aa[i] = "three"
  }
  if (A$group[i] %in% c("ABOS") ) {
    aa[i] = "four"
  }
}


A$col = aa
head(A)

A$stat = gsub(" ","",A$stat)
mi = brewer.pal(9,"Set1")
mi= mi[3:6]

head(A)

plots = list()

Asub = A %>% filter(name == unique(A$name)[1])


p <- ggplot(Asub, aes(x = group, y = dd)) + geom_boxplot(alpha = 1, aes(fill = col)) +
  geom_jitter(position = position_jitter(0.17), size = 0.1, alpha = 0.5) +
  labs(x = "", y = "") + 
  geom_text(data = Asub, aes(x = group, y = y, label = stat),check_overlap = TRUE,vjust = -0.05,size = geom.text.size,position = position_dodge(width=0.9)) + 
  guides(color = guide_legend(title = NULL), shape = guide_legend(title = NULL)) + 
  facet_wrap(. ~ name, scales = "free_y", ncol = 1) +scale_fill_manual(values = mi) +
  scale_y_continuous(breaks = c(0.04,0.08,0.12),limits = c(0.04,0.12))
p
p = p  + scale_x_discrete(limits = axis_order) + mytheme1
p

plots[[1]] = p

i = 2
Asub = A %>% filter(name == unique(A$name)[i])
p <- ggplot(Asub, aes(x = group, y = dd)) + geom_boxplot(alpha = 1, aes(fill = col)) +
  geom_jitter(position = position_jitter(0.17), size = 0.1, alpha = 0.5) +
  labs(x = "", y = "") + 
  geom_text(data = Asub, aes(x = group, y = y, label = stat),check_overlap = TRUE,vjust = -0.05,size = geom.text.size,position = position_dodge(width=0.9)) + 
  guides(color = guide_legend(title = NULL), shape = guide_legend(title = NULL)) + 
  facet_wrap(. ~ name, scales = "free_y", ncol = 1) +scale_fill_manual(values = mi) +
  scale_y_continuous(breaks = c(0.5,1.5,2.5,3.5,4.5),limits = c(0.5,4.5))
p
p = p  + scale_x_discrete(limits = axis_order) + mytheme1
p
plots[[i]] = p

i = 3
Asub = A %>% filter(name == unique(A$name)[i])
p <- ggplot(Asub, aes(x = group, y = dd)) + geom_boxplot(alpha = 1, aes(fill = col)) +
  geom_jitter(position = position_jitter(0.17), size = 0.1, alpha = 0.5) +
  labs(x = "", y = "") + 
  geom_text(data = Asub, aes(x = group, y = y, label = stat),check_overlap = TRUE,vjust = -0.05,size = geom.text.size,position = position_dodge(width=0.9)) + 
  guides(color = guide_legend(title = NULL), shape = guide_legend(title = NULL)) + 
  facet_wrap(. ~ name, scales = "free_y", ncol = 1) +scale_fill_manual(values = mi) +
  scale_y_continuous(breaks = c(1,2,3),limits = c(1,3.5))
p
p = p  + scale_x_discrete(limits = axis_order) + mytheme1
p
plots[[i]] = p

library(ggpubr)
p  = ggarrange(plotlist = plots,  legend="right",ncol = 1)
p



dir.create("./peper")

ggsave("./peper//soil_chemi.pdf", p, width = 8, height = 18)


#---微生物组分析#-------
library(ggClusterNet)
library(phyloseq)
axis_order = c("Control","A","O","S","AO","AS","OS","AOS")

ps = readRDS("./ps.rds")
ps
path = "./peper//"
dir.create(path)
# 过滤样本

##--------
ps_sub <- subset_samples(ps,!Group %in% c("Y"));ps_sub 
ps_sub <- subset_samples(ps_sub,Group %in%c("CK","A","C","D","AC","AD","CD","ACD"));ps_sub
ps_sub = filter_taxa(ps_sub, function(x) sum(x ) > 1 , TRUE);ps_sub #筛选序列数量大于1的

map = read.delim("./Diff_TOCmap.txt",row.names = 1)
head(map)
# colnames(map) = c("ID","gro1","gro2","gro3","Group")
sample_data(ps_sub) = map

map = as.data.frame(sample_data(ps_sub))
head(map)
ps = ps_sub

#---统计基本微生物数据#------------

#step_1：
where <- "where"
method <- "using 16S rRNA gene amplicon sequencing."
rep = 9
step_1 <- paste("We analyzed the composition of microbial communities in the",where,method)


#step_2 统计测序样本总量每个样品中的序列数量
a = sum(sample_sums(ps))
b = "high-quality sequences were obtained"
# 统计样本数量
b1 = paste("across the" ,length(sample_sums(ps)),"samples",sep = " ")
b1
#统计重复数量
repead <- paste("For this analysis, we collected",rep,"repeats for each samples.",sep = " ")
repead 
#合并句子
each_count <- paste(repead,b1,a,b," and an average read count per sample of ",round(mean(sample_sums(ps)),0),"(standard deviation (SD)",round(sd(sample_sums(ps)),2),").",sep = " ")
each_count 

# step3 
aa = vegan_otu(ps)
otu_table = as.data.frame(t(aa))
otu_table = as.matrix(otu_table)
otu_table [otu_table > 0] <-1

otu_table = t(otu_table)

OTU_sum <- colSums(otu_table)


d = length(OTU_sum[OTU_sum > 0])


c = paste("All sequences were clustered into", d, "operational taxonomic units (OTUs).",sep = " ")

sample_tax <- paste(c,"the numbers of OTU, generally ranged between ",
                    min(OTU_sum)," and ",max(OTU_sum)," per sample with an average of ",
                    round(mean(OTU_sum),0),"(SD ",round(sd(OTU_sum)),")",sep = "")
sample_tax 

###step 4 统计门水平的总体丰度信息
library("tidyverse")
Taxonomies <- ps %>%
  tax_glom(taxrank = "Phylum") %>% 
  transform_sample_counts(function(x) {x/sum(x)} )%>% 
  psmelt() %>%
  #filter(Abundance > 0.05) %>%
  arrange(Phylum)
iris_groups<- group_by(Taxonomies, Phylum)
ps0_sum <- dplyr::summarise(iris_groups, mean(Abundance), sd(Abundance))
ps0_sum[is.na(ps0_sum)] <- 0
colnames(ps0_sum) = c("ID","mean","sd")
ps0_sum <- dplyr::arrange(ps0_sum,desc(mean))
ps0_sum$mean <- ps0_sum$mean *100
ps0_sum <- as.data.frame(ps0_sum)
a = paste(ps0_sum[1,1],"(",round(ps0_sum[1,2],2),"%"," with sd ",round(ps0_sum[1,3],3),")",sep = " ")
b = paste(ps0_sum[2,1],"(",round(ps0_sum[2,2],2),"%"," with sd ",round(ps0_sum[2,3],3),")",sep = " ")
c = paste(ps0_sum[3,1],"(",round(ps0_sum[3,2],2),"%"," with sd ",round(ps0_sum[3,3],3),")",sep = " ")
d = paste(ps0_sum[4,1],"(",round(ps0_sum[4,2],2),"%"," with sd ",round(ps0_sum[4,3],3),")",sep = " ")
e = paste(ps0_sum[5,1],"(",round(ps0_sum[5,2],2),"%"," with sd ",round(ps0_sum[5,3],3),")",sep = " ")
f = paste(ps0_sum[6,1],"(",round(ps0_sum[6,2],2),"%"," with sd ",round(ps0_sum[6,3],3),")",sep = " ")
g =  paste(ps0_sum[7,1],"(",round(ps0_sum[7,2],2),"%"," with sd ",round(ps0_sum[7,3],3),")",sep = " ")
h = paste(ps0_sum[8,1],"(",round(ps0_sum[8,2],2),"%"," with sd ",round(ps0_sum[8,3],3),")",sep = " ")
i = paste(ps0_sum[9,1],"(",round(ps0_sum[9,2],2),"%"," with sd ",round(ps0_sum[9,3],3),")",sep = " ")
j =  paste(ps0_sum[10,1],"(",round(ps0_sum[10,2],2),"%"," with sd ",round(ps0_sum[10,3],3),")",sep = " ")

tax_sum = paste("The majority of OTU belonged to the phyla",a,b,c,d,e,f,g,h,i,"and",j,".",sep = " ")
tax_sum
##all_first 
line = paste(step_1,each_count ,sample_tax ,tax_sum,sep = "")
line



path = "ACD"
ps_sub <- subset_samples(ps_sub,Group %in%c("Control","A","O","S"));ps_sub
ps_sub = filter_taxa(ps_sub, function(x) sum(x ) > 1 , TRUE);ps_sub #筛选序列数量大于1的

otu = as.data.frame(otu_table(ps))
dim(otu)
tax  = as.data.frame(vegan_tax(ps))
map = as.data.frame(sample_data(ps))


#----多样性图表#--------

dir.create(path)
alppath = paste(path,"/alpha/",sep = "")
dir.create(alppath)


source("G:\\Shared_Folder\\Function_local\\R_function\\micro/alpha-diversity.R")
index = c("Shannon","Inv_Simpson","Pielou_evenness","Simpson_evenness" ,"Richness" ,"Chao1","ACE" )


# 使用分面出图，得益于我开发的R包EasyAovWlxPlot
alp = alpha(otu = otu,tax = tax,map = map,inde="Shannon",group = "Group",Plot = TRUE )
index= alp[[3]]
head(index)
data = cbind(data.frame(ID = index$ID,group = index$Group),index[10:16])
head(data)
# #使用案例
library(EasyStat)
result = MuiaovMcomper(data = data,num = c(3,5,8),method_Mc = "Tukey")
result
result$Shannon = gsub(" ","",result$Shannon)
result$Pielou_evenness = gsub(" ","",result$Pielou_evenness)
result$Chao1 = gsub(" ","",result$Chao1)
result1 = FacetMuiPlotresultBox(data = data,num = c(3,5,8),result = result,sig_show ="abc",ncol = 1 )
p = result1[[1]] + theme_bw()+ scale_x_discrete(limits = axis_order)
p  

# FileName <- paste(alppath,"/alpha_diversity_all_abc.csv", sep = "")
# write.csv(result,FileName,sep = "")



A = result1[[2]]
aa = c()
i = 1
A$group


for (i in 1:length(A$group)) {
  if (A$group[i] == "Control" ) {
    aa[i] = "Control"
  }
  if (A$group[i] %in% c("A","B","O","S") ) {
    aa[i] = "single"
  }
  if (A$group[i] %in% c("AB","AO","AS","BO","BS","OS") ) {
    aa[i] = "two"
  }
  if (A$group[i] %in% c("ABO","ABS","AOS","BOS") ) {
    aa[i] = "three"
  }
  if (A$group[i] %in% c("ABOS") ) {
    aa[i] = "four"
  }
}

A$col = aa

A$stat = gsub(" ","",A$stat)
mi = brewer.pal(9,"Set1")
mi= mi[3:6]
head(A)

plots = list()

Asub = A %>% filter(name == unique(A$name)[1])


p <- ggplot(Asub, aes(x = group, y = dd)) + geom_boxplot(alpha = 1, aes(fill = col)) +
  geom_jitter(position = position_jitter(0.17), size = 0.1, alpha = 0.5) +
  labs(x = "", y = "") + 
  geom_text(data = Asub, aes(x = group, y = y, label = stat),check_overlap = TRUE,vjust = -0.05,size = geom.text.size,position = position_dodge(width=0.9)) + 
  guides(color = guide_legend(title = NULL), shape = guide_legend(title = NULL)) + 
  facet_wrap(. ~ name, scales = "free_y", ncol = 1) +scale_fill_manual(values = mi) +
  scale_y_continuous(breaks = c(2,3,4,5,6),limits = c(2.5,6))
p
p = p  + scale_x_discrete(limits = axis_order) + mytheme1
p

plots[[1]] = p

i = 2
Asub = A %>% filter(name == unique(A$name)[i])
p <- ggplot(Asub, aes(x = group, y = dd)) + geom_boxplot(alpha = 1, aes(fill = col)) +
  geom_jitter(position = position_jitter(0.17), size = 0.1, alpha = 0.5) +
  labs(x = "", y = "") + 
  geom_text(data = Asub, aes(x = group, y = y, label = stat),check_overlap = TRUE,vjust = -0.05,size = geom.text.size,position = position_dodge(width=0.9)) + 
  guides(color = guide_legend(title = NULL), shape = guide_legend(title = NULL)) + 
  facet_wrap(. ~ name, scales = "free_y", ncol = 1) +scale_fill_manual(values = mi) +
  scale_y_continuous(breaks = c(0.5,0.75,1.0),limits = c(0.5,1.0))
p
p = p  + scale_x_discrete(limits = axis_order) + mytheme1
p
plots[[i]] = p

i = 3
Asub = A %>% filter(name == unique(A$name)[i])
p <- ggplot(Asub, aes(x = group, y = dd)) + geom_boxplot(alpha = 1, aes(fill = col)) +
  geom_jitter(position = position_jitter(0.17), size = 0.1, alpha = 0.5) +
  labs(x = "", y = "") + 
  geom_text(data = Asub, aes(x = group, y = y, label = stat),check_overlap = TRUE,vjust = -0.05,size = geom.text.size,position = position_dodge(width=0.9)) + 
  guides(color = guide_legend(title = NULL), shape = guide_legend(title = NULL)) + 
  facet_wrap(. ~ name, scales = "free_y", ncol = 1) +scale_fill_manual(values = mi) +
  scale_y_continuous(breaks = c(250,750,1250),limits = c(240,1350))
p
p = p  + scale_x_discrete(limits = axis_order) + mytheme1
p
plots[[i]] = p

library(ggpubr)
p  = ggarrange(plotlist = plots,  legend="right",ncol = 1)
p


# p <- ggplot(A, aes(x = group, y = dd)) + geom_boxplot(alpha = 1, aes(fill = col),outlier.alpha = 0) +
#   # geom_jitter(position = position_jitter(0.17), size = 0.1, alpha = 0.5) +
#   labs(x = "", y = "") + 
#   geom_text(data = A, aes(x = group, y = y, label = stat),check_overlap = TRUE,vjust = 0.55,size = geom.text.size,position = position_dodge(width=0.9)) + 
#   guides(color = guide_legend(title = NULL), shape = guide_legend(title = NULL)) +
#   facet_wrap(. ~ name, scales = "free_y", ncol = 1) +scale_fill_manual(values = mi)
# 
# p
# 
# p = p  + scale_x_discrete(limits = axis_order) + mytheme1
# p


FileName <- paste(alppath,"Facet_box", ".pdf", sep = "_")
ggsave(FileName, p, width = 8, height = 18)

#-----beta排序#---------


betapath = paste(path,"/beta/",sep = "")
dir.create(betapath)
ps

library(phyloseq)
library(vegan)


source("G:\\Shared_Folder\\Function_local\\R_function\\micro/BetaDiv.R")
source("G:\\Shared_Folder\\Function_local\\R_function\\micro/MicroTest.R")
source("G:\\Shared_Folder\\Function_local\\R_function\\micro/pairMicroTest.R")

library(phyloseq)
result = BetaDiv(ps = ps, group = "Group", dist = "bray", method = "PCoA", Micromet = "adonis", pvalue.cutoff = 0.05)

#不带标签默认出图
mi = brewer.pal(9,"Set1")
p1 = result[[1]] + mytheme1 + scale_color_manual(values = mi,guide = FALSE)  +
  scale_fill_manual(values = mi)
  
p1


#---修改颜色和形状

FileName <- paste(betapath,"/a2_bray_PCOA.pdf", sep = "")
ggsave(FileName, p1, width = 10, height = 8)
# FileName1 <- paste(betapath,"/a2_bray_PCOA.jpg", sep = "")
# ggsave(FileName1 , p1, width = 12, height = 8)



###微生物门类展示#---------

# BiocManager::install("ggalluvial")
# devtools::install_github("taowenmicro/EasyStat",force = TRUE)
library(EasyMicrobiome)
#-------门类水平展示
source("G:\\Shared_Folder\\Function_local\\R_function\\micro/barMainplot.R")
barpath = paste(path,"/barplot/",sep = "")
dir.create(barpath)
j = "Phylum"
result = barMainplot(ps = ps,j = "Phylum",rep = 9,axis_ord = NULL,label = FALSE ,sd = FALSE,Top = 10)

p4 =result[[1]]  + scale_fill_brewer(palette = "Paired") +  scale_x_discrete(limits = axis_order) + mytheme1
p4
p3 =result[[3]]  + scale_fill_brewer(palette = "Paired") +  scale_x_discrete(limits = axis_order) + mytheme1
p3

FileName1 <- paste(barpath,"/a2_",j,"_barflow",".pdf", sep = "")
ggsave(FileName1, p3, width = 12, height =8 )
FileName2 <- paste(barpath,"/a2_",j,"_barflow",".jpg", sep = "")
ggsave(FileName2, p3, width = 12, height =8 )

FileName1 <- paste(barpath,"/a2_",j,"_bar",".pdf", sep = "")
ggsave(FileName1, p4, width = 12, height =8 )
FileName2 <- paste(barpath,"/a2_",j,"_bar",".jpg", sep = "")
ggsave(FileName2, p4, width = 12, height =8 )



#----------附图环境一因子联合#------
library(EasyStat)
source("G:\\Shared_Folder\\Function_local\\R_function\\micro/rda-cca.R")

# env = read.csv("./env/env5change_nine.csv",row.names = 1)
env = read.csv("./env/env5_NEW_change_nine.csv",row.names = 1)
head(env)
# env$group = NULL
row.names(env) = env$ID
# env = data.frame(row.names = env$ID,Fe = env$Fe,Al = env$Al,Si = env$Si,Toc = env$TOCnew)
env = data.frame(row.names = env$ID,AFe = env$Fe,ASi = env$Si,DOC = env$TOCnew)

rdapath  = paste(path,"/rda_cca/",sep = "")
dir.create(rdapath)
result = RDA_CCA(ps = ps,env = env,path = path )
#提取图片
p = result[[1]]  + scale_fill_brewer(palette = "Paired")  + mytheme1
p
# 提取作图数据
dataplot = result[[2]]
# 提取带有标记的图片
p3 = result[[3]]

#提取理化提供群落关系的渐检验结果
aov = result[[4]]

##保存
plotnamea = paste(rdapath,"/RDA_envnew.pdf",sep = "")
ggsave(plotnamea, p, width = 10, height = 8)
plotnamea4 = paste(rdapath,"/RDA_envnew.jpg",sep = "")
ggsave(plotnamea4, p, width = 8, height = 6)


filenamea = paste(rdapath,"dataplotnew.txt",sep = "")
write.table(dataplot ,file=filenamea,sep="\t",col.names=NA)

filenamea = paste(rdapath,"aovnew.txt",sep = "")
write.table(aov,file=filenamea,sep="\t",col.names=NA)

plotnamea = paste(rdapath,"/RDA_envlabelnew.pdf",sep = "")
ggsave(plotnamea, p3, width = 18, height = 12)
plotnamea4 = paste(rdapath,"/RDA_envlabelnew.jpg",sep = "")
ggsave(plotnamea4, p3, width = 18, height = 12)


#---PCoA排序DOC分组#----

betapath = paste(path,"/betaDOC/",sep = "")
dir.create(betapath)


source("G:\\Shared_Folder\\Function_local\\R_function\\micro/BetaDiv.R")
source("G:\\Shared_Folder\\Function_local\\R_function\\micro/MicroTest.R")
source("G:\\Shared_Folder\\Function_local\\R_function\\micro/pairMicroTest.R")
head(map)
ps_sub <- subset_samples(ps,!Group %in% c("Y","CK","Control"));ps_sub

library(phyloseq)
result = BetaDiv(ps = ps_sub, group = "Group4", dist = "bray", method = "PCoA", Micromet = "adonis", pvalue.cutoff = 0.05)


p2 = result[[1]] + mytheme1 + scale_color_manual(values = mi,guide = FALSE)  +
  scale_fill_manual(values = mi)
p2

#---------精修图

plotdata <- result[[2]]
head(plotdata)


# 求均值
cent <- aggregate(cbind(x,y) ~Group, data = plotdata, FUN = mean)
cent
# 合并到样本坐标数据中
segs <- merge(plotdata, setNames(cent, c('Group','oNMDS1','oNMDS2')),
              by = 'Group', sort = FALSE)


p2$layers[[2]] = NULL
# library(ggcor)
library(ggsci)

p = p2 +geom_segment(data = segs,
                     mapping = aes(xend = oNMDS1, yend = oNMDS2,color = Group)) + # spiders
  geom_point(mapping = aes(x = x, y = y),data = cent, size = 5,pch = 24,color = "black",fill = "yellow") +
  scale_fill_lancet()+
  scale_color_lancet() 


p

FileName <- paste(betapath,"/a2_bray_PCOA.pdf", sep = "")
ggsave(FileName, p, width = 10, height = 8)


#------单个PH的可视化#--------
data = read.csv("./env/env5_NEW_change_nine.csv",row.names = 1)
head(data)

axis_order = c("Control","A","O","S","AO","AS","OS","AOS")

data$group = gsub("CK","Control",data$group)
data$group = gsub("Y","Raw soil",data$group)
data$group = as.factor(data$group)


data = data %>% filter(group %in% axis_order)
data$group = as.factor(data$group)
data = data[,c(1:2,6)]
head(data)
colnames(data)[3] = c("pH")

result = MuiaovMcomper(data = data,num = c(3))
result

result <- ord_sig(data = result,ID = NULL)
result$pH = gsub(" ","",result$pH)

data = data[data$group %in% axis_order,]
result = as.data.frame(result)
# result = result[match(axis_order,row.names(result)),]

result1 = FacetMuiPlotresultBox(data = data,num = c(3),result = result,sig_show ="abc",ncol = 1 )
p = result1[[1]] + theme_bw() + scale_x_discrete(limits = axis_order)
p

A = result1[[2]]

head(A)

aa = c()
i = 1
A$group


for (i in 1:length(A$group)) {
  if (A$group[i] == "Control" ) {
    aa[i] = "Control"
  }
  if (A$group[i] %in% c("A","B","O","S") ) {
    aa[i] = "single"
  }
  if (A$group[i] %in% c("AB","AO","AS","BO","BS","OS") ) {
    aa[i] = "two"
  }
  if (A$group[i] %in% c("ABO","ABS","AOS","BOS") ) {
    aa[i] = "three"
  }
  if (A$group[i] %in% c("ABOS") ) {
    aa[i] = "four"
  }
}

A$col = aa


p <- ggplot(A, aes(x = group, y = dd)) + geom_boxplot(alpha = 1, aes(fill = col)) +
  geom_jitter(position = position_jitter(0.17), size = 0.1, alpha = 0.5) +
  labs(x = "", y = "") + 
  geom_text(data = A, aes(x = group, y = y, label = stat),check_overlap = TRUE,vjust = -0.05,size = geom.text.size,position = position_dodge(width=0.9)) + 
  guides(color = guide_legend(title = NULL), shape = guide_legend(title = NULL)) + 
  facet_wrap(. ~ name, scales = "free_y", ncol = 1) +scale_fill_manual(values = mi)

p

p = p  + scale_x_discrete(limits = axis_order) + mytheme1
p


ggsave("./peper//PH.pdf", p, width = 10, height =8)

#---单组添加物质的土壤化学性质可视化#--------

axis_order = c("Control","A","O","S")

data = read.csv("./env//env5_NEW_change_nine.csv",row.names = 1)
head(data)
tail(data,20)

data = data[,c(1:2,5,7,10)]
head(data)
colnames(data)[3:5] = c("AFe(ug/g)","ASi(ug/g)","DOC(mg/kg)")

#--------
data = dplyr::filter(data, group %in% c("A","O","S","Control"))
data$group = as.factor(data$group)
result = MuiaovMcomper(data = data,num = c(3:5))
result

result <- ord_sig(data = result,ID = NULL)
result$AFe = gsub(" ","",result$AFe)
result$ASi = gsub(" ","",result$ASi)
result$DOC = gsub(" ","",result$DOC)


result1 = FacetMuiPlotresultBox(data = data,num = c(3:5),result = result,sig_show ="abc",ncol = 1 )
p = result1[[1]] + theme_bw()+ scale_x_discrete(limits = axis_order)
p

A = result1[[2]]

head(A)
# A$group = gsub("CK","Control",A$group)
aa = c()
i = 1
A$group

for (i in 1:length(A$group)) {
  if (A$group[i] == "Control" ) {
    aa[i] = "Control"
  }
  if (A$group[i] %in% c("A","B","O","S") ) {
    aa[i] = "single"
  }
  if (A$group[i] %in% c("AB","AO","AS","BO","BS","OS") ) {
    aa[i] = "two"
  }
  if (A$group[i] %in% c("ABO","ABS","AOS","BOS") ) {
    aa[i] = "three"
  }
  if (A$group[i] %in% c("ABOS") ) {
    aa[i] = "four"
  }
}

A$col = aa
head(A)

p <- ggplot(A, aes(x = group, y = dd)) + geom_boxplot(alpha = 1, aes(fill = group)) +
  geom_jitter(position = position_jitter(0.17), size = 0.1, alpha = 0.5) +
  labs(x = "", y = "") + 
  geom_text(data = A, aes(x = group, y = y, label = stat),check_overlap = TRUE,vjust = -0.05,size = geom.text.size,position = position_dodge(width=0.9)) + 
  guides(color = guide_legend(title = NULL), shape = guide_legend(title = NULL)) + 
  facet_wrap(. ~ name, scales = "free_y", ncol = 1) +scale_fill_manual(values = mi)

p

p = p  + scale_x_discrete(limits = axis_order) + mytheme1
p


ggsave("./peper/AOS_single.pdf", p, width = 8, height =18)

#---相关图表的修改#-

unique(map$Group)

unique(map$Group[map$Group4 == "Down"])
unique(map$Group[map$Group4 == "Up"])
unique(map$Group[map$Group4 == "Nose"])


#----------ACD可视化#--------------- 
path = "ACD_single"
ps_sub <- subset_samples(ps,Group %in%c("Control","A","O","S"));ps_sub
ps_sub = filter_taxa(ps_sub, function(x) sum(x ) > 1 , TRUE);ps_sub #筛选序列数量大于1的



otu = as.data.frame(otu_table(ps_sub))
dim(otu)
tax  = as.data.frame(vegan_tax(ps_sub))
map = as.data.frame(sample_data(ps_sub))


#----多样性图表#--------

dir.create(path)
alppath = paste(path,"/alpha/",sep = "")
dir.create(alppath)

axis_order = c("Control","A","O","S")

source("G:\\Shared_Folder\\Function_local\\R_function\\micro/alpha-diversity.R")
index = c("Shannon","Inv_Simpson","Pielou_evenness","Simpson_evenness" ,"Richness" ,"Chao1","ACE" )


# 使用分面出图，得益于我开发的R包EasyAovWlxPlot
alp = alpha(otu = otu,tax = tax,map = map,inde="Shannon",group = "Group",Plot = TRUE )
index= alp[[3]]
head(index)
data = cbind(data.frame(ID = index$ID,group = index$Group),index[10:16])
head(data)
# #使用案例
library(EasyStat)
result = MuiaovMcomper(data = data,num = c(3,5,8),method_Mc = "Tukey")
result
result$Shannon = gsub(" ","",result$Shannon)
result$Pielou_evenness = gsub(" ","",result$Pielou_evenness)
result$Chao1 = gsub(" ","",result$Chao1)
result1 = FacetMuiPlotresultBox(data = data,num = c(3,5,8),result = result,sig_show ="abc",ncol = 1 )
p = result1[[1]] + theme_bw()+ scale_x_discrete(limits = axis_order)
p  
aovMuiBoxP
# FileName <- paste(alppath,"/alpha_diversity_all_abc.csv", sep = "")
# write.csv(result,FileName,sep = "")



A = result1[[2]]
A

A <- ord_sig(data = A,ID = NULL)

aa = c()
i = 1
A$group


for (i in 1:length(A$group)) {
  if (A$group[i] == "Control" ) {
    aa[i] = "Control"
  }
  if (A$group[i] %in% c("A","B","O","S") ) {
    aa[i] = "single"
  }
  if (A$group[i] %in% c("AB","AO","AS","BO","BS","OS") ) {
    aa[i] = "two"
  }
  if (A$group[i] %in% c("ABO","ABS","AOS","BOS") ) {
    aa[i] = "three"
  }
  if (A$group[i] %in% c("ABOS") ) {
    aa[i] = "four"
  }
}

A$col = aa

A$stat = gsub(" ","",A$stat)
mi = brewer.pal(9,"Set1")
# mi= mi[3:6]

p <- ggplot(A, aes(x = group, y = dd)) + geom_boxplot(alpha = 1, aes(fill = col),outlier.alpha = 0) +
  # geom_jitter(position = position_jitter(0.17), size = 0.1, alpha = 0.5) +
  labs(x = "", y = "") + 
  geom_text(data = A, aes(x = group, y = y, label = stat),check_overlap = TRUE,vjust = 0.55,size = geom.text.size,position = position_dodge(width=0.9)) + 
  guides(color = guide_legend(title = NULL), shape = guide_legend(title = NULL)) +
  facet_wrap(. ~ name, scales = "free_y", ncol = 1) +scale_fill_manual(values = mi)

p

p = p  + scale_x_discrete(limits = axis_order) + mytheme1
p


FileName <- paste(alppath,"Facet_box", ".pdf", sep = "_")
ggsave(FileName, p, width = 6, height = 10)

#-----beta排序#---------


betapath = paste(path,"/beta/",sep = "")
dir.create(betapath)
ps

library(phyloseq)
library(vegan)


source("G:\\Shared_Folder\\Function_local\\R_function\\micro/BetaDiv.R")
source("G:\\Shared_Folder\\Function_local\\R_function\\micro/MicroTest.R")
source("G:\\Shared_Folder\\Function_local\\R_function\\micro/pairMicroTest.R")

library(phyloseq)
result = BetaDiv(ps = ps_sub, group = "Group", dist = "bray", method = "PCoA", Micromet = "adonis", pvalue.cutoff = 0.05)

#不带标签默认出图
mi = brewer.pal(9,"Set1")
p1 = result[[1]] + mytheme1 + scale_color_manual(values = mi,guide = FALSE)  +
  scale_fill_manual(values = mi)

p1


#---修改颜色和形状

FileName <- paste(betapath,"/a2_bray_PCOA.pdf", sep = "")
ggsave(FileName, p1, width = 10, height = 8)
# FileName1 <- paste(betapath,"/a2_bray_PCOA.jpg", sep = "")
# ggsave(FileName1 , p1, width = 12, height = 8)

#---flower plot#-----

flowpath = paste(path,"/flowplot/",sep = "")
dir.create(flowpath)
# source("G:\\Shared_Folder\\Function_local\\R_function\\micro/flowerplot.R")
# flowerplot(ps = ps,rep = 6,path =flowpath )
source("G:\\Shared_Folder\\Function_local\\R_function\\micro/ggflowerplot.R")
p0_1 <- ggflower(ps = ps_sub,
                 rep = 9,
                 group = "Group",
                 start = 1, # 风车效果
                 m1 = 2, # 花瓣形状，方形到圆形到棱形，数值逐渐减少。
                 a = 0.2, # 花瓣胖瘦
                 b = 1, # 花瓣距离花心的距离
                 lab.leaf = 1, # 花瓣标签到圆心的距离
                 col.cir = "yellow"
)
# p + scale_fill_brewer(palette = "Paired")
p0_2 <- p0_1
p0_2

FileName1 <- paste(flowpath,"ggflower.pdf", sep = "")
ggsave(FileName1, p0_1, width = 4, height = 4)
FileName2 <- paste(flowpath,"ggflower.jpg", sep = "")
ggsave(FileName2, p0_1, width = 4, height = 4 )


###微生物门类展示#---------

# BiocManager::install("ggalluvial")
# devtools::install_github("taowenmicro/EasyStat",force = TRUE)
library(EasyMicrobiome)
#-------门类水平展示
source("G:\\Shared_Folder\\Function_local\\R_function\\micro/barMainplot.R")
barpath = paste(path,"/barplot/",sep = "")
dir.create(barpath)
j = "Phylum"
result = barMainplot(ps = ps_sub,j = "Phylum",rep = 9,axis_ord = NULL,label = FALSE ,sd = FALSE,Top = 10)

p4 =result[[1]]  + scale_fill_brewer(palette = "Paired") +  scale_x_discrete(limits = axis_order) + mytheme1
p4
p3 =result[[3]]  + scale_fill_brewer(palette = "Paired") +  scale_x_discrete(limits = axis_order) + mytheme1
p3

FileName1 <- paste(barpath,"/a2_",j,"_barflow",".pdf", sep = "")
ggsave(FileName1, p3, width = 10, height =8 )
FileName2 <- paste(barpath,"/a2_",j,"_barflow",".jpg", sep = "")
ggsave(FileName2, p3, width = 10, height =8 )

FileName1 <- paste(barpath,"/a2_",j,"_bar",".pdf", sep = "")
ggsave(FileName1, p4, width = 10, height =8 )
FileName2 <- paste(barpath,"/a2_",j,"_bar",".jpg", sep = "")
ggsave(FileName2, p4, width = 10, height =8 )



#----------附图环境一因子联合#------
library(EasyStat)
source("G:\\Shared_Folder\\Function_local\\R_function\\micro/rda-cca.R")

# env = read.csv("./env/env5change_nine.csv",row.names = 1)
env = read.csv("./env/env5_NEW_change_nine.csv",row.names = 1)
head(env)
# env$group = NULL
row.names(env) = env$ID
# env = data.frame(row.names = env$ID,Fe = env$Fe,Al = env$Al,Si = env$Si,Toc = env$TOCnew)
env = data.frame(row.names = env$ID,AFe = env$Fe,ASi = env$Si,DOC = env$TOCnew)

rdapath  = paste(path,"/rda_cca/",sep = "")
dir.create(rdapath)
result = RDA_CCA(ps = ps_sub,env = env,path = path )
#提取图片
p = result[[1]]  + scale_fill_brewer(palette = "Paired")  + mytheme1
p
# 提取作图数据
dataplot = result[[2]]
# 提取带有标记的图片
p3 = result[[3]]

#提取理化提供群落关系的渐检验结果
aov = result[[4]]

##保存
plotnamea = paste(rdapath,"/RDA_envnew.pdf",sep = "")
ggsave(plotnamea, p, width = 10, height = 8)
plotnamea4 = paste(rdapath,"/RDA_envnew.jpg",sep = "")
ggsave(plotnamea4, p, width = 8, height = 6)


# filenamea = paste(rdapath,"dataplotnew.txt",sep = "")
# write.table(dataplot ,file=filenamea,sep="\t",col.names=NA)
# 
# filenamea = paste(rdapath,"aovnew.txt",sep = "")
# write.table(aov,file=filenamea,sep="\t",col.names=NA)
# 
# plotnamea = paste(rdapath,"/RDA_envlabelnew.pdf",sep = "")
# ggsave(plotnamea, p3, width = 18, height = 12)
# plotnamea4 = paste(rdapath,"/RDA_envlabelnew.jpg",sep = "")
# ggsave(plotnamea4, p3, width = 18, height = 12)


#---PCoA排序DOC分组#----

betapath = paste(path,"/betaDOC/",sep = "")
dir.create(betapath)


source("G:\\Shared_Folder\\Function_local\\R_function\\micro/BetaDiv.R")
source("G:\\Shared_Folder\\Function_local\\R_function\\micro/MicroTest.R")
source("G:\\Shared_Folder\\Function_local\\R_function\\micro/pairMicroTest.R")
head(map)
ps_sub <- subset_samples(ps,!Group %in% c("Y","CK","Control"));ps_sub

library(phyloseq)
result = BetaDiv(ps = ps_sub, group = "Group4", dist = "bray", method = "PCoA", Micromet = "adonis", pvalue.cutoff = 0.05)


p2 = result[[1]] + mytheme1 + scale_color_manual(values = mi,guide = FALSE)  +
  scale_fill_manual(values = mi)
p2

#---------精修图

plotdata <- result[[2]]
head(plotdata)


# 求均值
cent <- aggregate(cbind(x,y) ~Group, data = plotdata, FUN = mean)
cent
# 合并到样本坐标数据中
segs <- merge(plotdata, setNames(cent, c('Group','oNMDS1','oNMDS2')),
              by = 'Group', sort = FALSE)


p2$layers[[2]] = NULL
# library(ggcor)
library(ggsci)

p = p2 +geom_segment(data = segs,
                     mapping = aes(xend = oNMDS1, yend = oNMDS2,color = Group)) + # spiders
  geom_point(mapping = aes(x = x, y = y),data = cent, size = 5,pch = 24,color = "black",fill = "yellow") +
  scale_fill_lancet()+
  scale_color_lancet() 


p

FileName <- paste(betapath,"/a2_bray_PCOA.pdf", sep = "")
ggsave(FileName, p, width = 10, height = 8)


#----距离图表绘制#----------

#--读入文件
library(phyloseq)
library(vegan)
ps = readRDS("./ps.rds")

#--一下调整必须在map文件中包含两列，一列是时间梯度，使用数字表示，不能带有特殊表述，一列是treat，就是分组信息，这些分组信息必须是有时间梯度的。
map = as.data.frame(sample_data(ps))
head(map)


tre1 = c("Control","A","AO","AS","AOS")
tre2 = c("Control","O","AO","OS","AOS")
tre3 = c("Control","S","AS","OS","AOS")
datab  = data.frame(tre1 = tre1,tre2 = tre2,tre= tre3)

refer = "Control"

library(tidyverse)
plot = list()
for (i in 1:3) {
  tre = as.character(datab[[i]])
  map = as.data.frame(sample_data(ps))
  maps<- dplyr::filter(as.tibble(map),Group%in% tre)
  row.names(maps) = maps$ID.1
  ps_sub = ps
  sample_data( ps_sub ) = as.data.frame(maps);ps_sub 
  
  
  ps_rela  = transform_sample_counts(ps_sub , function(x) x / sum(x) );ps_rela 
  
  #提取otu表格，筛选数据
  otu= as.data.frame(t(vegan_otu(ps_rela)))
  
  head(otu)
  #-添加第一个时间的均值作为标准，求取均值列加入其后
  mapsub = as.data.frame(sample_data(ps_rela))
  #挑选time为1的样本
  ID = mapsub$ID[mapsub$Group == refer]
  #挑选样本
  head(otu)
  ID = as.character(ID)
  index = otu[,ID] # 筛选并计算均值
  head(index,n = 5)
  
  index_mean= rowMeans(index)
  otu$index_mean=index_mean # 均值添加到OTU表
  head(otu)
  # 计算包括终点均值的所有样品bray距离
  bray_curtis = vegan::vegdist(t(otu), method = "bray")
  bray_curtis= as.matrix(bray_curtis)
  
  #--------------整理数据-
  # 计算各组内部差异程度
  # 建立一个存储组内差异的数据框
  dat = t(as.data.frame(c("sampleA","sampleB","0.15","group","genosite")))
  colnames(dat) = c("sampleA","sampleB","distance","group","type")
  rownames(dat) = c("test")
  
  
  # 每个样品与final对应的距离
  # ii = 1
  for (ii in tre){
    group = rownames(mapsub[mapsub$Group %in% ii,])
    # m = 1
    for (m in 1:(length(group))) {
      x = c(group[m],"index_mean",bray_curtis[group[m],"index_mean"],ii,tre)
      dat=rbind(dat,x)
    }
  }
  dat = as.data.frame(dat[-1,], stringsAsFactors=F) # 删除首行框架数据
  
  # dat = dat[dat$distance != 0,] # 
  # 距离转换为3位数值
  dat$distance=round(as.numeric(as.character(dat$distance)), digits=3)
  # 分组添加levels，分组有时间顺序
  dat$group = factor(dat$group, levels=unique(dat$group))
  
  all = dat
  # 开始尝试散点图拟合
  require(splines)
  require(MASS)
  library(ggplot2)
  
  all1 = all
  all1$group = as.numeric(all1$group)
  
  head(all1)
  library(dplyr)
  iris_groups<- group_by(all, group)
  all3<- dplyr::summarise(iris_groups, mean(distance))
  colnames(all3) = c("group","mean")
  datawt = data.frame(ID  =all$sampleA,group = all$group,distance = all$distance)
  library(EasyStat)
  result= aovMcomper (data = datawt, i= 3,method_Mc = "Tukey")
  # 提取多重比较结果
  result[[1]]
  
  PlotresultBox = aovMuiBoxP(data = datawt, i= 3,sig_show ="abc",result = result[[1]])
  #提取图片
  p = PlotresultBox[[1]]
  p = p + geom_line(aes(x=group, y=mean,group = 1),data = all3) +
    geom_point(aes(x=group, y=mean),all3,alpha=1, size=4,pch = 21,fill ="#999999" )
  p <- p +  mytheme1 + scale_fill_manual(values = mi)
  plot[[i]] = p
  
}


library(ggpubr)
p  = ggarrange(plotlist = plot, common.legend = TRUE, legend="right", ncol = 3)
p

dispath = paste(path,"/disance/",sep = "")
dir.create(dispath)

FileName <- paste(dispath,"distance", ".pdf", sep = "_")
ggsave(FileName, p, width = 18, height =6)


#-------fig 5 #-----------
ps = readRDS("./ps.rds")
ps
path = "./peper/spe_OTU_cor//"
dir.create(path)
# 过滤样本

##--------
# 取子集
ps_sub <- subset_samples(ps,!Group %in% c("Y","CK"));ps_sub 
ps_sub <- subset_samples(ps,!Group %in% c("Y"));ps_sub 
ps_sub <- subset_samples(ps_sub,Group %in%c("CK","A","C","D","AC","AD","CD","ACD"));ps_sub 
map = read.delim("./Diff_TOCmap.txt",row.names = 1)
head(map)
colnames(map) = c("ID","gro1","gro2","gro3","Group")
sample_data(ps_sub) = map

map = as.data.frame(sample_data(ps_sub))
head(map)
ps = ps_sub
ps

otu = as.data.frame(otu_table(ps))
dim(otu)
tax  = as.data.frame(vegan_tax(ps))
map = as.data.frame(sample_data(ps))

#--第二次分析，CK不要了
ps_sub <- subset_samples(ps,!Group %in% c("Y","CK"));ps_sub
ps = ps_sub
ps

#-----------------特定菌-可视化丰度#---------------
library(EasyMicrobiome)
library(edgeR)

ps
otu_table = as.data.frame(vegan_otu(ps))

count = as.matrix(otu_table)
count <- t(count)
#数据整理形式
sub_design <- as.data.frame(sample_data(ps))
dim(sub_design)
sub_design$SampleType = as.character(sub_design$Group)
sub_design$SampleType <- as.factor(sub_design$Group)

# create DGE list
d = DGEList(counts=count, group=sub_design$SampleType)
d = calcNormFactors(d)#默认为TMM标准化

# 生成实验设计矩阵
design.mat = model.matrix(~ 0 + d$samples$group)
colnames(design.mat)=levels(sub_design$SampleType)
d2 = estimateGLMCommonDisp(d, design.mat)
d2 = estimateGLMTagwiseDisp(d2, design.mat)
fit = glmFit(d2, design.mat)

otu = as.matrix(fit$fitted.values)
dim(otu)
#-组合phyloseq对象
ps_sub1 <- phyloseq(otu_table(otu, taxa_are_rows=TRUE),
                    sample_data(ps_sub),
                    tax_table(ps_sub)
                    
)
ps_sub1

ps_sub <- ps_sub1 %>%
  subset_taxa(
    row.names(tax_table(ps_sub1 ))%in% c(
      # "OTU_14914",
      # "OTU_4523",
      # "OTU_3414",
      "OTU_192",
      "OTU_44",
      "OTU_470",
      "OTU_347",
      "OTU_3583",
      "OTU_1812",
      "OTU_245",
      "OTU_3357",
      "OTU_383",
      "OTU_537",
      "OTU_677")
  )
otu = as.data.frame(t(otu_table(ps_sub)))
head(otu)
map = as.data.frame(sample_data(ps_sub))
map = map[,5]
otuGroup = merge(otu,map,by= "row.names",all = F)
head(otuGroup)
colnames(otuGroup)[1] = "ID"
library("reshape2")
datap = melt(otuGroup, id=c("ID","Group"))
head(datap)
library(ggplot2)
p = ggplot() + geom_boxplot(data = datap,aes(x = Group,y = value,fill = Group)) +
  facet_wrap(~ variable,scales = "free",ncol = 4) 
p  = p +scale_fill_manual(values = mi) + mytheme1
p



ggsave("./peper/spe_OTU_cor//enrich_spe_OTU.pdf", p, width = 12, height = 11)

#---相关-enrich-#-----------

otuenv = merge(otu,env,by = "row.names",all = F)
head(otuenv)

otuenv = otuenv[,c(1:12,15)]
colnames(otuenv)[13] = "DOC"

row.names(otuenv) = otuenv$Row.names
otuenv$Row.names = NULL
otuenv$ID = row.names(otuenv)


require(ggplot2)
require(reshape2)


map = as.data.frame(sample_data(ps))
map = data.frame(ID = map$ID,group = map$Group)
head(map)
library(dplyr)

mtcars2 = melt(otuenv, id.vars=c('DOC',"ID"))
mtcars2 = mtcars2 %>% left_join(map, by = "ID")
head(mtcars2)

i = 1
plots = list()
library("ggExtra")

for (i in 1:length(unique(mtcars2$variable))) {
  mtcars3 <- filter(mtcars2 , variable == as.character(unique(mtcars2$variable)[i]))
  head(mtcars3)
  
  # library(dplyr)
  # iris_groups<- group_by(mtcars3, variable)
  # lab <- dplyr::summarise(iris_groups, max(DOC,na.rm = T), max(value,na.rm = T))
  # 
  lab <- mtcars3 %>%
    group_by(variable) %>%
    summarise(max(DOC), max(value))
  
  colnames(lab) = c("ID","DOC","value")
  
  lab = min(mtcars3$DOC)
  # ?stat_regline_equation
  p = ggplot2::ggplot() + 
    ggplot2::geom_point(data = mtcars3,aes(value,DOC,color = group)) + 
    ggpubr::stat_cor(data = mtcars3,aes(x = value,y = DOC),label.y=lab*0.85)+
    ggpubr::stat_regline_equation(data = mtcars3,aes(x = value,y = DOC),label.y=lab) +
    # facet_wrap(~variable, scales="free_x") +
    geom_smooth(data = mtcars3,aes(value,DOC), method=lm, se=T) +
    labs(x = as.character(unique(mtcars2$variable)[i])) + guides(color = F)
  p 
  
  p = p + scale_color_manual(values = mi) + mytheme1
  
  p = ggMarginal(p, groupColour = TRUE, groupFill = TRUE)  
  p
  plots[[i]] = p
}

#----使用list格式组合图形
library(ggpubr)
p  = ggarrange(plotlist = plots,ncol = 4,nrow = 3)
p

ggsave("./peper/spe_OTU_cor/cor_enrich.pdf",width = 12,height = 9)



#-- 另一批特定菌丰度展示
#---可能固定碳的卫微生物类群
library(tidyr)
ps_sub <- ps_sub1 %>%
  subset_taxa(
    row.names(tax_table(ps_sub1 ))%in% c(
      "OTU_50",
      "OTU_78",
      "OTU_149",
      "OTU_177",
      "OTU_285",
      "OTU_141",
      "OTU_333",
      "OTU_9170"
    )
  )
ps_sub

otu = as.data.frame(t(otu_table(ps_sub)))
head(otu)
map = as.data.frame(sample_data(ps_sub))
map = map[,5]
otuGroup = merge(otu,map,by= "row.names",all = F)
head(otuGroup)
colnames(otuGroup)[1] = "ID"
datap = melt(otuGroup, id=c("ID","Group"))
head(datap)
p = ggplot() + geom_boxplot(data = datap,aes(x = Group,y = value,fill = Group)) +
  facet_wrap(~ variable,scales = "free",ncol = 4) 
p 
p  = p +scale_fill_manual(values = mi) + mytheme1
p

ggsave("./peper/spe_OTU_cor//down_spe_OTU.pdf", p, width = 12, height = 8)



#--------相关-固碳#-------
require(ggplot2)
require(reshape2)
library(dplyr)

otuenv = merge(otu,env,by = "row.names",all = F)
head(otuenv)
otuenv = otuenv[,c(1:9,12)]
colnames(otuenv)[10] = "DOC"
row.names(otuenv) = otuenv$Row.names
otuenv$Row.names = NULL
otuenv$ID = row.names(otuenv)

map = as.data.frame(sample_data(ps))
map = data.frame(ID = map$ID,group = map$Group)
head(map)


mtcars2 = melt(otuenv, id.vars=c('DOC',"ID"))
mtcars2 = mtcars2 %>% left_join(map, by = "ID")
head(mtcars2)

i = 1
plots = list()
library("ggExtra")
for (i in 1:length(unique(mtcars2$variable))) {
  mtcars3 <- filter(mtcars2 , variable == as.character(unique(mtcars2$variable)[i]))
  head(mtcars3)
  
  # library(dplyr)
  # iris_groups<- group_by(mtcars3, variable)
  # lab <- dplyr::summarise(iris_groups, max(DOC,na.rm = T), max(value,na.rm = T))
  # 
  lab <- mtcars3 %>%
    group_by(variable) %>%
    summarise(max(DOC), max(value))
  
  colnames(lab) = c("ID","DOC","value")
  
  lab = min(mtcars3$DOC)
  # ?stat_regline_equation
  p = ggplot2::ggplot() + 
    ggplot2::geom_point(data = mtcars3,aes(value,DOC,color = group)) + 
    ggpubr::stat_cor(data = mtcars3,aes(x = value,y = DOC),label.y=lab*0.85)+
    ggpubr::stat_regline_equation(data = mtcars3,aes(x = value,y = DOC),label.y=lab) +
    # facet_wrap(~variable, scales="free_x") +
    geom_smooth(data = mtcars3,aes(value,DOC), method=lm, se=T) +
    labs(x = as.character(unique(mtcars2$variable)[i])) + guides(color = F)
  p 
  
  p = p + scale_color_manual(values = mi) + mytheme1
  
  p = ggMarginal(p, groupColour = TRUE, groupFill = TRUE)  
  p
  plots[[i]] = p
}

#----使用list格式组合图形
p  = ggarrange(plotlist = plots,ncol = 4,nrow = 2)
p

ggsave("./peper/spe_OTU_cor/cor_spe_d.pdf", p, width = 12, height = 6)

