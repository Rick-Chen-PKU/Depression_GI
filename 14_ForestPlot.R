#######################################################################
#             Depression and gastrointestinal disorders               #
#######################################################################
# > sessionInfo()
# R version 3.6.1 (2019-07-05)
# platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Ubuntu 20.04.4 LTS
# Last modified 13 Sep 2022

# This script create forestplot for two-sample MR results 


#---------------------------------------------------------------------#
#                            Housekeeping                             #----
#---------------------------------------------------------------------#

# Clear environment
rm(list=ls()) #Remove any existing objects in R 

# Set working directory 
setwd("D:/WorkDir/Rstudio/Forestplot")
library(forestplot)
library(dplyr)

if (!require("pacman")) install.packages("pacman")
pacman::p_load("forestplot",  "dplyr")


#---------------------------------------------------------------------#
#                          Forestplot                         	      #----
#---------------------------------------------------------------------#
setwd("D:/WorkDir/Rstudio/Forestplot")
rs_forest_main <- read.csv('forest_main_copy.csv', header = FALSE)
rs_forest_replication <- read.csv('forest_replication_copy.csv', header = FALSE)
rs_forest_combind <- read.csv('forest_combined_copy.csv', header = FALSE)


# 读入数据的时候一定要把header设置成FALSE，确保第一行不被当作列名称。
# pdf(file = 'forest2.pdf',width = 7,height = 6, onefile = FALSE)  # 设定onefile = FALSE可以避免多出一页空白图

rs_forest = rs_forest_main
pdf(file = paste0('my_forest_plot', 1,'.pdf'),width = 7,height = 6, onefile = FALSE)
subgps <- setdiff(seq(29),c(1,2,9,16,23))
rs_forest$V1[subgps] <- paste("  ",rs_forest$V1[subgps]) 
forestplot(
  labeltext = as.matrix(rs_forest[,c(1,5,6)]), #设置用于文本展示的列
  graph.pos = 2, #设置森林图的位置，此处设置为2，则出现在第2列
  mean = rs_forest$V2, #设置均值
  lower = rs_forest$V3, #设置均值的lowlimits限
  upper = rs_forest$V4, #设置均值的uplimits限
  is.summary=c(T,T,F,F,F,F,F,F,T,F,F,F,F,F,F,T,F,F,F,F,F,F,T,F,F,F,F,F,F),
  #该参数接受一个逻辑向量，用于定义数据中每一行是否是汇总值，若是，则在对应位置设置为TRUE，若否，则设置为FALSE；设置为TRUE的行则以粗体出现
  zero = 1, #设置参照值，此处我们展示的是HR值，故参照值是1，而不是0
  boxsize = 0.12, #设置点估计的方形大小
  lineheight = unit(5,'mm'), #设置图形中的行距
  colgap = unit(3,'mm'), #设置图形中的列间距
  lwd.zero = 1.3, #设置参考线的粗细
  lwd.ci = 1, #设置区间估计线的粗细
  col=fpColors(box="royalblue",line="darkblue", summary="royalblue", hrz_lines = "#444444"),
  # col=fpColors(box=c('black'),
  #              summary=c('black'),
  #              lines = c('black'),
  #              zero = c('red')),
  #legend = c("Total effect", "Direct effect", "Indirect effect"),
  #legend_args = fpLegend(pos = list("topright"),title = "Effects",r = unit(.1, "snpc"),gp = gpar(col = "#CCCCCC", lwd = 1.5)),
  #使用fpColors()函数定义图形元素的颜色，从左至右分别对应点估计方形，汇总值，区间估计线，参考线
  xlab="Odds ratio",#设置x轴标签
  lwd.xaxis=1.2,#设置X轴线的粗细
  # clip = c(0.50,2),#x轴范围
  xticks = c(0.70, 1.0, 2),#x轴刻度
  xlog = TRUE,  # 对数标尺
  align=c("l","l"),#字符对齐方式，c=居中，l=左对齐，r=右对齐
  graphwidth = unit(35,"mm"), #森林图宽度
  lty.ci = "solid",
  txt_gp=fpTxtGp(label=gpar(cex=0.85), 
                 ticks=gpar(cex=0.7), 
                 xlab=gpar(cex = 0.8)),#文本大小
  ci.vertices = T#森林图可信区间两端添加小竖线
  )
graphics.off()  



rs_forest = rs_forest_replication
pdf(file = paste0('my_forest_plot', 2,'.pdf'),width = 7,height = 6, onefile = FALSE)
subgps <- setdiff(seq(29),c(1,2,9,16,23))
rs_forest$V1[subgps] <- paste("  ",rs_forest$V1[subgps]) 
forestplot(
  labeltext = as.matrix(rs_forest[,c(5,6)]), #设置用于文本展示的列
  graph.pos = 1, #设置森林图的位置，此处设置为2，则出现在第2列
  mean = rs_forest$V2, #设置均值
  lower = rs_forest$V3, #设置均值的lowlimits限
  upper = rs_forest$V4, #设置均值的uplimits限
  is.summary=c(T,T,F,F,F,F,F,F,T,F,F,F,F,F,F,T,F,F,F,F,F,F,T,F,F,F,F,F,F),
  #该参数接受一个逻辑向量，用于定义数据中每一行是否是汇总值，若是，则在对应位置设置为TRUE，若否，则设置为FALSE；设置为TRUE的行则以粗体出现
  zero = 1, #设置参照值，此处我们展示的是HR值，故参照值是1，而不是0
  boxsize = 0.12, #设置点估计的方形大小
  lineheight = unit(5,'mm'), #设置图形中的行距
  colgap = unit(3,'mm'), #设置图形中的列间距
  lwd.zero = 1.3, #设置参考线的粗细
  lwd.ci = 1, #设置区间估计线的粗细
  col=fpColors(box="royalblue",line="darkblue", summary="royalblue", hrz_lines = "#444444"),
  # col=fpColors(box=c('black'),
  #              summary=c('black'),
  #              lines = c('black'),
  #              zero = c('red')),
  #legend = c("Total effect", "Direct effect", "Indirect effect"),
  #legend_args = fpLegend(pos = list("topright"),title = "Effects",r = unit(.1, "snpc"),gp = gpar(col = "#CCCCCC", lwd = 1.5)),
  #使用fpColors()函数定义图形元素的颜色，从左至右分别对应点估计方形，汇总值，区间估计线，参考线
  xlab="Odds ratio",#设置x轴标签
  lwd.xaxis=1.2,#设置X轴线的粗细
  # clip = c(0.50,2),#x轴范围
  xticks = c(0.70, 1.0, 2),#x轴刻度
  xlog = TRUE,  # 对数标尺
  align=c("l","l"),#字符对齐方式，c=居中，l=左对齐，r=右对齐
  graphwidth = unit(35,"mm"), #森林图宽度
  lty.ci = "solid",
  txt_gp=fpTxtGp(label=gpar(cex=0.85), 
                 ticks=gpar(cex=0.7), 
                 xlab=gpar(cex = 0.8)),#文本大小
  ci.vertices = T#森林图可信区间两端添加小竖线
)
graphics.off()  


rs_forest = rs_forest_combind
pdf(file = paste0('my_forest_plot', 3,'.pdf'),width = 7,height = 6, onefile = FALSE)
subgps <- setdiff(seq(29),c(1,2,9,16,23))
rs_forest$V1[subgps] <- paste("  ",rs_forest$V1[subgps]) 
forestplot(
  labeltext = as.matrix(rs_forest[,c(5,6)]), #设置用于文本展示的列
  graph.pos = 1, #设置森林图的位置，此处设置为2，则出现在第2列
  mean = rs_forest$V2, #设置均值
  lower = rs_forest$V3, #设置均值的lowlimits限
  upper = rs_forest$V4, #设置均值的uplimits限
  is.summary=c(T,T,F,F,F,F,F,F,T,F,F,F,F,F,F,T,F,F,F,F,F,F,T,F,F,F,F,F,F),
  #该参数接受一个逻辑向量，用于定义数据中每一行是否是汇总值，若是，则在对应位置设置为TRUE，若否，则设置为FALSE；设置为TRUE的行则以粗体出现
  zero = 1, #设置参照值，此处我们展示的是HR值，故参照值是1，而不是0
  boxsize = 0.12, #设置点估计的方形大小
  lineheight = unit(5,'mm'), #设置图形中的行距
  colgap = unit(3,'mm'), #设置图形中的列间距
  lwd.zero = 1.3, #设置参考线的粗细
  lwd.ci = 1, #设置区间估计线的粗细
  col=fpColors(box="royalblue",line="darkblue", summary="royalblue", hrz_lines = "#444444"),
  # col=fpColors(box=c('black'),
  #              summary=c('black'),
  #              lines = c('black'),
  #              zero = c('red')),
  #legend = c("Total effect", "Direct effect", "Indirect effect"),
  #legend_args = fpLegend(pos = list("topright"),title = "Effects",r = unit(.1, "snpc"),gp = gpar(col = "#CCCCCC", lwd = 1.5)),
  #使用fpColors()函数定义图形元素的颜色，从左至右分别对应点估计方形，汇总值，区间估计线，参考线
  xlab="Odds ratio",#设置x轴标签
  lwd.xaxis=1.2,#设置X轴线的粗细
  # clip = c(0.50,2),#x轴范围
  xticks = c(0.70, 1.0, 2),#x轴刻度
  xlog = TRUE,  # 对数标尺
  align=c("l","l"),#字符对齐方式，c=居中，l=左对齐，r=右对齐
  graphwidth = unit(35,"mm"), #森林图宽度
  lty.ci = "solid",
  txt_gp=fpTxtGp(label=gpar(cex=0.85), 
                 ticks=gpar(cex=0.7), 
                 xlab=gpar(cex = 0.8)),#文本大小
  ci.vertices = T#森林图可信区间两端添加小竖线
)
graphics.off()







