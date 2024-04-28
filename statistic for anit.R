##Yun Zhu 12/04/2024

rm(list = ls())

if(T){
  library(dplyr)
  library(readxl)
  library(tidyr)
  library(tidyverse)
  library(reshape2)
  library(ggprism)
  library(rstatix)
  library(lawstat)
  library(ggpubr)
}


##提取数据
data <- read_excel("data/ANIT生化结果.xlsx", col_types = "guess") %>% 
  select(c(1, 3:7)) %>% 
  filter(!DATE == "24h After ANIT") %>% 
  select(-c("DATE","Groups"))
colnames(data) <- c("ID","group","variable","value")

########################使用函数判断是否正态然后统计分析
##############
#########Creat the statistic function
tf <- function(id){
  df %>%
    filter(variable == id) %>% 
    t_test(value ~ group,var.equal = TRUE)
}#用于正态分布，方差整齐

welch_tf <- function(id){
  df %>%
    filter(variable == id) %>% 
    t_test(value ~ group,var.equal = FALSE)
}#用于正态分布，方差不齐/单组非正态分布（单组n>30），方差整齐

wf <- function(id){
  df %>%
    filter(variable == id) %>% 
    wilcox_test(value ~ group)
}#用于两组全非正态分布/单组非正态分布（单组n<30）



######creat t-test function for pvalue
tf_p <- function( id ){
  #提取p值，按顺序提取行序号
  tt <- as.data.frame(tf(id))[,c("group1","group2","p")] 
  #增加以输入id名称命名一列
  tt$variable <- id
  print(tt)
}

weltf_p <- function( id ){
  #提取p值，按顺序提取行序号
  tt <- as.data.frame(welch_tf(id))[,c("group1","group2","p")] 
  #增加以输入id名称命名一列
  tt$variable <- id
  print(tt)
}

######creat wilcox_test function for pvalue
wf_p <- function( id ){
  tt <- as.data.frame(wf(id))[,c("group1","group2","p")]
  tt$variable <- id
  print(tt)
}


for (o in 2:length(unique(data$group))) {
  pair <- c("APC", unique(data$group)[o])
  print(pair)
  
rm(df,df_pvalue)
#######Creat empty list
ls1 <- list()
df <- filter(data,group %in% pair)
##循环判断符合正态分布且方差齐,分别算出p值后保存list
for (i in 1:length(unique(data$variable))) {
  tryCatch({
    print(i)
    ## 计算
    j <- unique(data$variable)[i]
    # 提取检测的各组数值
    test_data1 <- df %>%
      filter(group == unique(df$group)[1]) %>%
      filter(variable == j) %>%
      select(value) %>%
      unlist() %>%
      as.numeric()
    test_data2 <- df %>%
      filter(group == unique(df$group)[2]) %>%
      filter(variable == j) %>%
      select(value) %>%
      unlist() %>%
      as.numeric()
    ## 正态分布检测
    Ntest1 <- shapiro.test(test_data1)$p.value > 0.05
    Ntest2 <- shapiro.test(test_data2)$p.value > 0.05
    
    inputv <- j
    print(paste0(inputv, " ", unique(data$group)[1], " Normalization: ", Ntest1))
    print(paste0(inputv, " ", unique(data$group)[2], " Normalization: ", Ntest2))
    
    if (sum(Ntest1, Ntest2) > 0) {
      if (sum(Ntest1, Ntest2) == 2) {
        ## 正态分布方差齐性检测
        Vtest1 <- df %>%
          filter(variable == j) %>%
          bartlett.test(data=., value ~ group) %>%
          .$p.value > 0.05
        print(paste0(inputv, " Variance is equal: ", Vtest1))
        if (Vtest1) {
          print("t test")
          ls1[[inputv]] <- tf_p(inputv)
        } else {
          print("welch t test")
          ls1[[inputv]] <- weltf_p(inputv)
        }
      }
      if (sum(Ntest1, Ntest2) == 1) {
        if (length(test_data1) >= 30) {
          print(paste0("Sample Size: ", length(test_data1)))
          ## 非正态方差齐性检测
          Vtest2 <- df %>%
            filter(variable == j) %>%
            levene_test(data=., value ~ group) %>%
            .$p > 0.05
          print(paste0(inputv, " Variance is equal: ", Vtest2))
          if (Vtest2) {
            print("welch t test")
            ls1[[inputv]] <- weltf_p(inputv)
          } else {
            print("wilcoxon")
            ls1[[inputv]] <- wf_p(inputv)
          }
        } else {
          print(paste0("Sample Size: ", length(test_data1)))
          print("wilcoxon")
          ls1[[inputv]] <- wf_p(inputv)
        }
      }
    } else {
      print("wilcoxon")
      ls1[[inputv]] <- wf_p(inputv)
    }
  }, error = function(e) {
    # 当发生错误时，打印错误消息并继续下一个迭代 
    cat("Error occurred:", conditionMessage(e), "\n")
  })
 }
}
######将保存的list转为dataframe,删除行名
df_pvalue <- do.call(rbind, ls1)  
rownames(df_pvalue) <- NULL


stat_test <- df_pvalue
##增加p label
stat_test <- add_significance(stat_test,"p")


ls <- list()##建立空白list

for (u in as.character(unique(df$variable))) {
  ##创建根据对比组设置的数据框
  tes <-  t(pair) %>% 
    as.data.frame()
  rownames(tes) <- NULL
  ##新建两列和stat.test匹配，variable同样匹配
  colnames(tes) <- c("group1","group2")
  tes$variable <- u
  ##选择一个磷脂的所有数值
  tt <- df %>% filter(variable == u ) 
  ##根据磷脂的最大值设置放置label的Y轴位置
  tes$y.position <- NA
  tes$y.position[1] <-max(tt$value)*1.12
  ls[[u]] <- tes##保存为list
}

##list转为dataframe
stat.yp <- do.call(rbind, ls) 
rownames(stat.yp) <- NULL

##合并统计结果和位置dataframe
stat.test2 <- left_join(stat_test,stat.yp,by=c("variable","group1","group2"))

df <- filter(df,variable %in% stat.test2$variable)
##根据原始值画出柱状图
p = ggbarplot(df, x = "group", y = "value", add = c("mean_se"),color = "group",size=1.5,
              palette = c("#FF0000","#339900"),facet.by = "variable")+
  theme_bw() +
  facet_wrap(~ variable,scales = "free")+
  labs(x = NULL, y = NULL)
p
##添加上p label
p1= p + stat_pvalue_manual(stat.test2,label = "p.signif")
p1
##美化处理
p2 = p1+geom_jitter(aes(shape = factor(group),color=factor(group)), 
                    width = 0.2, size = 2,height =0)+
  scale_shape_prism()+
  theme_prism(base_size = 14) +
  theme(legend.position = "none")
p2

ggsave(p2,file = paste0("output/compair with anit/",substr(pair[2], 1, 7),".png"),
       height = 15,width = 18) 







