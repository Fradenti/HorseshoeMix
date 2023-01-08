library(tidyverse)
source('code/R/Auxiliary_functions/DataGenerator.R')

# -------------------------------------------------------------------------


strata_mean <- function(X){
  x1 = cbind(
    rowMeans(X[,1:100]),
    rowMeans(X[,101:200]),
    rowMeans(X[,-c(1:200)]))
  
  return(x1)  
  
}   



# A, B, C
# nsim - burnin - thinning - n - p - nb
s1 <- c(5000, 5000, 5, 500, 250, 100)
s2 <- c(5000, 5000, 5, 500, 500, 100)
s3 <- c(5000, 5000, 5, 500, 750, 100)
Scenario <- list(s1,s2,s3)
seedlist <-round( seq(0,5900,length.out=30) )

Y1 <- Data.generator.regression(n = Scenario[[1]][4],
                                p = Scenario[[1]][5],
                                nb1 = Scenario[[1]][6],
                                seedlist = seedlist,
                                NSIM = 30,
                                sdX = 1)
Y2 <- Data.generator.regression(n = Scenario[[2]][4],
                                p = Scenario[[2]][5],
                                nb1 = Scenario[[2]][6],
                                seedlist = seedlist,
                                NSIM = 30,
                                sdX = 1)
Y3 <- Data.generator.regression(n = Scenario[[3]][4],
                                p = Scenario[[3]][5],
                                nb1 = Scenario[[3]][6],
                                seedlist = seedlist,
                                NSIM = 30,
                                sdX = 1)


# -------------------------------------------------------------------------


hsm1 <- readRDS("new_simu_HSmix_1.RDS")
hsm2 <- readRDS("new_simu_HSmix_2.RDS")
hsm3 <- readRDS("new_simu_HSmix_3.RDS")
length(hsm1)

# first scenario
s1_TRUTH_B <- do.call(rbind,map(Y1, ~.x$betat))
s1_HSM_B <- do.call(rbind,map(hsm1, ~.x$beta))
# second scenario
s2_TRUTH_B <- do.call(rbind,map(Y2, ~.x$betat))
s2_HSM_B <- do.call(rbind,map(hsm2, ~.x$beta))
# third scenario
s3_TRUTH_B <- do.call(rbind,map(Y3, ~.x$betat))
s3_HSM_B <- do.call(rbind,map(hsm3, ~.x$beta))

SE_hsm1 <- (s1_TRUTH_B - s1_HSM_B)^2
SE_hsm2 <- (s2_TRUTH_B - s2_HSM_B)^2
SE_hsm3 <- (s3_TRUTH_B - s3_HSM_B)^2


mean( SE_hsm1 )
mean( SE_hsm2 )
mean( SE_hsm3 )

sd(rowMeans(SE_hsm1)) 
sd(rowMeans(SE_hsm2))
sd(rowMeans(SE_hsm3))

   

SC1_B123_HSM <- data.frame(reshape2::melt(strata_mean(SE_hsm1)), sc = 1, type= "HSM")
SC2_B123_HSM <- data.frame(reshape2::melt(strata_mean(SE_hsm2)), sc = 2, type= "HSM")
SC3_B123_HSM <- data.frame(reshape2::melt(strata_mean(SE_hsm3)), sc = 3, type= "HSM")

HSM_123 <- rbind(SC1_B123_HSM,SC2_B123_HSM,SC3_B123_HSM)
# saveRDS(HSM_123,"new_syntesis_HSM.RDS")
#######



brhs <- readRDS("Bayereg_HS.RDS")
brhs_p <- readRDS("Bayereg_HSplus.RDS")
brlasso <- readRDS("Bayereg_Lasso.RDS")


length(brhs[[1]])
brhs[[1]][[1]]$mu.beta
map(brhs[[1]],~.x$mu.beta)
map(brhs[[2]],~.x$mu.beta)
map(brhs[[3]],~.x$mu.beta)



s1_HS_B  <- t(do.call(cbind,map(brhs[[1]], ~.x$mu.beta)))
s2_HS_B  <- t(do.call(cbind,map(brhs[[2]], ~.x$mu.beta)))
s3_HS_B  <- t(do.call(cbind,map(brhs[[3]], ~.x$mu.beta)))
s1_HSPlus_B <- t(do.call(cbind,map(brhs_p[[1]], ~.x$mu.beta)))
s2_HSPlus_B <- t(do.call(cbind,map(brhs_p[[2]], ~.x$mu.beta)))
s3_HSPlus_B <- t(do.call(cbind,map(brhs_p[[3]], ~.x$mu.beta)))
s1_LS_B  <- t(do.call(cbind,map(brlasso[[1]], ~.x$mu.beta)))
s2_LS_B  <- t(do.call(cbind,map(brlasso[[2]], ~.x$mu.beta)))
s3_LS_B  <- t(do.call(cbind,map(brlasso[[3]], ~.x$mu.beta)))



SE_HS1 <- (s1_TRUTH_B - s1_HS_B)^2
SE_HS2 <- (s2_TRUTH_B - s2_HS_B)^2
SE_HS3 <- (s3_TRUTH_B - s3_HS_B)^2

SE_HSPlus1 <- (s1_TRUTH_B - s1_HSPlus_B)^2
SE_HSPlus2 <- (s2_TRUTH_B - s2_HSPlus_B)^2
SE_HSPlus3 <- (s3_TRUTH_B - s3_HSPlus_B)^2

SE_LS1 <- (s1_TRUTH_B - s1_LS_B)^2
SE_LS2 <- (s2_TRUTH_B - s2_LS_B)^2
SE_LS3 <- (s3_TRUTH_B - s3_LS_B)^2



mean( SE_HS1 )
mean( SE_HS2 )
mean( SE_HS3 )
mean( SE_HSPlus1 )
mean( SE_HSPlus2 )
mean( SE_HSPlus3 )
mean( SE_LS1 )
mean( SE_LS2 )
mean( SE_LS3 )

sd(rowMeans(SE_HS1)) 
sd(rowMeans(SE_HS2))
sd(rowMeans(SE_HS3))
sd(rowMeans(SE_HSPlus1)) 
sd(rowMeans(SE_HSPlus2))
sd(rowMeans(SE_HSPlus3))
sd(rowMeans(SE_LS1)) 
sd(rowMeans(SE_LS2))
sd(rowMeans(SE_LS3))




SC1_B123_HS <- data.frame(reshape2::melt(strata_mean(SE_HS1)), sc = 1, type = "HS")
SC2_B123_HS <- data.frame(reshape2::melt(strata_mean(SE_HS2)), sc = 2, type = "HS")
SC3_B123_HS <- data.frame(reshape2::melt(strata_mean(SE_HS3)), sc = 3, type = "HS")
SC1_B123_HSPlus <- data.frame(reshape2::melt(strata_mean(SE_HSPlus1)), sc = 1, type = "HS+")
SC2_B123_HSPlus <- data.frame(reshape2::melt(strata_mean(SE_HSPlus2)), sc = 2, type = "HS+")
SC3_B123_HSPlus <- data.frame(reshape2::melt(strata_mean(SE_HSPlus3)), sc = 3, type = "HS+")
SC1_B123_LS <- data.frame(reshape2::melt(strata_mean(SE_LS1)), sc = 1, type="Lasso")
SC2_B123_LS <- data.frame(reshape2::melt(strata_mean(SE_LS2)), sc = 2, type="Lasso")
SC3_B123_LS <- data.frame(reshape2::melt(strata_mean(SE_LS3)), sc = 3, type="Lasso")




HS_123 <- rbind(SC1_B123_HS,SC2_B123_HS,SC3_B123_HS)
saveRDS(HS_123,"syntesis_HS.RDS")
HSPlus_123 <- rbind(SC1_B123_HSPlus,SC2_B123_HSPlus,SC3_B123_HSPlus)
saveRDS(HSPlus_123,"syntesis_HSPlus.RDS")
LS_123 <- rbind(SC1_B123_LS,SC2_B123_LS,SC3_B123_LS)
saveRDS(LS_123,"syntesis_LS.RDS")




# -------------------------------------------------------------------------
library(tidyverse)
HSM_123    <- readRDS("new_syntesis_HSM.RDS")
HS_123     <- readRDS("syntesis_HS.RDS")
HSPlus_123 <- readRDS("syntesis_HSPlus.RDS")
LS_123     <- readRDS("syntesis_LS.RDS")


library(patchwork)

D <- rbind(HSM_123,HS_123,HSPlus_123,LS_123)
D2 = D %>% mutate(type = factor(type,levels = c("HSM","HS","HS+","Lasso")),
                  sc2 = paste("Scenario",sc)) 

theme_set(theme_bw())

p1 <- ggplot(D2 %>% filter(sc==1))+
  geom_boxplot(aes(x=factor(Var2),y=value,col=type,))+
  facet_wrap(~sc2, nrow = 2) +
  scale_color_brewer("Model", palette = "Set1",
                    guide = guide_legend(
                      direction = "horizontal",
                      title.position = "top"
                    ))+
  theme(text= element_text(size=15),legend.direction="horizontal")+xlab("Block")+
  ylab("Stratified MSE")

p2 <- ggplot(D2 %>% filter(sc==2 & type != "Lasso"))+
  geom_boxplot(aes(x=factor(Var2),y=value,col=type))+
  facet_wrap(~sc2, nrow = 2) +
  scale_color_brewer("Model", palette = "Set1")+
  theme(text= element_text(size=15),legend.position = "none")+xlab("Block")+
  ylab("Stratified MSE")

p3 <- ggplot(D2 %>% filter(sc==3 & type != "Lasso"))+
  geom_boxplot(aes(x=factor(Var2),y=value,col=type))+
  facet_wrap(~sc2, nrow = 2) +
  scale_color_brewer("Model", palette = "Set1")+
  theme(text= element_text(size=15),legend.position = "none")+xlab("Block")+
  ylab("Stratified MSE")


(p1+p2+p3+guide_area())+plot_layout(guides="collect")


# BW option

GC <- c("#FFFFFF",
        "#EBEBEB",
        "#C8C8C8",
        "#525252")


p1 <- ggplot(D2 %>% filter(sc==1))+
  geom_boxplot(aes(x=factor(Var2),y=value,fill=type))+
  facet_wrap(~sc2, nrow = 2) +
  scale_fill_manual("Model",values = GC,
                  guide = guide_legend(
                          direction = "horizontal",
                          title.position = "top",ncol = 4
                        ))+
  theme(text= element_text(size=18),legend.direction="horizontal",
        legend.title = element_text(size=18),
        legend.text = element_text(size=18))+
  xlab("Block")+
  ylab("Stratified MSE")
p1

p2 <- ggplot(D2 %>% filter(sc==2 & type != "Lasso"))+
  geom_boxplot(aes(x=factor(Var2),y=value,fill=type))+
  facet_wrap(~sc2, nrow = 2) +
  scale_fill_manual("Model",values = GC[1:3])+
  theme(text= element_text(size=18),legend.position = "none")+xlab("Block")+
  ylab("Stratified MSE")

p3 <- ggplot(D2 %>% filter(sc==3 & type != "Lasso"))+
  geom_boxplot(aes(x=factor(Var2),y=value,fill=type))+
  facet_wrap(~sc2, nrow = 2) +
  scale_fill_manual("Model",values = GC[1:3])+
  theme(text= element_text(size=18),legend.position = "none")+xlab("Block")+
  ylab("Stratified MSE")


(p1+p2+p3+guide_area())+plot_layout(guides="collect")



labeldat1 = D2 %>% filter(sc==1 ) %>%
  group_by(Var2 = factor(Var2), type) %>%
  summarize(ypos = min(value) - .0002 ) 

labeldat2 = D2 %>% filter(sc==2 & type != "Lasso") %>%
  group_by(Var2 = factor(Var2), type) %>%
  summarize(ypos = min(value) - .0004 ) 

labeldat3 = D2 %>% filter(sc==3 & type != "Lasso") %>%
  group_by(Var2 = factor(Var2), type) %>%
  summarize(ypos = min(value) - .0004 ) 


p1 <- ggplot()+
  geom_boxplot(data = D2 %>% filter(sc==1), 
               mapping = aes(x=factor(Var2),y=value,fill=type),
               position = position_dodge(width = .8))+
  facet_wrap(~sc2, nrow = 2) +
  geom_text(data=labeldat1,mapping = aes(x=factor(Var2),y=ypos,group = type,label=type),
             position = position_dodge(width = .8))+
  theme(text= element_text(size=15),
        legend.direction="horizontal",
        legend.position = "none")+xlab("Block")+
  scale_fill_manual("Model",values = GC)+
  ylab("Stratified MSE")
p1
p2 <- ggplot(D2 %>% filter(sc==2 & type != "Lasso"))+
  geom_boxplot(aes(x=factor(Var2),y=value,fill=type))+
  facet_wrap(~sc2, nrow = 2) +
  geom_text(data=labeldat2,mapping = aes(x=factor(Var2),y=ypos,group = type,label=type),
            position = position_dodge(width = .8))+
  scale_fill_manual("Model",values = GC[1:3])+
  theme(text= element_text(size=15),legend.position = "none")+xlab("Block")+
  ylab("Stratified MSE")

p3 <- ggplot(D2 %>% filter(sc==3 & type != "Lasso"))+
  geom_boxplot(aes(x=factor(Var2),y=value,fill=type))+
  facet_wrap(~sc2, nrow = 2) +
  geom_text(data=labeldat3,mapping = aes(x=factor(Var2),y=ypos,group = type,label=type),
            position = position_dodge(width = .8))+
  scale_fill_manual("Model",values = GC[1:3])+
  theme(text= element_text(size=15),legend.position = "none")+xlab("Block")+
  ylab("Stratified MSE")


(p1/(p2+p3))
