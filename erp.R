require(lme4)
require(ggplot2)
#require(plyr)
require(dplyr)
require(ggsignif)
require(ggpubr)
require(lmerTest)
require(emmeans)
require(rsq)
require(tidyverse)
require(broom)
require(broom.mixed)
require(data.table)
require(tibble)
require(memisc)




data <- read.csv('/Users/agnesgao/Downloads/ag_data_main_effects_all_subs_all_words_data_merged (1).csv', header=T)

data$naming<-NA
data$naming[which(data$resp_relate == 1)] <- 'Same Related'
data$naming[which(data$resp_relate == 2)] <- 'Different Related'
data$naming[which(data$resp_relate == 3)] <- 'Different Rerelated'
data$naming[which(data$resp_relate == 4)] <- 'No Response'
data$naming[which(data$resp_relate == 5)] <- 'Nonwords'
data$naming<-as.factor(data$naming)

#group length: 1-4 short; >5 long
mid <- c("MID")
replacement_values <- c("LONG")

#Use replace() to replace the mids
data$leng2 <- replace(data$leng2, data$leng2 %in% mid, replacement_values)


#get average between 300-500ms
data$n400 <- rowMeans(subset(data[,c(132:182)]),na.rm = T)

#get average between 200-300ms
data$n250 <- rowMeans(subset(data[,c(121:146)]),na.rm = T)

#get average 80 - 150ms
data$p1 <- rowMeans(subset(data[,c(71:88)]),na.rm = T) #80 - 150

#get average  160ms - 240ms
data$p2 <- rowMeans(subset(data[,c(91:111)]),na.rm = T)

#rank re-order for naming and anteriority
data$naming <- factor(data$naming, levels=c("Same Related","Different Related","Different Rerelated"))
data$Anteriority <-factor(data$Anteriority, levels=c("Frontal","Central","Posterior"))

#N400 
#midline
model1<-lmer(n400 ~ naming*Anteriority + (1|PartID) + (1|electrode) + (1|wordID), 
             data=subset(data, subset = naming %in% c("Same Related","Different Related","Different Rerelated") & Hemisphere == "Z"))
summary(model1)

mLME1 <- emmeans(model1, pairwise~naming|Anteriority, mode = "satterthwaite", lmerTest.limit = 240000,adjust='tukey')
summary(mLME1)$contrasts
#lat
model1_1<-lmer(n400 ~ naming*Anteriority + (1|PartID) + (1|electrode) + (1|wordID), 
             data=subset(data, subset = naming %in% c("Same Related","Different Related","Different Rerelated") & Hemisphere != "Z"))
summary(model1_1)

mLME1_1 <- emmeans(model1_1, pairwise~naming|Anteriority, mode = "satterthwaite", lmerTest.limit = 255096,adjust='tukey')
summary(mLME1_1)$contrasts


eff_size(mLME1_1, sigma = sigma(model1_1), edf = df.residual(model1_1))



#n250
#midline
model1_2<-lmer(n250 ~ naming*Anteriority + (1|PartID) + (1|electrode) + (1|wordID), 
             data=subset(data, subset = naming %in% c("Same Related","Different Related","Different Rerelated") & Hemisphere == "Z"))
summary(model1_2)

write.csv( summary(model1_2)$coefficients,"/Users/agnesgao/Documents/pr/n250_mid_main.csv" )

mLME1_2 <- emmeans(model1_2, pairwise~naming|Anteriority, mode = "satterthwaite", lmerTest.limit = 240000,adjust='tukey')
summary(mLME1_2)$contrasts


#lat
model1_3<-lmer(n250 ~ naming*Anteriority + (1|PartID) + (1|electrode) + (1|wordID), 
               data=subset(data, subset = naming %in% c("Same Related","Different Related","Different Rerelated") & Hemisphere != "Z"))
summary(model1_3)

write.csv( summary(model1_3)$coefficients,"/Users/agnesgao/Documents/pr/n250_lat_main.csv" )

mLME1_3 <- emmeans(model1_3, pairwise~naming|Anteriority, mode = "satterthwaite", lmerTest.limit = 255096,adjust='tukey')
summary(mLME1_3)$contrasts


#p2
#midline
model1_4<-lmer(p2 ~ naming*Anteriority + (1|PartID) + (1|electrode) + (1|wordID), 
               data=subset(data, subset = naming %in% c("Same Related","Different Related","Different Rerelated") & Hemisphere == "Z"))
summary(model1_4)

write.csv( summary(model1_4)$coefficients,"/Users/agnesgao/Documents/pr/p2_mid_main.csv" )

mLME1_4 <- emmeans(model1_4, pairwise~naming|Anteriority, mode = "satterthwaite", lmerTest.limit = 240000,adjust='tukey')
summary(mLME1_4)$contrasts


#lat
model1_5<-lmer(p2 ~ naming*Anteriority + (1|PartID) + (1|electrode) + (1|wordID), 
               data=subset(data, subset = naming %in% c("Same Related","Different Related","Different Rerelated") & Hemisphere != "Z"))
summary(model1_5)

write.csv( summary(model1_5)$coefficients,"/Users/agnesgao/Documents/pr/p2_lat_main.csv" )

mLME1_5 <- emmeans(model1_5, pairwise~naming, mode = "satterthwaite", lmerTest.limit = 255096,adjust='tukey')
summary(mLME1_5)$contrasts


#p1
#midline
model1_7<-lmer(p1 ~ naming*Anteriority + (1|PartID) + (1|electrode) + (1|wordID), 
               data=subset(data, subset = naming %in% c("Same Related","Different Related","Different Rerelated") & Hemisphere == "Z"))
summary(model1_7)

write.csv( summary(model1_7)$coefficients,"/Users/agnesgao/Documents/pr/p1_mid_main.csv" )

mLME1_7 <- emmeans(model1_7, pairwise~naming, mode = "satterthwaite", lmerTest.limit = 255096,adjust='tukey')
summary(mLME1_7)$contrasts

#lat
model1_6<-lmer(p1 ~ naming*Anteriority + (1|PartID) + (1|electrode) + (1|wordID), 
               data=subset(data, subset = naming %in% c("Same Related","Different Related","Different Rerelated") & Hemisphere != "Z"))
summary(model1_6)

write.csv( summary(model1_6)$coefficients,"/Users/agnesgao/Documents/pr/p1_lat_main.csv" )

mLME1_6 <- emmeans(model1_6, pairwise~naming, mode = "satterthwaite", lmerTest.limit = 255096,adjust='tukey')
summary(mLME1_6)$contrasts

#p2.conc
#lat
model.p2.conc.lat<-lmer(p2 ~ naming*concg2*Anteriority + (1|PartID) + (1|electrode) + (1|wordID), 
               data=subset(data, subset = naming %in% c("Same Related","Different Related","Different Rerelated") & Hemisphere != "Z"))
summary(model.p2.conc.lat)

write.csv( summary(model.p1.conc.lat)$coefficients,"/Users/agnesgao/Documents/pr/p1_lat_conc.csv" )

mLME1_8 <- emmeans(model.p2.conc.lat, pairwise~concg2|naming, mode = "satterthwaite", lmerTest.limit = 255096,adjust='tukey')
summary(mLME1_8)$contrasts

#mid
model.p1.conc.mid<-lmer(p1 ~ naming*concg2*Anteriority + (1|PartID) + (1|electrode) + (1|wordID), 
                        data=subset(data, subset = naming %in% c("Same Related","Different Related","Different Rerelated") & Hemisphere != "Z"))
summary(model.p1.conc.mid)

write.csv( summary(model.p1.conc.mid)$coefficients,"/Users/agnesgao/Documents/pr/p1_mid_conc.csv" )

mLME1_9 <- emmeans(model.p1.conc.mid, pairwise~concg2|naming*Anteriority, mode = "satterthwaite", lmerTest.limit = 255096,adjust='tukey')
summary(mLME1_9)$contrasts


###p2 ortho
#lat
model.p2.ortho.lat<-lmer(p2 ~ naming*orthog2*Anteriority + (1|PartID) + (1|electrode) + (1|wordID), 
                        data=subset(data, subset = naming %in% c("Same Related","Different Related","Different Rerelated") & Hemisphere != "Z"))
summary(model.p2.ortho.lat)

write.csv( summary(model.p2.ortho.lat)$coefficients,"/Users/agnesgao/Documents/pr/p2_lat_ortho.csv" )

mLME1_10 <- emmeans(model.p2.ortho.lat, pairwise~orthog2|naming, mode = "satterthwaite", lmerTest.limit = 255096,adjust='tukey')
summary(mLME1_10)$contrasts

#mid
model.p2.ortho.mid<-lmer(p2 ~ naming*orthog2*Anteriority + (1|PartID) + (1|electrode) + (1|wordID), 
                        data=subset(data, subset = naming %in% c("Same Related","Different Related","Different Rerelated") & Hemisphere == "Z"))
summary(model.p2.ortho.mid)

write.csv( summary(model.p2.ortho.mid)$coefficients,"/Users/agnesgao/Documents/pr/p2_mid_ortho.csv" )

mLME1_11 <- emmeans(model.p2.ortho.mid, pairwise~orthog2|naming*Anteriority, mode = "satterthwaite", lmerTest.limit = 255096,adjust='tukey')
summary(mLME1_11)$contrasts



#conc
#MIDLINE
model.c1<-lmer(n400 ~ naming*concg2*Anteriority + (1|PartID) + (1|electrode) + (1|wordID), 
               data=subset(data, subset = naming %in% c("Same Related","Different Related","Different Rerelated") & Hemisphere == "Z"))
summary(model.c1)


write.csv( summary(model.c1)$coefficients,"/Users/agnesgao/Documents/pr/conc_n4_mid.csv" )

mLME.c1 <- emmeans(model.c1, pairwise~concg2|naming*Anteriority, mode = "satterthwaite", lmerTest.limit = 255096,adjust='tukey')
summary(mLME.c1)$contrasts
#LAT
model.c2<-lmer(n400 ~ naming*concg2*Anteriority + (1|PartID) + (1|Hemisphere), 
               data=subset(data, subset = naming %in% c("Same Related","Different Related","Different Rerelated") & Hemisphere != "Z"))
summary(model.c2)
write.csv( summary(model.c2)$coefficients,"/Users/agnesgao/Documents/pr/conc_n4_lat.csv" )

mLME.c2 <- emmeans(model.c2, pairwise~concg2|naming*Anteriority, mode = "satterthwaite", lmerTest.limit = 255096, adjust='tukey')
summary(mLME.c2)$contrasts

eff_size(mLME.c2, sigma = sigma(model.c2), edf = df.residual(model.c2))

#ortho*

model1_2<-lmer(n400 ~ naming*orthog2*Anteriority + (1|PartID) + (1|electrode) + (1|wordID), 
               data=subset(data, subset = naming %in% c("Same Related","Different Related","Different Rerelated") & Hemisphere == "Z"))
summary(model1_2)
mLME1_2 <- emmeans(model1_2, pairwise~orthog2|naming, mode = "satterthwaite", lmerTest.limit = 240000, adjust='sidak')
summary(mLME1_2)$contrasts


#length*

model1_3<-lmer(n400 ~ naming*leng2*Anteriority + (1|PartID) + (1|electrode) + (1|wordID), 
               data=subset(data, subset = naming %in% c("Same Related","Different Related","Different Rerelated") & Hemisphere != "Z"))
summary(model1_3)
mLME1_3 <- emmeans(model1_3, pairwise~leng2|naming*Anteriority, mode = "satterthwaite", lmerTest.limit = 240000, adjust='sidak')
summary(mLME1_3)$contrasts

#lat
model2<-lmer(n400 ~ naming*Anteriority + (1|PartID) + (1|electrode) + (1|wordID), 
             data=subset(data, subset = naming %in% c("Same Related","Different Related","Different Rerelated")  & Hemisphere != "Z"))
summary(model2)

mLME2 <- emmeans(model2, pairwise~naming, mode = "satterthwaite", lmerTest.limit = 240000, adjust='sidak')
summary(mLME2)$contrasts

#conc
model2_1<-lmer(n400 ~ naming*concg2*Anteriority + (1|PartID) + (1|electrode) + (1|wordID), 
               data=subset(data, subset = naming %in% c("Same Related","Different Related","Different Rerelated") & Hemisphere != "Z"))
summary(model2_1)
mLME2_1 <- emmeans(model2_1, pairwise~concg2|naming, mode = "satterthwaite", lmerTest.limit = 240000, adjust='sidak')
summary(mLME2_1)$contrasts

#ortho

model2.m<-lmer(n400 ~ naming*orthog2*Anteriority + (1|PartID) + (1|electrode) + (1|wordID), 
               data=subset(data, subset = naming %in% c("Same Related","Different Related","Different Rerelated") & Hemisphere == "Z"))
summary(model2.m)

write.csv( summary(model2.m)$coefficients,"/Users/agnesgao/Documents/pr/ortho_n4_mid.csv" )


mLME2.m <- emmeans(model2.m, pairwise~orthog2|naming, mode = "satterthwaite", lmerTest.limit = 240000, adjust='sidak')
summary(mLME2.m)$contrasts

model2_2<-lmer(n400 ~ naming*orthog2*Anteriority + (1|PartID) + (1|electrode) + (1|wordID), 
               data=subset(data, subset = naming %in% c("Same Related","Different Related","Different Rerelated") & Hemisphere != "Z"))
summary(model2_2)

write.csv( summary(model2_2)$coefficients,"/Users/agnesgao/Documents/pr/ortho_n4_lat.csv" )


mLME2_2 <- emmeans(model2_2, pairwise~orthog2|naming*Anteriority, mode = "satterthwaite", lmerTest.limit = 255096, adjust='tukey')
summary(mLME2_2)$contrasts


eff_size(mLME2_2, sigma = sigma(model2_2), edf = df.residual(model2_2))

#length

model2_3<-lmer(n400 ~ naming*leng2*Anteriority + (1|PartID) + (1|electrode) + (1|wordID), 
               data=subset(data, subset = naming %in% c("Same Related","Different Related","Different Rerelated") & Hemisphere != "Z"))
summary(model2_3)
mLME2_3 <- emmeans(model2_3, pairwise~leng2|naming, mode = "satterthwaite", lmerTest.limit = 240000, adjust='sidak')
summary(mLME2_3)$contrasts

#phono

model2_4<-lmer(n400 ~ naming*PHONO2_G*Anteriority + (1|PartID) + (1|electrode) + (1|wordID), 
               data=subset(data, subset = naming %in% c("Same Related","Different Related","Different Rerelated") & Hemisphere != "Z"))
summary(model2_4)
mLME2_4 <- emmeans(model2_4, pairwise~PHONO2_G|naming*Anteriority, mode = "satterthwaite", lmerTest.limit = 240000, adjust='sidak')
summary(mLME2_4)$contrasts

#N250
#midline
model3<-lmer(n250 ~ naming*Anteriority + (1|PartID) + (1|electrode) + (1|wordID), 
             data=subset(data, subset = naming %in% c("Same Related","Different Related","Different Rerelated") & Hemisphere == "Z"))
summary(model3)

mLME3 <- emmeans(model3, pairwise~naming|Anteriority, mode = "satterthwaite", lmerTest.limit = 240000, adjust='sidak')
summary(mLME3)$contrasts

#lat

model4<-lmer(n250 ~ naming*Anteriority + (1|PartID) + (1|electrode) + (1|wordID), 
             data=subset(data, subset = naming %in% c("Same Related","Different Related","Different Rerelated")  & Hemisphere != "Z"))
summary(model4)

mLME4 <- emmeans(model4, pairwise~naming|Anteriority, mode = "satterthwaite", lmerTest.limit = 240000, adjust='sidak')
summary(mLME4)$contrasts

#conc
model4_1<-lmer(n250 ~ naming*concg2*Anteriority + (1|PartID) + (1|electrode) + (1|wordID), 
               data=subset(data, subset = naming %in% c("Same Related","Different Related","Different Rerelated") & Hemisphere != "Z"))
summary(model4_1)
mLME4_1 <- emmeans(model4_1, pairwise~concg2|naming*Anteriority, mode = "satterthwaite", lmerTest.limit = 240000, adjust='sidak')
summary(mLME4_1)$contrasts

#ortho***

model4_2<-lmer(n250 ~ naming*orthog2*Anteriority + (1|PartID) + (1|electrode) + (1|wordID), 
               data=subset(data, subset = naming %in% c("Same Related","Different Related","Different Rerelated") & Hemisphere != "Z"))
summary(model4_2)

write.csv( summary(model4_2)$coefficients,"/Users/agnesgao/Documents/pr/orth_n250_lat.csv" )
mLME4_2 <- emmeans(model4_2, pairwise~orthog2|naming, mode = "satterthwaite", lmerTest.limit = 260000, adjust='tukey')
summary(mLME4_2)$contrasts

model4_3<-lmer(n250 ~ naming*orthog2*Anteriority + (1|PartID) + (1|electrode) + (1|wordID), 
               data=subset(data, subset = naming %in% c("Same Related","Different Related","Different Rerelated") & Hemisphere == "Z"))
summary(model4_3)
write.csv( summary(model4_3)$coefficients,"/Users/agnesgao/Documents/pr/orth_n250_mid.csv" )

#length

model4_3<-lmer(n400 ~ naming*leng2*Anteriority + (1|PartID) + (1|electrode) + (1|wordID), 
               data=subset(data, subset = naming %in% c("Same Related","Different Related","Different Rerelated") & Hemisphere != "Z"))
summary(model4_3)
write.csv( summary(model4_3)$coefficients,"/Users/agnesgao/Documents/pr/len_n4_lat.csv" )
mLME4_3 <- emmeans(model4_3, pairwise~leng2|naming*Anteriority, mode = "satterthwaite", lmerTest.limit = 260000, adjust='tukey')
summary(mLME4_3)$contrasts

model4_4<-lmer(n400 ~ naming*leng2*Anteriority + (1|PartID) + (1|electrode) + (1|wordID), 
               data=subset(data, subset = naming %in% c("Same Related","Different Related","Different Rerelated") & Hemisphere == "Z"))
summary(model4_4)
write.csv( summary(model4_4)$coefficients,"/Users/agnesgao/Documents/pr/len_n4_mid.csv" )
mLME4_4 <- emmeans(model4_4, pairwise~leng2|naming, mode = "satterthwaite", lmerTest.limit = 240000, adjust='sidak')
summary(mLME4_4)$contrasts

#phono - p2
#lat
model4_4<-lmer(p2 ~ naming*PHONO2_G*Anteriority + (1|PartID) + (1|electrode) + (1|wordID), 
               data=subset(data, subset = naming %in% c("Same Related","Different Related","Different Rerelated") & Hemisphere != "Z"))
summary(model4_4)

write.csv( summary(model4_4)$coefficients,"/Users/agnesgao/Documents/pr/pho_p2_lat.csv" )

mLME4_4 <- emmeans(model4_4, pairwise~PHONO2_G|naming*Anteriority, mode = "satterthwaite", lmerTest.limit = 260000, adjust='tukey')
summary(mLME4_4)$contrasts

#mid
model4_5<-lmer(p2 ~ naming*PHONO2_G*Anteriority + (1|PartID) + (1|electrode) + (1|wordID), 
               data=subset(data, subset = naming %in% c("Same Related","Different Related","Different Rerelated") & Hemisphere == "Z"))
summary(model4_5)

write.csv( summary(model4_5)$coefficients,"/Users/agnesgao/Documents/pr/pho_p2_mid.csv" )

mLME4_5 <- emmeans(model4_5, pairwise~PHONO2_G|naming*Anteriority, mode = "satterthwaite", lmerTest.limit = 260000, adjust='tukey')
summary(mLME4_5)$contrasts


########## P1 for length ################ PURELY VISUAL, NO PREDICTION EFFECT
#lat
model5<-lmer(p1 ~ naming*leng2*Anteriority + (1|PartID) + (1|electrode) + (1|wordID), 
               data=subset(data, subset = naming %in% c("Same Related","Different Related","Different Rerelated") & Hemisphere != "Z"))
summary(model5)
write.csv( summary(model5)$coefficients,"/Users/agnesgao/Documents/pr/len_la.csv" )
mLME5 <- emmeans(model5, pairwise~leng2|naming*Anteriority, mode = "satterthwaite", lmerTest.limit = 260000, adjust='tukey')
summary(mLME5)$contrasts
#midline
model5_1<-lmer(p1 ~ naming*leng2*Anteriority + (1|PartID) + (1|electrode) + (1|wordID), 
             data=subset(data, subset = naming %in% c("Same Related","Different Related","Different Rerelated") & Hemisphere == "Z"))
summary(model5_1)

summary(model5_1)$coefficients 
write.csv( summary(model5_1)$coefficients,"/Users/agnesgao/Documents/pr/len.csv" )

mLME5_1 <- emmeans(model5_1, pairwise~leng2|naming, mode = "satterthwaite", lmerTest.limit = 260000, adjust='tukey')
summary(mLME5_1)$contrasts

####

########### P2 ################
#lat
model6<-lmer(p2 ~ naming*leng2*Anteriority + (1|PartID) + (1|electrode) + (1|wordID), 
             data=subset(data, subset = naming %in% c("Same Related","Different Related","Different Rerelated") & Hemisphere != "Z"))
summary(model6)
mLME6 <- emmeans(model6, pairwise~leng2|naming, mode = "satterthwaite", lmerTest.limit = 240000, adjust='sidak')
summary(mLME6)$contrasts

#midline
model6_1<-lmer(p2 ~ naming*leng2*Anteriority + (1|PartID) + (1|electrode) + (1|wordID), 
               data=subset(data, subset = naming %in% c("Same Related","Different Related","Different Rerelated") & Hemisphere == "Z"))
summary(model6_1)
mLME6_1 <- emmeans(model6_1, pairwise~leng2|naming, mode = "satterthwaite", lmerTest.limit = 240000, adjust='sidak')
summary(mLME6_1)$contrasts


############################by item IGNORE #################################
word <- subset(data,subset = condition %in% c("Predicted","Unpredicted","Unrelated")& Hemisphere != "Z")%>% 
  group_by(wordID,condition,electrode,Hemisphere, Anteriority,PartID,target) %>% summarise(
    w_n400 = mean(n400,na.rm = T),
    w_n250 = mean(n250,na.rm = T))

word$wordID <-as.factor(word$wordID)
word$wordID <-as.factor(word$target)


for ( i in unique(word$target) ){
  #create a subset data 
  data_sub <- subset(word, target == i)
  output <- assign(paste("lmer",i,sep=""),tidy(lmer(w_n400 ~ condition*Anteriority*Hemisphere + (1|electrode) + (1|PartID),data=data_sub),"fixed"))
}


############################ERP plot#####################################

coi <- c("CZ")

data$PartID<-as.factor(data$PartID)
data$naming <-as.factor(data$naming)
condition_avg <- subset(data, subset = naming %in% c("Same Related","Different Related","Different Rerelated")) %>%
  filter(electrode %in% coi)%>%
  group_by(naming) %>% summarise(across(starts_with("X"), ~ mean(.x, na.rm = TRUE)))





plot <- condition_avg %>%
  pivot_longer(
    cols = starts_with("X"),
    names_to = "Timep",
    values_to = "Amp"
  ) 
plot$Timep <- gsub("X","",as.character(plot$Timep))
plot$Timep <- as.numeric(plot$Timep)


plot <- plot %>% 
  mutate(
       dif = Timep - 50,
       diff = (dif*4)-4,
       rank = (dense_rank(diff) - row_number()[dif== 0])) %>% 
  rename(
    Time = diff
  )



p<-ggplot(plot, aes(x = Time, y = Amp ,group=naming)) + 
  geom_line(aes(color=naming)) + 
  scale_x_continuous(breaks = c(seq(-200, 800,200))) +
  scale_y_continuous(breaks = c(seq(-20, 5, 10)))+
  scale_y_reverse()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) + 
  theme(
    legend.position = c(0.9, 0.1)
    #axis.text.x  = element_blank(),
    #axis.text.y  = element_blank()
  )+
  labs(
    x = "Time (ms)",
    y = "Amplitude (µV)"
  )
p

shift_axis <- function(p, x=0,y=0){
  g <- ggplotGrob(p)
  dummy <- data.frame(x=x,y=y)
  ax_y <- g[["grobs"]][g$layout$name == "axis-b"][[1]]
  ax_x <- g[["grobs"]][g$layout$name == "axis-l"][[1]]
  p + annotation_custom(grid::grobTree(ax_y, vp = grid::viewport(y=1, height=sum(ax_y$height))), 
                        ymax=y, ymin=y) +
    annotation_custom(grid::grobTree(ax_x, vp = grid::viewport(x=1, width = sum(ax_x$height))), 
                      xmax=x, xmin=x) +
    geom_hline(aes(yintercept=y), data = dummy) +
    geom_vline(aes(xintercept=x), data = dummy) +
    theme(axis.text.x = element_blank(), 
          axis.ticks.x=element_blank(),
          axis.text.y = element_blank(), 
          axis.ticks.y=element_blank())

  
}


shift_axis(p,x=0,y = 0)

#####################PLOT FEATURES################
#####CONC#######

coi <- c("CZ")

condition_avg_conc <- subset(data, subset = naming %in% c("Same Related","Different Related","Different Rerelated")) %>%
  filter(electrode %in% coi)%>%
  group_by(naming,concg2) %>% 
  summarise(across(starts_with("X"), ~ mean(.x, na.rm = TRUE)))


conc_diff <-subset(data, subset = naming %in% c("Same Related","Different Related","Different Rerelated")) %>%
  group_by(naming,concg2) %>% 
  summarise(electrode,n400)


write.csv(conc_diff, "/Users/agnesgao/Downloads/conc.csv", row.names=T)
  
                              
  
conc <- subset(data, subset = naming %in% c("Same Related","Different Related","Different Rerelated")) %>% 
  group_by(concg2,naming) %>% summarise(
    mean = mean(n400, na.rm = T),
    sd = sd(n400, na.rm = T),
    N = length(n400)) %>% mutate(se = sd/sqrt(N))

plot_conc <- condition_avg_conc %>%
  pivot_longer(
    cols = starts_with("X"),
    names_to = "Timep",
    values_to = "Amp"
  ) 
plot_conc$Timep <- gsub("X","",as.character(plot_conc$Timep))
plot_conc$Timep <- as.numeric(plot_conc$Timep)


plot_conc <- plot_conc %>% 
  mutate(
    dif = Timep - 50,
    diff = (dif*4)-4,
    rank = (dense_rank(diff) - row_number()[dif== 0])) %>% 
  rename(
    Time = diff
  )



#conc_diff <- plot_conc %>%
 # group_by(naming) %>% 
  #summarise (diff = Amp[concg2 ==" HIGH"] - Amp[concg2 =="low"])


p<-ggplot(plot_conc, aes(x = Time, y = Amp, group=interaction(naming,concg2))) + 
  geom_line(aes(color=naming,linetype = concg2)) + 
  scale_x_continuous(breaks = c(seq(-200, 800,200))) +
  scale_y_continuous(breaks = c(seq(-20, 5, 10)))+
  scale_y_reverse()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) + 
  theme(
    legend.position = c(0.9, 0.1),
    axis.text.x  = element_blank(),
    axis.text.y  = element_blank()
  )+
  labs(
    x = "Time (ms)",
    y = "Amplitude (µV)"
  )
p

shift_axis <- function(p, x=0,y=0){
  g <- ggplotGrob(p)
  dummy <- data.frame(x=x,y=y)
  ax_y <- g[["grobs"]][g$layout$name == "axis-b"][[1]]
  ax_x <- g[["grobs"]][g$layout$name == "axis-l"][[1]]
  p + annotation_custom(grid::grobTree(ax_y, vp = grid::viewport(y=1, height=sum(ax_y$height))), 
                        ymax=y, ymin=y) +
    annotation_custom(grid::grobTree(ax_x, vp = grid::viewport(x=1, width = sum(ax_x$height))), 
                      xmax=x, xmin=x) +
    geom_hline(aes(yintercept=y), data = dummy) +
    geom_vline(aes(xintercept=x), data = dummy) +
    theme(axis.text.x = element_blank(), 
          axis.ticks.x=element_blank(),
          axis.text.y = element_blank(), 
          axis.ticks.y=element_blank())
  
  
}


shift_axis(p,x=0,y = 0)




#####################ortho plot###########
coi <- c("PZ")

condition_avg_orth <- subset(data, subset = naming %in% c("Same Related","Different Related","Different Rerelated")) %>%
  filter(electrode %in% coi)%>%
  group_by(naming,orthog2) %>% 
  summarise(across(starts_with("X"), ~ mean(.x, na.rm = TRUE)))


plot_orth <- condition_avg_orth %>%
  pivot_longer(
    cols = starts_with("X"),
    names_to = "Timep",
    values_to = "Amp"
  ) 
plot_orth$Timep <- gsub("X","",as.character(plot_orth$Timep))
plot_orth$Timep <- as.numeric(plot_orth$Timep)


plot_orth <- plot_orth %>% 
  mutate(
    dif = Timep - 50,
    diff = (dif*4)-4,
    rank = (dense_rank(diff) - row_number()[dif== 0])) %>% 
  rename(
    Time = diff
  )



#conc_diff <- plot_conc %>%
# group_by(naming) %>% 
#summarise (diff = Amp[concg2 ==" HIGH"] - Amp[concg2 =="low"])


p<-ggplot(plot_orth, aes(x = Time, y = Amp, group=interaction(naming,orthog2))) + 
  geom_line(aes(color=naming,linetype = orthog2)) + 
  scale_x_continuous(breaks = c(seq(-200, 800,200))) +
  scale_y_continuous(breaks = c(seq(-20, 5, 10)))+
  scale_y_reverse()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) + 
  theme(
    legend.position = c(0.9, 0.1),
    axis.text.x  = element_blank(),
    axis.text.y  = element_blank()
  )+
  labs(
    x = "Time (ms)",
    y = "Amplitude (µV)"
  )
p

shift_axis <- function(p, x=0,y=0){
  g <- ggplotGrob(p)
  dummy <- data.frame(x=x,y=y)
  ax_y <- g[["grobs"]][g$layout$name == "axis-b"][[1]]
  ax_x <- g[["grobs"]][g$layout$name == "axis-l"][[1]]
  p + annotation_custom(grid::grobTree(ax_y, vp = grid::viewport(y=1, height=sum(ax_y$height))), 
                        ymax=y, ymin=y) +
    annotation_custom(grid::grobTree(ax_x, vp = grid::viewport(x=1, width = sum(ax_x$height))), 
                      xmax=x, xmin=x) +
    geom_hline(aes(yintercept=y), data = dummy) +
    geom_vline(aes(xintercept=x), data = dummy) +
    theme(axis.text.x = element_blank(), 
          axis.ticks.x=element_blank(),
          axis.text.y = element_blank(), 
          axis.ticks.y=element_blank())
  
  
}


shift_axis(p,x=0,y = 0)

######################phono plot###############
coi <- c("CZ")

condition_avg_phono <- subset(data, subset = naming %in% c("Same Related","Different Related","Different Rerelated")) %>%
  filter(electrode %in% coi)%>%
  group_by(naming,PHONO2_G) %>% 
  summarise(across(starts_with("X"), ~ mean(.x, na.rm = TRUE)))


plot_phono <- condition_avg_phono %>%
  pivot_longer(
    cols = starts_with("X"),
    names_to = "Timep",
    values_to = "Amp"
  ) 
plot_phono$Timep <- gsub("X","",as.character(plot_phono$Timep))
plot_phono$Timep <- as.numeric(plot_phono$Timep)


plot_phono <- plot_phono %>% 
  mutate(
    dif = Timep - 50,
    diff = (dif*4)-4,
    rank = (dense_rank(diff) - row_number()[dif== 0])) %>% 
  rename(
    Time = diff
  )



#conc_diff <- plot_conc %>%
# group_by(naming) %>% 
#summarise (diff = Amp[concg2 ==" HIGH"] - Amp[concg2 =="low"])


p<-ggplot(plot_phono, aes(x = Time, y = Amp, group=interaction(naming,PHONO2_G))) + 
  geom_line(aes(color=naming,linetype = PHONO2_G)) + 
  scale_x_continuous(breaks = c(seq(-200, 800,200))) +
  scale_y_continuous(breaks = c(seq(-20, 5, 10)))+
  scale_y_reverse()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) + 
  theme(
    legend.position = c(0.9, 0.1),
    axis.text.x  = element_blank(),
    axis.text.y  = element_blank()
  )+
  labs(
    x = "Time (ms)",
    y = "Amplitude (µV)"
  )
p

shift_axis <- function(p, x=0,y=0){
  g <- ggplotGrob(p)
  dummy <- data.frame(x=x,y=y)
  ax_y <- g[["grobs"]][g$layout$name == "axis-b"][[1]]
  ax_x <- g[["grobs"]][g$layout$name == "axis-l"][[1]]
  p + annotation_custom(grid::grobTree(ax_y, vp = grid::viewport(y=1, height=sum(ax_y$height))), 
                        ymax=y, ymin=y) +
    annotation_custom(grid::grobTree(ax_x, vp = grid::viewport(x=1, width = sum(ax_x$height))), 
                      xmax=x, xmin=x) +
    geom_hline(aes(yintercept=y), data = dummy) +
    geom_vline(aes(xintercept=x), data = dummy) +
    theme(axis.text.x = element_blank(), 
          axis.ticks.x=element_blank(),
          axis.text.y = element_blank(), 
          axis.ticks.y=element_blank())
  
  
}


shift_axis(p,x=0,y = 0)

##################leng plot######################

coi <- c("O2")

condition_avg_len <- subset(data, subset = naming %in% c("Same Related","Different Related","Different Rerelated")) %>%
  filter(electrode %in% coi)%>%
  group_by(naming,leng2) %>% 
  summarise(across(starts_with("X"), ~ mean(.x, na.rm = TRUE)))


plot_len <- condition_avg_len %>%
  pivot_longer(
    cols = starts_with("X"),
    names_to = "Timep",
    values_to = "Amp"
  ) 
plot_len$Timep <- gsub("X","",as.character(plot_len$Timep))
plot_len$Timep <- as.numeric(plot_len$Timep)


plot_len <- plot_len %>% 
  mutate(
    dif = Timep - 50,
    diff = (dif*4)-4,
    rank = (dense_rank(diff) - row_number()[dif== 0])) %>% 
  rename(
    Time = diff
  )
  
  #-200 - 200
  plot_len_sub <- plot_len %>% filter(Timep %in% (3:101))


p<-ggplot(plot_len_sub, aes(x = Time, y = Amp, group=interaction(naming,leng2))) + 
  geom_line(aes(color=naming,linetype = leng2)) + 
  scale_linetype_manual(values = c("SHORT"= "solid", 
                                   "LONG"= "longdash"))+
  scale_x_continuous(breaks = c(seq(-200, 400, 100))) +
  scale_y_continuous(breaks = c(seq(-20, 5, 10)))+
  scale_y_reverse()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) + 
  theme(
    legend.position = c(0.1, 0.2),
    axis.text.x  = element_blank(),
    axis.text.y  = element_blank()
  )+
  labs(
    x = "Time (ms)",
    y = "Amplitude (µV)"
  )
p

shift_axis <- function(p, x=0,y=0){
  g <- ggplotGrob(p)
  dummy <- data.frame(x=x,y=y)
  ax_y <- g[["grobs"]][g$layout$name == "axis-b"][[1]]
  ax_x <- g[["grobs"]][g$layout$name == "axis-l"][[1]]
  p + annotation_custom(grid::grobTree(ax_y, vp = grid::viewport(y=1, height=sum(ax_y$height))), 
                        ymax=y, ymin=y) +
    annotation_custom(grid::grobTree(ax_x, vp = grid::viewport(x=1, width = sum(ax_x$height))), 
                      xmax=x, xmin=x) +
    geom_hline(aes(yintercept=y), data = dummy) +
    geom_vline(aes(xintercept=x), data = dummy) +
    theme(axis.text.x = element_blank(), 
          axis.ticks.x=element_blank(),
          axis.text.y = element_blank(), 
          axis.ticks.y=element_blank())
  
  
}


shift_axis(p,x=0,y = 0)



#####################bar plot##################

bar.p1 <- subset(data, subset = naming %in% c("Same Related","Different Related","Different Rerelated")) %>% 
  group_by(naming) %>% summarise(
    mean = mean(p1, na.rm = T),
    sd = sd(p1, na.rm = T),
    N = length(p1)) %>% mutate(se = sd/sqrt(N))

ggplot(data = bar.p1, aes(x = naming, y = mean, fill=naming,palette = "jco")) + 
  geom_bar(stat="identity") + 
  geom_errorbar(width = .3, aes(ymin = mean - se, ymax = mean + se)) + 
  theme_classic()+
  geom_errorbar(width = .3, aes(ymin = mean - se, ymax = mean + se),position = position_dodge(), colour="black") + 
  labs( x = "Condition",
        y = "Mean Amplitude (µV)")

bar.p1.conc <- subset(data, subset = naming %in% c("Same Related","Different Related","Different Rerelated")) %>% 
  group_by(concg2,naming) %>% summarise(
    mean = mean(p1, na.rm = T),
    sd = sd(p1, na.rm = T),
    N = length(p1)) %>% mutate(se = sd/sqrt(N))

ggplot(data = bar.p1.conc, aes(x = concg2, y = mean, fill=naming,palette = "jco")) + 
  geom_bar(stat="identity") + 
  geom_errorbar(width = .3, aes(ymin = mean - se, ymax = mean + se)) + 
  theme_classic()+
  geom_errorbar(width = .3, aes(ymin = mean - se, ymax = mean + se),position = position_dodge(), colour="black") + 
  labs( x = "Concreteness",
        y = "Mean Amplitude (µV)")+
  facet_wrap(~naming)


bar.p2 <- subset(data, subset = naming %in% c("Same Related","Different Related","Different Rerelated")) %>% 
  group_by(naming) %>% summarise(
    mean = mean(p2, na.rm = T),
    sd = sd(p2, na.rm = T),
    N = length(p2)) %>% mutate(se = sd/sqrt(N))

ggplot(data = bar.p2, aes(x = naming, y = mean, fill=naming,palette = "jco")) + 
  geom_bar(stat="identity") + 
  geom_errorbar(width = .3, aes(ymin = mean - se, ymax = mean + se)) + 
  theme_classic()+
  geom_errorbar(width = .3, aes(ymin = mean - se, ymax = mean + se),position = position_dodge(), colour="black") + 
  labs( x = "Condition",
        y = "Mean Amplitude (µV)")



bar.n400 <- subset(data, subset = naming %in% c("Same Related","Different Related","Different Rerelated")) %>% 
  group_by(naming) %>% summarise(
    mean = mean(n400, na.rm = T),
    sd = sd(n400, na.rm = T),
    N = length(n400)) %>% mutate(se = sd/sqrt(N))

ggplot(data = bar.n400, aes(x = naming, y = mean, fill=naming,palette = "jco")) + 
  geom_bar(stat="identity") + 
  geom_errorbar(width = .3, aes(ymin = mean - se, ymax = mean + se)) + 
  theme_classic()+
  geom_errorbar(width = .3, aes(ymin = mean - se, ymax = mean + se),position = position_dodge(), colour="black") + 
  labs( x = "Condition",
        y = "N400 Mean Amplitude (µV)")

bar.conc.n400 <- subset(data, subset = naming %in% c("Same Related","Different Related","Different Rerelated")) %>% 
  group_by(concg2,naming) %>% summarise(
    mean = mean(n400, na.rm = T),
    sd = sd(n400, na.rm = T),
    N = length(n400)) %>% mutate(se = sd/sqrt(N))

ggplot(data = bar.conc.n400, aes(x = concg2, y = mean, fill=naming,palette = "jco")) + 
  geom_bar(stat="identity") + 
  geom_errorbar(width = .3, aes(ymin = mean - se, ymax = mean + se)) + 
  theme_classic()+
  geom_errorbar(width = .3, aes(ymin = mean - se, ymax = mean + se),position = position_dodge(), colour="black") + 
  labs( x = "Concreteness",
        y = "N400 Amplitude (µV)")+
  facet_wrap(~naming)

bar.conc.n250 <- subset(data, subset = naming %in% c("Same Related","Different Related","Different Rerelated")) %>% 
  group_by(concg2,naming) %>% summarise(
    mean = mean(n250, na.rm = T),
    sd = sd(n250, na.rm = T),
    N = length(n250)) %>% mutate(se = sd/sqrt(N))

ggplot(data = bar.conc.n250, aes(x = concg2, y = mean, fill=naming,palette = "jco")) + 
  geom_bar(stat="identity") + 
  geom_errorbar(width = .3, aes(ymin = mean - se, ymax = mean + se)) + 
  theme_classic()+
  geom_errorbar(width = .3, aes(ymin = mean - se, ymax = mean + se),position = position_dodge(), colour="black") + 
  labs( x = "Concreteness",
        y = "N250 Amplitude (µV)")+
  facet_wrap(~naming)



bar.ortho.n400 <- subset(data, subset = naming %in% c("Same Related","Different Related","Different Rerelated")) %>% 
  group_by(orthog2,naming) %>% summarise(
    mean = mean(n400, na.rm = T),
    sd = sd(n400, na.rm = T),
    N = length(n400)) %>% mutate(se = sd/sqrt(N))

ggplot(data = bar.ortho.n400, aes(x = orthog2, y = mean, fill=naming,palette = "jco")) + 
  geom_bar(stat="identity") + 
  geom_errorbar(width = .3, aes(ymin = mean - se, ymax = mean + se)) + 
  theme_classic()+
  geom_errorbar(width = .3, aes(ymin = mean - se, ymax = mean + se),position = position_dodge(), colour="black") + 
  labs( x = "Orthographic Neighborhood Size",
        y = "N400 Amplitude (µV)")+
  facet_wrap(~naming)


bar.ortho.p2 <- subset(data, subset = naming %in% c("Same Related","Different Related","Different Rerelated")) %>% 
  group_by(orthog2,naming) %>% summarise(
    mean = mean(p2, na.rm = T),
    sd = sd(p2, na.rm = T),
    N = length(p2)) %>% mutate(se = sd/sqrt(N))

ggplot(data = bar.ortho.p2, aes(x = orthog2, y = mean, fill=naming,palette = "jco")) + 
  geom_bar(stat="identity") + 
  geom_errorbar(width = .3, aes(ymin = mean - se, ymax = mean + se)) + 
  theme_classic()+
  geom_errorbar(width = .3, aes(ymin = mean - se, ymax = mean + se),position = position_dodge(), colour="black") + 
  labs( x = "Orthographic Neighborhood Size",
        y = "Mean Amplitude (µV)")+
  facet_wrap(~naming)


bar.phono.p2 <- subset(data, subset = naming %in% c("Same Related","Different Related","Different Rerelated")) %>% 
  group_by(PHONO2_G,naming) %>% summarise(
    mean = mean(p2, na.rm = T),
    sd = sd(p2, na.rm = T),
    N = length(p2)) %>% mutate(se = sd/sqrt(N))

ggplot(data = bar.phono.p2, aes(x = PHONO2_G, y = mean, fill=naming,palette = "jco")) + 
  geom_bar(stat="identity") + 
  geom_errorbar(width = .3, aes(ymin = mean - se, ymax = mean + se)) + 
  theme_classic()+
  geom_errorbar(width = .3, aes(ymin = mean - se, ymax = mean + se),position = position_dodge(), colour="black") + 
  labs( x = "Phonological Neighborhood Size",
        y = "Mean Amplitude (µV)")+
  facet_wrap(~naming)


bar.len.p1 <- subset(data, subset = naming %in% c("Same Related","Different Related","Different Rerelated")) %>% 
  group_by(leng2,naming) %>% summarise(
    mean = mean(p1, na.rm = T),
    sd = sd(p1, na.rm = T),
    N = length(p1)) %>% mutate(se = sd/sqrt(N))

ggplot(data = bar.len.p1, aes(x = leng2, y = mean, fill=naming,palette = "jco")) + 
  geom_bar(stat="identity") + 
  geom_errorbar(width = .3, aes(ymin = mean - se, ymax = mean + se)) + 
  theme_classic()+
  geom_errorbar(width = .3, aes(ymin = mean - se, ymax = mean + se),position = position_dodge(), colour="black") + 
  labs( x = "Word Length (char)",
        y = "P1 Amplitude (µV)")+
  facet_wrap(~naming)

bar.len.n400 <- subset(data, subset = naming %in% c("Same Related","Different Related","Different Rerelated")) %>% 
  group_by(leng2,naming) %>% summarise(
    mean = mean(n400, na.rm = T),
    sd = sd(n400, na.rm = T),
    N = length(n400)) %>% mutate(se = sd/sqrt(N))

ggplot(data = bar.len.n400, aes(x = leng2, y = mean, fill=naming,palette = "jco")) + 
  geom_bar(stat="identity") + 
  geom_errorbar(width = .3, aes(ymin = mean - se, ymax = mean + se)) + 
  theme_classic()+
  geom_errorbar(width = .3, aes(ymin = mean - se, ymax = mean + se),position = position_dodge(), colour="black") + 
  labs( x = "Word Length (char)",
        y = "N400 Amplitude (µV)")+
  facet_wrap(~naming)

########################## DW #####################
#p1
#main
p1.dw <-subset(data, subset = naming %in% c("Same Related","Different Related","Different Rerelated")) %>%
  group_by(naming, electrode) %>% 
  summarise(mean=mean(p1))

p1.dw.pred <-subset(p1.dw, subset = naming %in% c("Same Related"))
p1.dw.unpred <-subset(p1.dw, subset = naming %in% c("Different Related"))
p1.dw.unre <-subset(p1.dw, subset = naming %in% c("Different Rerelated"))


#order electrode
elec <- c("FP1","FP2","F3","F4","F7","F8","FC1","FC2","FC5","FC6","C3","C4","T3",
          "T4","CP1","CP2","CP5","CP6","P3","P4","T5","T6","O1","O2","AFZ","FZ","CZ","PZ","POZ")


p1.dw.pred$electrode <- factor(p1.dw.pred$electrode, levels=elec)
p1.dw.pred<-p1.dw.pred[order(p1.dw.pred$electrode),]

p1.dw.unpred$electrode <- factor(p1.dw.unpred$electrode, levels=elec)
p1.dw.unpred<-p1.dw.unpred[order(p1.dw.unpred$electrode),]

p1.dw.unre$electrode <- factor(p1.dw.unre$electrode, levels=elec)
p1.dw.unre<-p1.dw.unre[order(p1.dw.unre$electrode),]



p1.both.dw<-read.csv('/Users/agnesgao/Documents/p1_both_dw.csv', header=F,row.names = NULL)

fwrite(data.frame(p1.both.dw),"/Users/agnesgao/Downloads/p1_both_dw.csv",sep=",",row.names = FALSE,col.names = FALSE)


#n400
#all electrode

#unpredicted - predicted


n4.dw.main <-subset(data, subset = naming %in% c("Same Related","Different Related","Different Rerelated")) %>%
  group_by(naming, electrode) %>% 
  summarise(mean=mean(n400))

n4.dw.pred <-subset(n4.dw.main, subset = naming %in% c("Same Related"))
n4.dw.unpred <-subset(n4.dw.main, subset = naming %in% c("Different Related"))
n4.dw.unre <-subset(n4.dw.main, subset = naming %in% c("Different Rerelated"))

#order electrode
elec <- c("FP1","FP2","F3","F4","F7","F8","FC1","FC2","FC5","FC6","C3","C4","T3",
          "T4","CP1","CP2","CP5","CP6","P3","P4","T5","T6","O1","O2","AFZ","FZ","CZ","PZ","POZ")


n4.dw.pred$electrode <- factor(n4.dw.pred$electrode, levels=elec)
n4.dw.pred<-n4.dw.pred[order(n4.dw.pred$electrode),]

n4.dw.unpred$electrode <- factor(n4.dw.unpred$electrode, levels=elec)
n4.dw.unpred<-n4.dw.unpred[order(n4.dw.unpred$electrode),]

n4.dw.unre$electrode <- factor(n4.dw.unre$electrode, levels=elec)
n4.dw.unre<-n4.dw.unre[order(n4.dw.unre$electrode),]




fwrite(data.frame(namerelatedw),"/Users/agnesgao/Downloads/namerelate.csv",sep=",",row.names = FALSE,col.names = FALSE)

namerelatedw<-read.csv('/Users/agnesgao/Documents/namingandrelate_dw.csv', header=F,row.names = NULL)


#CONC
n4.dw.conc <-subset(data, subset = naming %in% c("Same Related","Different Related","Different Rerelated")) %>%
  group_by(naming, concg2,electrode) %>% 
  summarise(mean=mean(n400))

n4.dw.conc.pred <-subset(n4.dw.conc, subset = naming %in% c("Same Related"))
n4.dw.conc.unpred <-subset(n4.dw.conc, subset = naming %in% c("Different Related"))
n4.dw.conc.unre <-subset(n4.dw.conc, subset = naming %in% c("Different Rerelated"))


#order electrode
elec <- c("FP1","FP2","F3","F4","F7","F8","FC1","FC2","FC5","FC6","C3","C4","T3",
          "T4","CP1","CP2","CP5","CP6","P3","P4","T5","T6","O1","O2","AFZ","FZ","CZ","PZ","POZ")


n4.dw.conc.pred$electrode <- factor(n4.dw.conc.pred$electrode, levels=elec)
n4.dw.conc.pred<-n4.dw.conc.pred[order(n4.dw.conc.pred$electrode),]

n4.dw.conc.unpred$electrode <- factor(n4.dw.conc.unpred$electrode, levels=elec)
n4.dw.conc.unpred<-n4.dw.conc.unpred[order(n4.dw.conc.unpred$electrode),]

n4.dw.conc.unre$electrode <- factor(n4.dw.conc.unre$electrode, levels=elec)
n4.dw.conc.unre<-n4.dw.conc.unre[order(n4.dw.conc.unre$electrode),]

fwrite(data.frame(conc.unre.dw),"/Users/agnesgao/Downloads/conc_unre.csv",sep=",",row.names = FALSE,col.names = FALSE)

conc.unre.dw<-read.csv('/Users/agnesgao/Documents/conc_unre.csv', header=F,row.names = NULL)

#CONC-n1
p1.dw.conc <-subset(data, subset = naming %in% c("Same Related","Different Related","Different Rerelated")) %>%
  group_by(naming, concg2,electrode) %>% 
  summarise(mean=mean(p1))

p1.dw.conc.pred <-subset(p1.dw.conc, subset = naming %in% c("Same Related"))
p1.dw.conc.unpred <-subset(p1.dw.conc, subset = naming %in% c("Different Related"))
p1.dw.conc.unre <-subset(p1.dw.conc, subset = naming %in% c("Different Rerelated"))


#order electrode
elec <- c("FP1","FP2","F3","F4","F7","F8","FC1","FC2","FC5","FC6","C3","C4","T3",
          "T4","CP1","CP2","CP5","CP6","P3","P4","T5","T6","O1","O2","AFZ","FZ","CZ","PZ","POZ")


p1.dw.conc.pred$electrode <- factor(p1.dw.conc.pred$electrode, levels=elec)
p1.dw.conc.pred<-p1.dw.conc.pred[order(p1.dw.conc.pred$electrode),]

p1.dw.conc.unpred$electrode <- factor(p1.dw.conc.unpred$electrode, levels=elec)
p1.dw.conc.unpred<-p1.dw.conc.unpred[order(p1.dw.conc.unpred$electrode),]

p1.dw.conc.unre$electrode <- factor(p1.dw.conc.unre$electrode, levels=elec)
p1.dw.conc.unre<-p1.dw.conc.unre[order(p1.dw.conc.unre$electrode),]


p1.dw.conc.pred<-read.csv('/Users/agnesgao/Documents/conc_p1_pred.csv', header=F,row.names = NULL)

fwrite(data.frame(p1.dw.conc.pred),"/Users/agnesgao/Downloads/conc_p1_pred.csv",sep=",",row.names = FALSE,col.names = FALSE)



#ortho
n4.dw.orth <-subset(data, subset = naming %in% c("Same Related","Different Related","Different Rerelated")) %>%
  group_by(naming, orthog2,electrode) %>% 
  summarise(mean=mean(n400))

n4.dw.orth.pred <-subset(n4.dw.orth, subset = naming %in% c("Same Related"))
n4.dw.orth.unpred <-subset(n4.dw.orth, subset = naming %in% c("Different Related"))
n4.dw.orth.unre <-subset(n4.dw.orth, subset = naming %in% c("Different Rerelated"))


#order electrode
elec <- c("FP1","FP2","F3","F4","F7","F8","FC1","FC2","FC5","FC6","C3","C4","T3",
          "T4","CP1","CP2","CP5","CP6","P3","P4","T5","T6","O1","O2","AFZ","FZ","CZ","PZ","POZ")


n4.dw.orth.pred$electrode <- factor(n4.dw.orth.pred$electrode, levels=elec)
n4.dw.orth.pred<-n4.dw.orth.pred[order(n4.dw.orth.pred$electrode),]

n4.dw.orth.unpred$electrode <- factor(n4.dw.orth.unpred$electrode, levels=elec)
n4.dw.orth.unpred<-n4.dw.orth.unpred[order(n4.dw.orth.unpred$electrode),]

n4.dw.orth.unre$electrode <- factor(n4.dw.orth.unre$electrode, levels=elec)
n4.dw.orth.unre<-n4.dw.orth.unre[order(n4.dw.orth.unre$electrode),]



orth.unre.dw<-read.csv('/Users/agnesgao/Documents/orth_unre.csv', header=F,row.names = NULL)

fwrite(data.frame(orth.unre.dw),"/Users/agnesgao/Downloads/orth_unre.csv",sep=",",row.names = FALSE,col.names = FALSE)

#ortho-p2
p2.dw.orth <-subset(data, subset = naming %in% c("Same Related","Different Related","Different Rerelated")) %>%
  group_by(naming, orthog2,electrode) %>% 
  summarise(mean=mean(p2))

p2.dw.orth.pred <-subset(p2.dw.orth, subset = naming %in% c("Same Related"))
p2.dw.orth.unpred <-subset(p2.dw.orth, subset = naming %in% c("Different Related"))
p2.dw.orth.unre <-subset(p2.dw.orth, subset = naming %in% c("Different Rerelated"))


#order electrode
elec <- c("FP1","FP2","F3","F4","F7","F8","FC1","FC2","FC5","FC6","C3","C4","T3",
          "T4","CP1","CP2","CP5","CP6","P3","P4","T5","T6","O1","O2","AFZ","FZ","CZ","PZ","POZ")


p2.dw.orth.pred$electrode <- factor(p2.dw.orth.pred$electrode, levels=elec)
p2.dw.orth.pred<-p2.dw.orth.pred[order(p2.dw.orth.pred$electrode),]

p2.dw.orth.unpred$electrode <- factor(p2.dw.orth.unpred$electrode, levels=elec)
p2.dw.orth.unpred<-p2.dw.orth.unpred[order(p2.dw.orth.unpred$electrode),]

p2.dw.orth.unre$electrode <- factor(p2.dw.orth.unre$electrode, levels=elec)
p2.dw.orth.unre<-p2.dw.orth.unre[order(p2.dw.orth.unre$electrode),]



p2.orth.pred.dw<-read.csv('/Users/agnesgao/Documents/p2_ortho_pred.csv', header=F,row.names = NULL)

fwrite(data.frame(p2.orth.pred.dw),"/Users/agnesgao/Downloads/ortho_p2_pred.csv",sep=",",row.names = FALSE,col.names = FALSE)

#phono
p2.dw.pho <-subset(data, subset = naming %in% c("Same Related","Different Related","Different Rerelated")) %>%
  group_by(naming, PHONO2_G,electrode) %>% 
  summarise(mean=mean(p2))

p2.dw.pho.pred <-subset(p2.dw.pho, subset = naming %in% c("Same Related"))
p2.dw.pho.unpred <-subset(p2.dw.pho, subset = naming %in% c("Different Related"))
p2.dw.pho.unre <-subset(p2.dw.pho, subset = naming %in% c("Different Rerelated"))


#order electrode
elec <- c("FP1","FP2","F3","F4","F7","F8","FC1","FC2","FC5","FC6","C3","C4","T3",
          "T4","CP1","CP2","CP5","CP6","P3","P4","T5","T6","O1","O2","AFZ","FZ","CZ","PZ","POZ")


p2.dw.pho.pred$electrode <- factor(p2.dw.pho.pred$electrode, levels=elec)
p2.dw.pho.pred<-p2.dw.pho.pred[order(p2.dw.pho.pred$electrode),]

p2.dw.pho.unpred$electrode <- factor(p2.dw.pho.unpred$electrode, levels=elec)
p2.dw.pho.unpred<-p2.dw.pho.unpred[order(p2.dw.pho.unpred$electrode),]

p2.dw.pho.unre$electrode <- factor(p2.dw.pho.unre$electrode, levels=elec)
p2.dw.pho.unre<-p2.dw.pho.unre[order(p2.dw.pho.unre$electrode),]



pho.pred.dw<-read.csv('/Users/agnesgao/Documents/p2_phono_pred.csv', header=F,row.names = NULL)

fwrite(data.frame(pho.pred.dw),"/Users/agnesgao/Downloads/p2_pho_pred.csv",sep=",",row.names = FALSE,col.names = FALSE)


#n4
#len
n4.dw.len <-subset(data, subset = naming %in% c("Same Related","Different Related","Different Rerelated")) %>%
  group_by(naming, leng2,electrode) %>% 
  summarise(mean=mean(n400))

n4.dw.len.pred <-subset(n4.dw.len, subset = naming %in% c("Same Related"))
n4.dw.len.unpred <-subset(n4.dw.len, subset = naming %in% c("Different Related"))
n4.dw.len.unre <-subset(n4.dw.len, subset = naming %in% c("Different Rerelated"))


#order electrode
elec <- c("FP1","FP2","F3","F4","F7","F8","FC1","FC2","FC5","FC6","C3","C4","T3",
          "T4","CP1","CP2","CP5","CP6","P3","P4","T5","T6","O1","O2","AFZ","FZ","CZ","PZ","POZ")


n4.dw.len.pred$electrode <- factor(n4.dw.len.pred$electrode, levels=elec)
n4.dw.len.pred<-n4.dw.len.pred[order(n4.dw.len.pred$electrode),]

n4.dw.len.unpred$electrode <- factor(n4.dw.len.unpred$electrode, levels=elec)
n4.dw.len.unpred<-n4.dw.len.unpred[order(n4.dw.len.unpred$electrode),]

n4.dw.len.unre$electrode <- factor(n4.dw.len.unre$electrode, levels=elec)
n4.dw.len.unre<-n4.dw.len.unre[order(n4.dw.len.unre$electrode),]



len.pred.dw<-read.csv('/Users/agnesgao/Documents/len_n4_pred.csv', header=F,row.names = NULL)

fwrite(data.frame(len.pred.dw),"/Users/agnesgao/Downloads/len_n4_pred.csv",sep=",",row.names = FALSE,col.names = FALSE)


#p1
#len
p1.dw.len <-subset(data, subset = naming %in% c("Same Related","Different Related","Different Rerelated")) %>%
  group_by(naming, leng2,electrode) %>% 
  summarise(mean=mean(p1))

p1.dw.len.pred <-subset(p1.dw.len, subset = naming %in% c("Same Related"))
p1.dw.len.unpred <-subset(p1.dw.len, subset = naming %in% c("Different Related"))
p1.dw.len.unre <-subset(p1.dw.len, subset = naming %in% c("Different Rerelated"))


#order electrode
elec <- c("FP1","FP2","F3","F4","F7","F8","FC1","FC2","FC5","FC6","C3","C4","T3",
          "T4","CP1","CP2","CP5","CP6","P3","P4","T5","T6","O1","O2","AFZ","FZ","CZ","PZ","POZ")


p1.dw.len.pred$electrode <- factor(p1.dw.len.pred$electrode, levels=elec)
p1.dw.len.pred<-p1.dw.len.pred[order(p1.dw.len.pred$electrode),]

p1.dw.len.unpred$electrode <- factor(p1.dw.len.unpred$electrode, levels=elec)
p1.dw.len.unpred<-p1.dw.len.unpred[order(p1.dw.len.unpred$electrode),]

p1.dw.len.unre$electrode <- factor(p1.dw.len.unre$electrode, levels=elec)
p1.dw.len.unre<-p1.dw.len.unre[order(p1.dw.len.unre$electrode),]



len.unre.dw<-read.csv('/Users/agnesgao/Documents/len_unre.csv', header=F,row.names = NULL)

fwrite(data.frame(len.unre.dw),"/Users/agnesgao/Downloads/len_unre.csv",sep=",",row.names = FALSE,col.names = FALSE)

##p2 main
#p1
#len
p2.dw <-subset(data, subset = naming %in% c("Same Related","Different Related","Different Rerelated")) %>%
  group_by(naming, electrode) %>% 
  summarise(mean=mean(p2))

p2.dw.pred <-subset(p2.dw, subset = naming %in% c("Same Related"))
p2.dw.unpred <-subset(p2.dw, subset = naming %in% c("Different Related"))
p2.dw.unre <-subset(p2.dw, subset = naming %in% c("Different Rerelated"))


#order electrode
elec <- c("FP1","FP2","F3","F4","F7","F8","FC1","FC2","FC5","FC6","C3","C4","T3",
          "T4","CP1","CP2","CP5","CP6","P3","P4","T5","T6","O1","O2","AFZ","FZ","CZ","PZ","POZ")


p2.dw.pred$electrode <- factor(p2.dw.pred$electrode, levels=elec)
p2.dw.pred<-p2.dw.pred[order(p2.dw.pred$electrode),]

p2.dw.unpred$electrode <- factor(p2.dw.unpred$electrode, levels=elec)
p2.dw.unpred<-p2.dw.unpred[order(p2.dw.unpred$electrode),]

p2.dw.unre$electrode <- factor(p2.dw.unre$electrode, levels=elec)
p2.dw.unre<-p2.dw.unre[order(p2.dw.unre$electrode),]



p2.both.dw<-read.csv('/Users/agnesgao/Documents/p2_both_dw.csv', header=F,row.names = NULL)

fwrite(data.frame(p2.both.dw),"/Users/agnesgao/Downloads/p2_both_dw.csv",sep=",",row.names = FALSE,col.names = FALSE)

##### plot dw ##############
#length
coi <- c("PZ")

condition_avg_lendw <- subset(data, subset = naming %in% c("Same Related","Different Related","Different Rerelated")) %>%
  filter(electrode %in% coi)%>%
  group_by(naming,leng2) %>% 
  summarise(across(starts_with("X"), ~ mean(.x, na.rm = TRUE)))


plot_lendw <- condition_avg_lendw %>%
  pivot_longer(
    cols = starts_with("X"),
    names_to = "Timep",
    values_to = "Amp"
  ) 
plot_lendw$Timep <- gsub("X","",as.character(plot_lendw$Timep))
plot_lendw$Timep <- as.numeric(plot_lendw$Timep)


plot_lendw <- plot_lendw %>% 
  mutate(
    dif = Timep - 50,
    diff = (dif*4)-4,
    rank = (dense_rank(diff) - row_number()[dif== 0])) %>% 
  rename(
    Time = diff
  )

#-200 - 200
plot_lendw_sub <- plot_lendw %>% filter(Timep %in% (1:101))

len_diff <- subset(plot_lendw, subset = naming %in% c("Same Related","Different Related","Different Rerelated")) %>% 
 group_by(naming, Timep,diff) %>% 
summarise (diff_amp = Amp[leng2 =="LONG"] - Amp[leng2 =="SHORT"])


p<-ggplot(len_diff, aes(x = diff, y = diff_amp, group=naming)) + 
  geom_line(aes(color=naming)) + 
  scale_x_continuous(breaks = c(seq(-200, 400,100))) +
  scale_y_continuous(breaks = c(seq(-20, 5, 10)))+
  scale_y_reverse()+
  scale_color_manual(values=c("#CC0000", "#9900CC","#FF9900"))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) + 
  theme(
    legend.position = c(0.9, 0.1),
    axis.text.x  = element_blank(),
    axis.text.y  = element_blank()
  )+
  labs(
    x = "Time (ms)",
    y = "Amplitude (µV)"
  )
p

shift_axis <- function(p, x=0,y=0){
  g <- ggplotGrob(p)
  dummy <- data.frame(x=x,y=y)
  ax_y <- g[["grobs"]][g$layout$name == "axis-b"][[1]]
  ax_x <- g[["grobs"]][g$layout$name == "axis-l"][[1]]
  p + annotation_custom(grid::grobTree(ax_y, vp = grid::viewport(y=1, height=sum(ax_y$height))), 
                        ymax=y, ymin=y) +
    annotation_custom(grid::grobTree(ax_x, vp = grid::viewport(x=1, width = sum(ax_x$height))), 
                      xmax=x, xmin=x) +
    geom_hline(aes(yintercept=y), data = dummy) +
    geom_vline(aes(xintercept=x), data = dummy) +
    theme(axis.text.x = element_blank(), 
          axis.ticks.x=element_blank(),
          axis.text.y = element_blank(), 
          axis.ticks.y=element_blank())
  
  
}


shift_axis(p,x=0,y = 0)


#CONC
coi <- c("CZ")

condition_avg_concdw <- subset(data, subset = naming %in% c("Same Related","Different Related","Different Rerelated")) %>%
  filter(electrode %in% coi)%>%
  group_by(naming,concg2) %>% 
  summarise(across(starts_with("X"), ~ mean(.x, na.rm = TRUE)))


plot_concdw <- condition_avg_concdw %>%
  pivot_longer(
    cols = starts_with("X"),
    names_to = "Timep",
    values_to = "Amp"
  ) 
plot_concdw$Timep <- gsub("X","",as.character(plot_concdw$Timep))
plot_concdw$Timep <- as.numeric(plot_concdw$Timep)


plot_concdw <- plot_concdw %>% 
  mutate(
    dif = Timep - 50,
    diff = (dif*4)-4,
    rank = (dense_rank(diff) - row_number()[dif== 0])) %>% 
  rename(
    Time = diff
  )


conc_diff <- subset(plot_concdw, subset = naming %in% c("Same Related","Different Related","Different Rerelated")) %>% 
  group_by(naming, Timep,diff) %>% 
  summarise (diff_amp = Amp[concg2 =="HIGH"] - Amp[concg2 =="LOW"])


p<-ggplot(conc_diff, aes(x = diff, y = diff_amp, group=naming)) + 
  geom_line(aes(color=naming)) + 
  scale_x_continuous(breaks = c(seq(-200, 996,200))) +
  scale_y_continuous(breaks = c(seq(-20, 5, 10)))+
  scale_color_manual(values=c("#CC0000", "#9900CC","#FF9900"))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) + 
  theme(
    legend.position = c(1, 1)
    #axis.text.x  = element_blank(),
    #axis.text.y  = element_blank()
  )+
  labs(
    x = "Time (ms)",
    y = "Amplitude (µV)"
  )
p

shift_axis <- function(p, x=0,y=0){
  g <- ggplotGrob(p)
  dummy <- data.frame(x=x,y=y)
  ax_y <- g[["grobs"]][g$layout$name == "axis-b"][[1]]
  ax_x <- g[["grobs"]][g$layout$name == "axis-l"][[1]]
  p + annotation_custom(grid::grobTree(ax_y, vp = grid::viewport(y=1, height=sum(ax_y$height))), 
                        ymax=y, ymin=y) +
    annotation_custom(grid::grobTree(ax_x, vp = grid::viewport(x=1, width = sum(ax_x$height))), 
                      xmax=x, xmin=x) +
    geom_hline(aes(yintercept=y), data = dummy) +
    geom_vline(aes(xintercept=x), data = dummy) 
    #theme(axis.text.x = element_blank(), 
          #axis.ticks.x=element_blank(),
     #     axis.text.y = element_blank(), 
          #axis.ticks.y=element_blank())
  
  
}


shift_axis(p,x=0,y = 0)



#ortho
coi <- c("PZ")

condition_avg_orthodw <- subset(data, subset = naming %in% c("Same Related","Different Related","Different Rerelated")) %>%
  filter(electrode %in% coi)%>%
  group_by(naming,orthog2) %>% 
  summarise(across(starts_with("X"), ~ mean(.x, na.rm = TRUE)))


plot_orthodw <- condition_avg_orthodw %>%
  pivot_longer(
    cols = starts_with("X"),
    names_to = "Timep",
    values_to = "Amp"
  ) 
plot_orthodw$Timep <- gsub("X","",as.character(plot_orthodw$Timep))
plot_orthodw$Timep <- as.numeric(plot_orthodw$Timep)


plot_orthodw <- plot_orthodw %>% 
  mutate(
    dif = Timep - 50,
    diff = (dif*4)-4,
    rank = (dense_rank(diff) - row_number()[dif== 0])) %>% 
  rename(
    Time = diff
  )


ortho_diff <- subset(plot_orthodw, subset = naming %in% c("Same Related","Different Related","Different Rerelated")) %>% 
  group_by(naming, Timep,Time) %>% 
  summarise (dw = Amp[orthog2 =="HIGH"] - Amp[orthog2 =="LOW"])


p<-ggplot(ortho_diff, aes(x = Time, y = dw, group=naming)) + 
  geom_line(aes(color=naming)) + 
  scale_x_continuous(breaks = c(seq(-200, 996,200))) +
  scale_y_continuous(breaks = c(seq(-20, 5, 10)))+
  scale_y_reverse()+
  scale_color_manual(values=c("#CC0000", "#9900CC","#FF9900"))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) + 
  theme(
    legend.position = c(0.9, 0.9),
    axis.text.x  = element_blank(),
    axis.text.y  = element_blank()
  )+
  labs(
    x = "Time (ms)",
    y = "Amplitude (µV)"
  )
p

shift_axis <- function(p, x=0,y=0){
  g <- ggplotGrob(p)
  dummy <- data.frame(x=x,y=y)
  ax_y <- g[["grobs"]][g$layout$name == "axis-b"][[1]]
  ax_x <- g[["grobs"]][g$layout$name == "axis-l"][[1]]
  p + annotation_custom(grid::grobTree(ax_y, vp = grid::viewport(y=1, height=sum(ax_y$height))), 
                        ymax=y, ymin=y) +
    annotation_custom(grid::grobTree(ax_x, vp = grid::viewport(x=1, width = sum(ax_x$height))), 
                      xmax=x, xmin=x) +
    geom_hline(aes(yintercept=y), data = dummy) +
    geom_vline(aes(xintercept=x), data = dummy) +
    theme(axis.text.x = element_blank(), 
          axis.ticks.x=element_blank(),
          axis.text.y = element_blank(), 
          axis.ticks.y=element_blank())
  
  
}


shift_axis(p,x=0,y = 0)


#phono
coi <- c("PZ")

condition_avg_phodw <- subset(data, subset = naming %in% c("Same Related","Different Related","Different Rerelated")) %>%
  filter(electrode %in% coi)%>%
  group_by(naming,PHONO2_G) %>% 
  summarise(across(starts_with("X"), ~ mean(.x, na.rm = TRUE)))


plot_phodw <- condition_avg_phodw %>%
  pivot_longer(
    cols = starts_with("X"),
    names_to = "Timep",
    values_to = "Amp"
  ) 
plot_phodw$Timep <- gsub("X","",as.character(plot_phodw$Timep))
plot_phodw$Timep <- as.numeric(plot_phodw$Timep)


plot_phodw <- plot_phodw %>% 
  mutate(
    dif = Timep - 50,
    diff = (dif*4)-4,
    rank = (dense_rank(diff) - row_number()[dif== 0])) %>% 
  rename(
    Time = diff
  )


pho_diff <- subset(plot_phodw, subset = naming %in% c("Same Related","Different Related","Different Rerelated")) %>% 
  group_by(naming, Timep,Time) %>% 
  summarise (dw = Amp[PHONO2_G =="HIGH"] - Amp[PHONO2_G =="LOW"])


p<-ggplot(pho_diff, aes(x = Time, y = dw, group=naming)) + 
  geom_line(aes(color=naming)) + 
  scale_x_continuous(breaks = c(seq(-200, 996,200))) +
  scale_y_continuous(breaks = c(seq(-20, 5, 10)))+
  scale_y_reverse()+
  scale_color_manual(values=c("#CC0000", "#9900CC","#FF9900"))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) + 
  theme(
    legend.position = c(0.9, 0.9),
    axis.text.x  = element_blank(),
    axis.text.y  = element_blank()
  )+
  labs(
    x = "Time (ms)",
    y = "Amplitude (µV)"
  )
p

shift_axis <- function(p, x=0,y=0){
  g <- ggplotGrob(p)
  dummy <- data.frame(x=x,y=y)
  ax_y <- g[["grobs"]][g$layout$name == "axis-b"][[1]]
  ax_x <- g[["grobs"]][g$layout$name == "axis-l"][[1]]
  p + annotation_custom(grid::grobTree(ax_y, vp = grid::viewport(y=1, height=sum(ax_y$height))), 
                        ymax=y, ymin=y) +
    annotation_custom(grid::grobTree(ax_x, vp = grid::viewport(x=1, width = sum(ax_x$height))), 
                      xmax=x, xmin=x) +
    geom_hline(aes(yintercept=y), data = dummy) +
    geom_vline(aes(xintercept=x), data = dummy) +
    theme(axis.text.x = element_blank(), 
          axis.ticks.x=element_blank(),
          axis.text.y = element_blank(), 
          axis.ticks.y=element_blank())
  
  
}


shift_axis(p,x=0,y = 0)
