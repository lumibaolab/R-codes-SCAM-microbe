#all codes
library(lme4)
library(emmeans)
library(lmerTest)
library(multcompView)
library(multcomp)
library(plyr)
library(emmeans)
library(pbkrtest)
library(broom.mixed)
#=================================================================#
#OVERALL DATA, full model test across all Schoenoplectus plants

trait.d<-read.csv(file="All_data_plant_measurement_less_samp_corrected.csv")
trait.d[trait.d==0]<-NA
trait.d$Genotype<-as.factor(trait.d$Genotype)
options(contrasts = c("contr.sum","contr.poly"))
#GLB
lm.mod.grn2<-lmer(log(greenlf_biom)~Salinity*Mbiome.trt*Site+Cohort+initial_wt_grams+(1|Genotype), 
                  data=trait.d,na.action=na.omit)
anova(lm.mod.grn2, dff="Satterthwaite")
#AG
agb.mod0<-lmer(log(AG_biom)~Salinity*Mbiome.trt*Site+Cohort+initial_wt_grams+(1|Genotype), 
               data=trait.d,na.action=na.omit)
anova(agb.mod0, dff="Satterthwaite") # get the analysis of variance table
#BG
bgb.mod0<-lmer(log(BG_biom)~Salinity*Mbiome.trt*Site+Cohort+initial_wt_grams+(1|Genotype), 
               data=trait.d,na.action=na.omit)
anova(bgb.mod0, dff="Satterthwaite")
#Total Biomass
trait.d$total_biom<-trait.d$AG_biom+trait.d$BG_biom
total_b.mod01<-lmer(log(total_biom)~Salinity*Mbiome.trt*Site+Cohort+initial_wt_grams+(1|Genotype), 
                    data=trait.d,na.action=na.omit)
anova(total_b.mod01, dff="Satterthwaite")
#Total Biomass
trait.d$total_biom<-trait.d$AG_biom+trait.d$BG_biom
total_b.mod01<-lmer(log(total_biom)~Salinity*Mbiome.trt*Site+Cohort+initial_wt_grams+(1|Genotype), 
                    data=trait.d,na.action=na.omit)
anova(total_b.mod01, dff="Satterthwaite")

#root-to-shoot ratio
trait.d$rtst<-trait.d$BG_biom/trait.d$AG_biom
rtst.mod1<-lmer(log(rtst)~Salinity*Mbiome.trt*Site+Cohort+initial_wt_grams+(1|Genotype), 
                data=trait.d,na.action=na.omit)
anova(rtst.mod1, dff="Satterthwaite")
#significant interaction so
#Post-hoc tests:
emm<-emmeans(rtst.mod1, pairwise~Salinity:Mbiome.trt:Site)
emm
summary(emm, infer = TRUE, null = log(35), type = "response")
all.ht.cld<-multcomp::cld(emm$emmean,alpha=0.05, Letters=letters, adjust="tukey",type="response")
all.ht.cld

#height
height.mod01<-lmer(log(final_aveheight)~Salinity*Mbiome.trt*Site+Cohort+initial_wt_grams+(1|Genotype), 
                   data=trait.d,na.action=na.omit)
anova(height.mod01, dff="Satterthwaite")
#Post-hoc tests:
emm<-emmeans(height.mod01,pairwise~Salinity:Mbiome.trt:Site,adjust="tukey", alpha=0.05, Letters=letters,type="response")
emm
summary(emm, infer = TRUE, null = log(35), type = "response")
all.ht.cld<-multcomp::cld(emm$emmean,alpha=0.05, Letters=letters, adjust="tukey",type="response")
all.ht.cld
#visualize nature of interaction
emmip(height.mod01, Salinity ~ Mbiome.trt | Site)

trait.d$plantsize<-trait.d$stemdiam_ave*trait.d$final_aveheight
size.mod0<-lmer(log(plantsize)~Salinity*Mbiome.trt*Site+Cohort+initial_wt_grams+(1|Genotype), 
                data=trait.d,na.action=na.omit)
anova(size.mod0, dff="Satterthwaite")

#Stem number 
trait.d<-read.csv(file="All_data_plant_measurement_less_samp_corrected.csv")
trait.d[trait.d==0]<-NA
trait.d$Genotype<-as.factor(trait.d$Genotype)
options(contrasts = c("contr.sum","contr.poly"))
#transform
trait.d$sqstem_num<-sqrt(trait.d$stem_num)
stem.num.mod1<-lmer(sqrt(trait.d$stem_num)~Salinity*Mbiome.trt*Site+Cohort+initial_wt_grams+(1|Genotype), 
                    data=trait.d,na.action=na.omit)
anova(stem.num.mod1, dff="Satterthwaite") # get the analysis of variance table
#Post-hoc tests:
#need to regrid because sqrt transformation is used
sn.rg <- ref_grid(stem.num.mod1)
emm.sn.regrid<-emmeans(regrid(sn.rg), pairwise~Salinity:Mbiome.trt:Site,adjust="tukey", alpha=0.05, transform = "log",type="response")
emm.sn.regrid
all.sn.cld<-multcomp::cld(emm.sn.regrid$emmean,alpha=0.05, Letters=letters, adjust="tukey",type="response")

#SD
#area of each pot = 36 cm2
trait.d$stemden<-trait.d$stem_num/36
stemden.mod1<-lmer(log(stemden)~Salinity*Mbiome.trt*Site+Cohort+initial_wt_grams+(1|Genotype), 
                   data=trait.d,na.action=na.omit)
anova(stemden.mod1, dff="Satterthwaite")
#Post-hoc tests:
emm3<-emmeans(stemden.mod1,pairwise~Salinity:Mbiome.trt:Site,adjust="tukey", alpha=0.05, Letters=letters,type="response")
emm3
#cld(emm2$emmean)
all.ht.cld<-multcomp::cld(emm3$emmean,alpha=0.05, Letters=letters, adjust="tukey",type="response")
all.ht.cld

#Stem diameter
diam.mod0<-lmer(log(stemdiam_ave)~Salinity*Mbiome.trt*Site+Cohort+initial_wt_grams+(1|Genotype), 
                data=trait.d,na.action=na.omit)
anova(diam.mod0, dff="Satterthwaite")
#Post-hoc tests:
emm4<-emmeans(diam.mod0,pairwise~Salinity:Mbiome.trt:Site,adjust="tukey", alpha=0.05, Letters=letters,type="response")
emm4
#cld(emm4$emmean)
all.ht.cld<-multcomp::cld(emm4$emmean,alpha=0.05, Letters=letters, adjust="tukey",type="response")
all.ht.cld

#*****************************#
#getting BLUPs / reaction norm
#SELLMAN
#Plot only those where Salinity|Genotype model was significant
trait.d<-read.csv(file="All_data_plant_measurement_less_samp_corrected.csv")
sellman<-trait.d[trait.d$Site=="Sellman",]
#replace 0 with NAs
sellman[sellman==0]<-NA
sellman$Genotype<-as.factor(sellman$Genotype)
options(contrasts = c("contr.sum","contr.poly"))

#GLB
lm.mod.grn1<-lmer(log(greenlf_biom)~Salinity*Mbiome.trt+Cohort+initial_wt_grams+(1|Genotype), 
                  data=sellman,na.action=na.omit)
#vary slope then compare to model with varying intercept
lm.mod.grn2<-lmer(log(greenlf_biom)~Salinity*Mbiome.trt+Cohort+initial_wt_grams+(Salinity|Genotype), 
                  data=sellman,na.action=na.omit)
lm.mod.grn3<-lmer(log(greenlf_biom)~Salinity*Mbiome.trt+Cohort+initial_wt_grams+(Mbiome.trt|Salinity:Genotype), 
                  data=sellman,na.action=na.omit)
anova(lm.mod.grn1,lm.mod.grn3,lm.mod.grn2, test="FIML") # test if genotype is significant

#AG
agb.mod0<-lmer(log(AG_biom)~Salinity*Mbiome.trt+Cohort+initial_wt_grams+(1|Genotype), 
               data=sellman,na.action=na.omit)
agb.mod1<-lmer(log(AG_biom)~Salinity*Mbiome.trt+Cohort+initial_wt_grams+(Salinity|Genotype), 
               data=sellman,na.action=na.omit)
agb.mod2<-lmer(log(AG_biom)~Salinity*Mbiome.trt+Cohort+initial_wt_grams+(Mbiome.trt|Salinity:Genotype), 
               data=sellman,na.action=na.omit)
anova(agb.mod0,agb.mod1,agb.mod2, method="FIML")
#model 2 is significant
#get BLUPS
blups_agb<- broom.mixed::tidy(agb.mod2, effects="ran_vals") #back transform
blups_agb #get estimate and std.error
write.table(blups_agb, file="AG_blups_sellman_28Oct.csv",sep=',')
bgb.mod0<-lmer(log(BG_biom)~Salinity*Mbiome.trt+Cohort+initial_wt_grams+(1|Genotype), 
               data=sellman,na.action=na.omit)
bgb.mod2<-lmer(log(BG_biom)~Salinity*Mbiome.trt+Cohort+initial_wt_grams+(Salinity|Genotype), 
               data=sellman,na.action=na.omit)
bgb.mod1<-lmer(log(BG_biom)~Salinity*Mbiome.trt+Cohort+initial_wt_grams+(Mbiome.trt|Salinity:Genotype), 
               data=sellman,na.action=na.omit)
anova(bgb.mod0,bgb.mod1,bgb.mod2)
#get BLUPS
blups_BG<- broom.mixed::tidy(bgb.mod1, effects="ran_vals") #back transform
blups_BG #get estimate and std.error
write.table(blups_BG, file="BG_blups_sellman_28Oct.csv",sep=',')

#plot results BLUPs
sellman$total_biom<-sellman$AG_biom+sellman$BG_biom
total_b.mod0<-lmer(log(total_biom)~Salinity*Mbiome.trt+Cohort+initial_wt_grams+(1|Genotype), 
                   data=sellman,na.action=na.omit)
total_b.mod1<-lmer(log(total_biom)~Salinity*Mbiome.trt+Cohort+initial_wt_grams+(Salinity|Genotype), 
                   data=sellman,na.action=na.omit)
total_b.mod2<-lmer(log(total_biom)~Salinity*Mbiome.trt+Cohort+initial_wt_grams+(Mbiome.trt|Salinity:Genotype), 
                   data=sellman,na.action=na.omit)
anova(total_b.mod1,total_b.mod2,total_b.mod0)

#model 2 is significant
#get BLUPS
blups_totalbiom<- broom.mixed::tidy(total_b.mod2, effects="ran_vals") #back transform
blups_totalbiom #get estimate and std.error
write.table(blups_totalbiom, file="total_biom_blups_sellman_new_28Oct.csv",sep=',')

sellman$rtst<-sellman$BG_biom/sellman$AG_biom

rtst.mod1<-lmer(log(rtst)~Salinity*Mbiome.trt+Cohort+initial_wt_grams+(1|Genotype), 
                data=sellman,na.action=na.omit)
rtst.mod2<-lmer(log(rtst)~Salinity*Mbiome.trt+Cohort+initial_wt_grams+(Salinity|Genotype), 
                data=sellman,na.action=na.omit)
rtst.mod3<-lmer(log(rtst)~Salinity*Mbiome.trt+Cohort+initial_wt_grams+(Mbiome.trt|Salinity:Genotype), 
                data=sellman,na.action=na.omit)
anova(rtst.mod1,rtst.mod2,rtst.mod3)

#height
height.mod00<-lmer(log(final_aveheight)~Salinity*Mbiome.trt+Cohort+initial_wt_grams+(1|Genotype), 
                   data=sellman,na.action=na.omit)
height.mod01<-lmer(log(final_aveheight)~Salinity*Mbiome.trt+Cohort+initial_wt_grams+(Salinity|Genotype), 
                   data=sellman,na.action=na.omit)
height.mod02<-lmer(log(final_aveheight)~Salinity*Mbiome.trt+Cohort+initial_wt_grams+(Mbiome.trt|Salinity:Genotype), 
                   data=sellman,na.action=na.omit)
anova(height.mod00,height.mod01,height.mod02)

sellman$plantsize<-sellman$stemdiam_ave*sellman$final_aveheight
size.mod0<-lmer(log(plantsize)~Salinity*Mbiome.trt+Cohort+initial_wt_grams+(1|Genotype), 
                data=sellman,na.action=na.omit)
size.mod1<-lmer(log(plantsize)~Salinity*Mbiome.trt+Cohort+initial_wt_grams+(Salinity|Genotype), 
                data=sellman,na.action=na.omit)
size.mod2<-lmer(log(plantsize)~Salinity*Mbiome.trt+Cohort+initial_wt_grams+(Mbiome.trt|Salinity:Genotype), 
                data=sellman,na.action=na.omit)
anova(size.mod0,size.mod1,size.mod2)
#model 2 is significant
#get BLUPS
blups_size.sm<- broom.mixed::tidy(size.mod2, effects="ran_vals") #back transform
blups_size.sm #get estimate and std.error
write.table(blups_size.sm, file="size_blups_sellman_28Oct.csv",sep=',')
diam.mod0<-lmer(log(stemdiam_ave)~Salinity*Mbiome.trt+Cohort+initial_wt_grams+(1|Genotype), 
                data=sellman,na.action=na.omit)
diam.mod1<-lmer(log(stemdiam_ave)~Salinity*Mbiome.trt+Cohort+initial_wt_grams+(Salinity|Genotype), 
                data=sellman,na.action=na.omit)
diam.mod1<-lmer(log(stemdiam_ave)~Salinity*Mbiome.trt+Cohort+initial_wt_grams+(Mbiome.trt|Salinity:Genotype), 
                data=sellman,na.action=na.omit)
anova(diam.mod0,diam.mod1)
#SD
#area of each pot = 36 cm2
sellman$stemden<-sellman$stem_num/36
stemden.mod0<-lmer(log(stemden)~Salinity+Mbiome.trt+Cohort+initial_wt_grams+(1|Genotype), 
                   data=sellman,na.action=na.omit)
stemden.mod1<-lmer(log(stemden)~Salinity*Mbiome.trt+Cohort+initial_wt_grams+(Salinity|Genotype), 
                   data=sellman,na.action=na.omit)
stemden.mod2<-lmer(log(stemden)~Salinity*Mbiome.trt+Cohort+initial_wt_grams+(Mbiome.trt|Salinity:Genotype), 
                   data=sellman,na.action=na.omit)
anova(stemden.mod0,stemden.mod1,stemden.mod2)

#**************#
# CORN
trait.d<-read.csv(file="All_data_plant_measurement_less_samp_corrected.csv")
corn.d<-trait.d[trait.d$Site=="Corn",]
#replace 0 with NAs
corn.d[corn.d==0]<-NA
corn.d$Genotype<-as.factor(corn.d$Genotype)

#GLB
lm.mod.grn1<-lmer(log(greenlf_biom)~Salinity*Mbiome.trt*Cohort+initial_wt_grams+(1|Genotype), 
                  data=corn.d,na.action=na.omit)
#vary slope then compare to model with varying intercept
lm.mod.grn2<-lmer(log(greenlf_biom)~Salinity*Mbiome.trt+Cohort+initial_wt_grams+(Salinity|Genotype), 
                  data=corn.d,na.action=na.omit)
lm.mod.grn3<-lmer(log(greenlf_biom)~Salinity*Mbiome.trt*Cohort+initial_wt_grams+(Mbiome.trt|Salinity:Genotype), 
                  data=corn.d,na.action=na.omit)
anova(lm.mod.grn1,lm.mod.grn3,lm.mod.grn2) # test if genotype 
#AG
agb.mod0<-lmer(log(AG_biom)~Salinity*Mbiome.trt+Cohort+initial_wt_grams+(1|Genotype), 
               data=corn.d,na.action=na.omit)
agb.mod1<-lmer(log(AG_biom)~Salinity*Mbiome.trt+Cohort+initial_wt_grams+(Salinity|Genotype), 
               data=corn.d,na.action=na.omit)
agb.mod2<-lmer(log(AG_biom)~Salinity*Mbiome.trt*Cohort+initial_wt_grams+(Mbiome.trt|Salinity:Genotype), 
               data=corn.d,na.action=na.omit)
anova(agb.mod0,agb.mod1,agb.mod2)

#BG
bgb.mod0<-lmer(log(BG_biom)~Salinity*Mbiome.trt+Cohort+initial_wt_grams+(1|Genotype), 
               data=corn.d,na.action=na.omit)
bgb.mod1<-lmer(log(BG_biom)~Salinity*Mbiome.trt+Cohort+initial_wt_grams+(Salinity|Genotype), 
               data=corn.d,na.action=na.omit)
bgb.mod2<-lmer(log(BG_biom)~Salinity*Mbiome.trt+Cohort+initial_wt_grams+(Mbiome.trt|Salinity:Genotype), 
               data=corn.d,na.action=na.omit)
anova(bgb.mod0,bgb.mod1,bgb.mod2)
#model 2 is significant but deviation is zero so no BLUPS
#Total Biomass
corn.d$total_biom<-corn.d$AG_biom+corn.d$BG_biom

total_b.mod0<-lmer(log(total_biom)~Salinity*Mbiome.trt+Cohort+initial_wt_grams+(1|Genotype), 
                   data=corn.d,na.action=na.omit)
total_b.mod1<-lmer(log(total_biom)~Salinity*Mbiome.trt+Cohort+initial_wt_grams+(Salinity|Genotype), 
                   data=corn.d,na.action=na.omit)
total_b.mod2<-lmer(log(total_biom)~Salinity*Mbiome.trt+Cohort+initial_wt_grams+(Mbiome.trt|Salinity:Genotype), 
                   data=corn.d,na.action=na.omit)
anova(total_b.mod0,total_b.mod1,total_b.mod2)

#root-to-shoot ratio
corn.d$rtst<-corn.d$BG_biom/corn.d$AG_biom

rtst.mod1<-lmer(log(rtst)~Salinity*Mbiome.trt+Cohort+initial_wt_grams+(1|Genotype), 
                data=corn.d,na.action=na.omit)
rtst.mod2<-lmer(log(rtst)~Salinity*Mbiome.trt+Cohort+initial_wt_grams+(Salinity|Genotype), 
                data=corn.d,na.action=na.omit)
rtst.mod3<-lmer(log(rtst)~Salinity*Mbiome.trt+Cohort+initial_wt_grams+(Mbiome.trt|Salinity:Genotype), 
                data=corn.d,na.action=na.omit)
anova(rtst.mod1,rtst.mod2,rtst.mod3)
#Plant height
height.mod00<-lmer(log(final_aveheight)~Salinity*Mbiome.trt+Cohort+initial_wt_grams+(1|Genotype), 
                   data=corn.d,na.action=na.omit)
height.mod01<-lmer(log(final_aveheight)~Salinity*Mbiome.trt+Cohort+initial_wt_grams+(Salinity|Genotype), 
                   data=corn.d,na.action=na.omit)
height.mod02<-lmer(log(final_aveheight)~Salinity*Mbiome.trt+Cohort+initial_wt_grams+(Mbiome.trt|Salinity:Genotype), 
                   data=corn.d,na.action=na.omit)
anova(height.mod00,height.mod01,height.mod02)
#get BLUPS
blups_height.corn<- broom.mixed::tidy(height.mod02, effects="ran_vals") #back transform
blups_height.corn #get estimate and std.error
write.table(blups_height.corn, file="height_blups_corn_28Oct_2.csv",sep=',', row.names=FALSE)
#PS
corn.d$plantsize<-corn.d$stemdiam_ave*corn.d$final_aveheight
size.mod0<-lmer(log(plantsize)~Salinity*Mbiome.trt+Cohort+initial_wt_grams+(1|Genotype), 
                data=corn.d,na.action=na.omit)
size.mod1<-lmer(log(plantsize)~Salinity*Mbiome.trt+Cohort+initial_wt_grams+(Salinity|Genotype), 
                data=corn.d,na.action=na.omit)
size.mod2<-lmer(log(plantsize)~Salinity*Mbiome.trt+Cohort+initial_wt_grams+(Mbiome.trt|Salinity:Genotype), 
                data=corn.d,na.action=na.omit)
anova(size.mod0,size.mod1,size.mod2)
#get BLUPS
blups_size.corn<- broom.mixed::tidy(size.mod2, effects="ran_vals") #back transform
blups_size.corn #get estimate and std.error
write.table(blups_size.corn, file="plantsize_blups_corn_28Oct.csv",sep=',', row.names=FALSE)

#Stem diameter
diam.mod0<-lmer(log(stemdiam_ave)~Salinity*Mbiome.trt+Cohort+initial_wt_grams+(1|Genotype), 
                data=corn.d,na.action=na.omit)
diam.mod1<-lmer(log(stemdiam_ave)~Salinity*Mbiome.trt+Cohort+initial_wt_grams+(Salinity|Genotype), 
                data=corn.d,na.action=na.omit)
diam.mod2<-lmer(log(stemdiam_ave)~Salinity*Mbiome.trt+Cohort+initial_wt_grams+(Mbiome.trt|Salinity:Genotype), 
                data=corn.d,na.action=na.omit)
anova(diam.mod0,diam.mod1,diam.mod2)
#get BLUPS
blups_diam.corn<- broom.mixed::tidy(diam.mod2, effects="ran_vals") #back transform
blups_diam.corn #get estimate and std.error
write.table(blups_diam.corn, file="diam_blups_corn_28Oct.csv",sep=',', row.names=FALSE)
#Stem number 
#transform
corn.d$sqstem_num<-sqrt(corn.d$stem_num)
stem.num.mod1<-lmer(sqstem_num~Salinity*Mbiome.trt*Cohort+initial_wt_grams+(1|Genotype), 
                    data=corn.d,na.action=na.omit)
stem.num.mod2<-lmer(sqstem_num~Salinity*Mbiome.trt+Cohort+initial_wt_grams+(Salinity|Genotype), 
                    data=corn.d,na.action=na.omit)
stem.num.mod3<-lmer(sqstem_num~Salinity*Mbiome.trt+Cohort+initial_wt_grams+(Mbiome.trt|Salinity:Genotype), 
                    data=corn.d,na.action=na.omit)
anova(stem.num.mod1,stem.num.mod2,stem.num.mod3)
#SD
#area of each pot = 36 cm2
corn.d$stemden<-corn.d$stem_num/36
#no singularity
stemden.mod0<-lmer(log(stemden)~Salinity*Mbiome.trt+Cohort+initial_wt_grams+(1|Genotype), 
                   data=corn.d,na.action=na.omit)
stemden.mod1<-lmer(log(stemden)~Salinity*Mbiome.trt+Cohort+initial_wt_grams+(Salinity|Genotype), 
                   data=corn.d,na.action=na.omit)
stemden.mod2<-lmer(log(stemden)~Salinity*Mbiome.trt+Cohort+initial_wt_grams+(Mbiome.trt|Salinity:Genotype), 
                   data=corn.d,na.action=na.omit)
anova(stemden.mod0,stemden.mod1,stemden.mod2)

#random likelihood test for genotype at site level
#SELLMAN, do the same for corn and across all dataset
#=================================================================#
trait.d<-read.csv(file="All_data_plant_measurement_less_samp_corrected.csv")
sellman<-trait.d[trait.d$Site=="Sellman",]
sellman[sellman==0]<-NA
sellman$Genotype<-as.factor(sellman$Genotype)
options(contrasts = c("contr.sum","contr.poly"))

#GLB
lm.mod.grn2<-lmer(log(greenlf_biom)~Salinity*Mbiome.trt*Cohort+initial_wt_grams+(1|Genotype), 
                  data=sellman,na.action=na.omit)
lm.mod.grn3<-lm(log(greenlf_biom)~Salinity*Mbiome.trt*Cohort+initial_wt_grams, 
                data=sellman,na.action=na.omit)
anova(lm.mod.grn2,lm.mod.grn3)
#AG
agb.mod0<-lmer(log(AG_biom)~Salinity*Mbiome.trt*Cohort+initial_wt_grams+(1|Genotype), 
               data=sellman,na.action=na.omit)
agb.mod1<-lm(log(AG_biom)~Salinity*Mbiome.trt*Cohort+initial_wt_grams, 
             data=sellman,na.action=na.omit)
anova(agb.mod0, agb.mod1)

#BG
bgb.mod0<-lmer(log(BG_biom)~Salinity*Mbiome.trt*Cohort+initial_wt_grams+(1|Genotype), 
               data=sellman,na.action=na.omit)
bgb.mod1<-lm(log(BG_biom)~Salinity*Mbiome.trt*Cohort+initial_wt_grams, 
             data=sellman,na.action=na.omit)
anova(bgb.mod0, bgb.mod1)
#Total Biomass
sellman$total_biom<-sellman$AG_biom+sellman$BG_biom
total_b.mod0<-lmer(log(total_biom)~Salinity*Mbiome.trt*Cohort+initial_wt_grams+(1|Genotype), 
                   data=sellman,na.action=na.omit)
total_b.mod1<-lm(log(total_biom)~Salinity*Mbiome.trt*Cohort+initial_wt_grams, 
                 data=sellman,na.action=na.omit)
anova(total_b.mod0, total_b.mod1)
#root-to-shoot ratio
sellman$rtst<-sellman$BG_biom/sellman$AG_biom
rtst.mod1<-lmer(log(rtst)~Salinity*Mbiome.trt*Cohort+initial_wt_grams+(1|Genotype), 
                data=sellman,na.action=na.omit)
rtst.mod2<-lm(log(rtst)~Salinity*Mbiome.trt*Cohort+initial_wt_grams, 
              data=sellman,na.action=na.omit)
anova(rtst.mod1, rtst.mod2)

#Plant height
height.mod00<-lmer(log(final_aveheight)~Salinity*Mbiome.trt*Cohort+initial_wt_grams+(1|Genotype), 
                   data=sellman,na.action=na.omit)
height.mod1<-lm(log(final_aveheight)~Salinity*Mbiome.trt*Cohort+initial_wt_grams, 
                data=sellman,na.action=na.omit)
anova(height.mod00, height.mod1)
#Plant size
sellman$plantsize<-sellman$stemdiam_ave*sellman$final_aveheight
size.mod2<-lmer(log(plantsize)~Salinity*Mbiome.trt*Cohort+initial_wt_grams+(1|Genotype), 
                data=sellman,na.action=na.omit)
size.mod3<-lm(log(plantsize)~Salinity*Mbiome.trt*Cohort+initial_wt_grams, 
              data=sellman,na.action=na.omit)
anova(size.mod2, size.mod3)

#Stem diameter
diam.mod0<-lmer(log(stemdiam_ave)~Salinity*Mbiome.trt*Cohort+initial_wt_grams+(1|Genotype), 
                data=sellman,na.action=na.omit)
diam.mod1<-lm(log(stemdiam_ave)~Salinity*Mbiome.trt*Cohort+initial_wt_grams, 
              data=sellman,na.action=na.omit)
anova(diam.mod0, diam.mod1)
#Stem number 
#transform
sellman$sqstem_num<-sqrt(sellman$stem_num)
# singularity
stem.num.mod1<-lmer(sqstem_num~Salinity*Mbiome.trt*Cohort+initial_wt_grams+(1|Genotype), 
                    data=sellman,na.action=na.omit)
stem.num.mod2<-lm(sqstem_num~Salinity*Mbiome.trt*Cohort+initial_wt_grams, 
                  data=sellman,na.action=na.omit)
anova(stem.num.mod1, stem.num.mod2)
#Stem density
#area of each pot = 36 cm2
sellman$stemden<-sellman$stem_num/36
# singularity
stemden.mod1<-lmer(log(stemden)~Salinity*Mbiome.trt*Cohort+initial_wt_grams+(1|Genotype), 
                   data=sellman,na.action=na.omit)
stemden.mod2<-lm(log(stemden)~Salinity*Mbiome.trt*Cohort+initial_wt_grams, 
                 data=sellman,na.action=na.omit)
anova(stemden.mod1, stemden.mod2)

#microbe effect test 
#change into wide format
df.p<-read.csv(file="All_data_plant_measurement_less_samp_corrected.csv")
df.p$id<-paste(df.p$Genotype,df.p$Clone,sep="-")
df.p$Sal_ino<-paste(df.p$Salinity,df.p$Mbiome.trt,sep="-")
write.table(df.p, file="plant_trait_long.csv", sep=',', row.names=FALSE)
data_wide <- dcast(df.p, id+ Site+Salinity ~ Mbiome.trt, value.var="final_aveheight")
#subtract for microbiome effect
data_wide$mic_eff<-data_wide$Live-data_wide$Sterile
write.table(data_wide, file="plant_trait_wide_height.csv", sep=',', row.names=FALSE)
#plant size
df<-read.csv(file="plant_trait_long.csv")
df$plantsize<-df$stem_num*df$final_aveheight
df.ps<- dcast(df, id+ Site+Salinity+Cohort ~ Mbiome.trt, value.var="plantsize")
#subtract for microbiome effect
df.ps$mic_eff<-df.ps$Live-df.ps$Sterile
write.table(df.ps, file="plant_trait_wide_size.csv", sep=',', row.names=FALSE)
#stem number
df<-read.csv(file="plant_trait_long.csv")
df.sd<- dcast(df, id+ Site+Salinity+Cohort ~ Mbiome.trt, value.var="stem_num")
#subtract for microbiome effect
df.sd$mic_eff<-df.sd$Live-df.sd$Sterile
write.table(df.sd, file="plant_trait_wide_stem_num.csv", sep=',', row.names=FALSE)
#stem_density
df<-read.csv(file="plant_trait_long.csv")
df$stemden<-df$stem_num/36
df.sd<- dcast(df, id+ Site+Salinity+Cohort ~ Mbiome.trt, value.var="stemden")
#subtract for microbiome effect
df.sd$mic_eff<-df.sd$Live-df.sd$Sterile
write.table(df.sd, file="plant_trait_wide_density.csv", sep=',', row.names=FALSE)

#manually parse out the genotype
#GET emmeans
#height
df.all<-read.csv(file="plant_trait_wide_height.csv")
df.all$Genotype<-as.factor(df.all$Genotype)
options(contrasts = c("contr.sum","contr.poly"))
ph.mod<-lmer(mic_eff~Salinity*Site*Cohort+(1|Genotype), 
             data=df.all,na.action=na.omit)
#test for main effects
anova(ph.mod)
#Post-hoc tests:
all.ht.emm<-emmeans(ph.mod, pairwise~Salinity:Cohort:Site)
all.ht.emm
summary(all.ht.emm, infer = TRUE, null = log(35), type = "response")
all.ht.cld<-multcomp::cld(all.ht.emm$emmean,alpha=0.05, Letters=letters, adjust="tukey")
all.ht.cld
#SN
df.sn<-read.csv(file="plant_trait_wide_stem_num.csv")
df.sn$Genotype<-as.factor(df.sn$Genotype)
options(contrasts = c("contr.sum","contr.poly"))
sn.mod<-lmer(mic_eff~Salinity*Site*Cohort+(1|Genotype), 
             data=df.sn,na.action=na.omit)
#test for main effects
anova(sn.mod)
#Post-hoc tests: get emmeans/lsmeans to plot
sn.emmeans<-eameans(sn.mod, pairwise~Salinity:Cohort:Site)
sn.emmeans
#SD density
df.sd<-read.csv(file="plant_trait_wide_density.csv")
df.sd$Genotype<-as.factor(df.sd$Genotype)
options(contrasts = c("contr.sum","contr.poly"))
den.mod<-lmer(mic_eff~Salinity*Site*Cohort+(1|Genotype), 
              data=df.sd,na.action=na.omit)
#test for main effects
anova(den.mod)
#Post-hoc tests: 
den.emmeans<-emmeans(den.mod, pairwise~Salinity:Cohort:Site)
den.emmeans


#+++++++++++#
#PLSR, do for both fungi and bacteria
fundiv<-read.csv(file="FinalData_scam_root_fungi.csv")
#replace 0 with NAs
fundiv[fundiv==0]<-NA
fundiv$Genotype<-as.factor(fundiv$Genotype)

#perform plsr
fundiv$total_biom<-fundiv$AG_biom+fundiv$BG_biom
fundiv$stemden<-fundiv$stem_num/36
fundiv$RS<-fundiv$BG_biom/fundiv$AG_biom
fundiv$PS<-fundiv$stem_num*fundiv$final_aveheight
fun.div.pls<-plsr(enspie~scale(final_aveheight)+scale(AG_biom) +scale(BG_biom)+
                    scale(greenlf_biom)+scale(RS)+scale(PS)+
                    scale(stem_num)+scale(stemdiam_ave)+scale(stemden),
                  data=fundiv, validation="CV")
# Find the number of dimensions with lowest cross validation error
cv<-RMSEP(fun.div.pls)
best.dims<-which.min(cv$val[estimate = "adjCV", , ]) - 1
best.dims
# Rerun the model
new.fundiv.plsr<-plsr(enspie~scale(final_aveheight)+scale(AG_biom) +scale(BG_biom)+
                        scale(greenlf_biom)+scale(RS)+scale(PS)+
                        scale(stem_num)+scale(stemdiam_ave)+scale(stemden), 
                      data=fundiv, ncomp = 2, validation="CV", scale=FALSE,jackknife = TRUE)
summary(new.fundiv.plsr)
jack.test(new.fundiv.plsr,ncomp=1:2, use.mean = TRUE)
#extract coefficient
coefficients<-coef(new.fundiv.plsr, ncomp=1:2)
coefficients
#extract explained variance attributed to components
compnames(new.fundiv.plsr, comps = 1:2, explvar = TRUE)
#comps = 1:2, 
factor=(row.names(coefficients))
coef<-data.frame(factor)
coef
# change name
levels(coef$factor)[levels(coef$factor) =="scale(AG_biom)"]<-"AG"
levels(coef$factor)[levels(coef$factor) =="scale(BG_biom)"]<-"BG"
levels(coef$factor)[levels(coef$factor) =="scale(greenlf_biom)"]<-"GB"
#levels(coef$factor)[levels(coef$factor) =="scale(total_biom)"]<-"Total"
levels(coef$factor)[levels(coef$factor) =="scale(final_aveheight)"]<-"PH"
levels(coef$factor)[levels(coef$factor) =="scale(stem_num)"]<-"SN"
levels(coef$factor)[levels(coef$factor) =="scale(stemdiam_ave)"]<-"SDi"
levels(coef$factor)[levels(coef$factor) =="scale(stemden)"]<-"SD"
levels(coef$factor)[levels(coef$factor) =="scale(RS)"]<-"R:S"
levels(coef$factor)[levels(coef$factor) =="scale(PS)"]<-"PS"
coef
#plot correlation loadings 
tiff("fungi_plsr_corr.tiff", width = 110, height = 110, units = 'mm', res = 300)
plot(new.fundiv.plsr, plottype = "correlation", 
     ploty=FALSE, labels=coef$factor, cex=0.85)
par(new=TRUE)
plot(new.fundiv.plsr, plottype = "correlation", 
     ploty=TRUE, plotx=FALSE, pch=17, col="darkred")

library(vegan)
library(ggplot2)
library(cowplot)

#+++++++++++++++++++++++++++++++++++++++#
#(+) sol-microbe treated plants only endosphere bacteria 
#do the same for funig
endo.bac<-read.csv(file="FinalData_scam_root_bacteria.csv")
live.bac.endo<-endo.bac[endo.bac$Mbiome.trt=="Live"  ,]
#Setting up contrasts:
options(contrasts=c("contr.sum","contr.poly"))
#test lme/lm
bact.endo.live.mod1<-lmer(log(enspie)~Salinity*Cohort*Site+(1|Genotype), data=live.bac.endo, na.action=na.omit)
anova(bact.endo.live.mod1, dff="Satterthwaite")
summary(bact.endo.live.mod1)
#random effects
bact.endo.live.mod1<-lmer(log(enspie)~Salinity*Cohort*Site+(1|Genotype), data=live.bac.endo, na.action=na.omit)
bact.endo.live.mod2 <-lm(log(enspie)~Salinity*Cohort*Site, data=live.bac.endo, na.action=na.omit)
anova(bact.endo.live.mod1,bact.endo.live.mod2)

#(-) soil-microbe plants
endo.bac<-read.csv(file="FinalData_scam_root_bacteria.csv")
sterile.endo.bac<-endo.bac[endo.bac$Mbiome.trt=="Sterile" ,]
#Setting up contrasts:
options(contrasts=c("contr.sum","contr.poly"))
#test 
bact.endo.sterile.mod1<-lmer(log(enspie)~Salinity*Cohort*Site+(1|Genotype), data=sterile.endo.bac, na.action=na.omit)
anova(bact.endo.sterile.mod1, dff="Satterthwaite")
summary(bact.endo.sterile.mod1)
Anova(bact.endo.sterile.mod1, type=3)
#random effects
bact.endo.sterile.mod1<-lmer(log(enspie)~Salinity*Cohort*Site+(1|Genotype), data=sterile.endo.bac, na.action=na.omit)
bact.endo.sterile.mod2 <-lm(log(enspie)~Salinity*Cohort*Site, data=sterile.endo.bac, na.action=na.omit)
anova(bact.endo.sterile.mod1,bact.endo.sterile.mod2)

#Root bacteria community composition
#do the same for fungi
#within vs between site; across all data
endo.bac<-read.csv(file="FinalData_scam_root_bacteria.csv")
#make as dummy variable
# make dummy variable for Salinity
endo.bac$salinity<-NA
endo.bac$salinity[endo.bac$Salinity=="Low_0ppt"] <-0
endo.bac$salinity[endo.bac$Salinity=="High_15ppt"] <-1
with(endo.bac, table(salinity))
# make dummy variable for microbiome treatment
endo.bac$microbiome<-NA
endo.bac$microbiome[endo.bac$Mbiome.trt=="Sterile"] <-0
endo.bac$microbiome[endo.bac$Mbiome.trt=="Live"] <-1
with(endo.bac, table(microbiome))
# make dummy variable for microbiome treatment
endo.bac$cohort<-NA
endo.bac$cohort[endo.bac$Cohort=="Ancestral"] <-0
endo.bac$cohort[endo.bac$Cohort=="Modern"] <-1
with(endo.bac, table(cohort))
# make dummy variable for microbiome treatment
endo.bac$site<-NA
endo.bac$site[endo.bac$Site=="Corn"] <-0
endo.bac$site[endo.bac$Site=="Sellman"] <-1
with(endo.bac, table(site))

#make species matrix
sp.cols.rt <- grep("otu_", names(endo.bac))

# Make Bray Curtis dissimilarity matrix 
# Square root transformation, Wisconsin double standardization,#this emphasizes the environmental variables
scam.root.mat<-vegdist((endo.bac[,sp.cols.rt]), method="bray", binary=FALSE, metaMDS=TRUE, sqrt.dist=TRUE)

dbRDA_all.root.1 <- capscale(scam.root.mat ~ site+salinity+cohort+microbiome, 
                             data=endo.bac, na.option=na.omit)
#Null Model
#generate one model with NOTHING to explain the braycurtis dm matrix
dbRDA_all.root.0 <- capscale(scam.root.mat~1,data=endo.bac)

#use forward selection to choose which elements of the full model explain a significant amount of variaiton in the unifrac dm by comparing it to the null
dbRDA_all.root <- ordistep(dbRDA_all.root.0,scope=formula(dbRDA_all.root.1),
                           direction="forward",Pin=.1,trace=T,pstep=500000)
dbRDA_all.root
#test significance of analysis
anova(dbRDA_all.root)
#test 
anova(dbRDA_all.root, by="margin")
anova(dbRDA_all.root, by="terms", permu=200) # same above,test for sign. environ. variables
# test axes for significance
anova(dbRDA_all.root, by="axis", perm.max=500) 
#barplot inertia for each PCo axes
screeplot(dbRDA_all.root)
#retrieve coefficient
coef(dbRDA_all.root)
#get residual from the model
resid(dbRDA_all.root)


#PLOT
#any rid of NA values?
anyNA(endo.bac)
#make model sumary
B <- summary(dbRDA_all.root)
#plot
A.1 <- scores(dbRDA_all.root)
A.2 <- A.1$sites
A.3 <- cbind(A.2, endo.bac)
#scores for arows
A.4 <- data.frame(scores(dbRDA_all.root, display = "bp"))
#subset A4 for labeling
A.4 <- A.4[sort(rownames(A.4)),]
A4.sub1 <- A.4[2,]
A4.sub2 <- A.4[1,]
A4.sub3 <- A.4[3,]
A4.sub4 <- A.4[4,]
#make plot
p <- ggplot(data = A.3, aes(x = CAP1, y = CAP2))
p.dbrda.root.bac <- p +
  geom_point(data = A.3, alpha = 2/5, size=4,
             aes(shape = Site, color = Mbiome.trt, fill=Mbiome.trt, stroke = 1)) +
  theme_bw() + scale_color_manual(values=c("black", "red"))+
  scale_fill_manual(values=c("black", "red"))+
  xlab(label = paste("CAP1 (", round(B$concont$importance[2,1]*100, digits = 1), "%)", sep="")) +
  ylab(label = paste("CAP2 (", round(B$concont$importance[2,2]*100, digits = 1), "%)", sep="")) +
  scale_shape_manual(name = "Site", labels = c("Corn", "Sellman"), 
                     values = c(1,17))
p.dbrda.root.bac
#ggsave(filename="Fig_dbrda_root_bacteria.tiff", plot=p.dbrda.root.bac, dpi=400, units=c("mm"), width=120, height=100)

#(LIVE)+) soil-microbe treatment only
#within vs between site; across all data
endo.bac<-read.csv(file="All_root_bacteria_data_rarefied.csv")
live.bac<-endo.bac[endo.bac$Mbiome.trt=="Live"  ,]
#make as dummy variable
# make dummy variable for Salinity
live.bac$salinity<-NA
live.bac$salinity[live.bac$Salinity=="Low_0ppt"] <-0
live.bac$salinity[live.bac$Salinity=="High_15ppt"] <-1
with(live.bac, table(salinity))
# make dummy variable for microbiome treatment
live.bac$microbiome<-NA
live.bac$microbiome[live.bac$Mbiome.trt=="Sterile"] <-0
live.bac$microbiome[live.bac$Mbiome.trt=="Live"] <-1
with(live.bac, table(microbiome))
# make dummy variable for microbiome treatment
live.bac$cohort<-NA
live.bac$cohort[live.bac$Cohort=="Ancestral"] <-0
live.bac$cohort[live.bac$Cohort=="Modern"] <-1
with(live.bac, table(cohort))
# make dummy variable for microbiome treatment
live.bac$site<-NA
live.bac$site[live.bac$Site=="Corn"] <-0
live.bac$site[live.bac$Site=="Sellman"] <-1
with(live.bac, table(site))

#make species matrix
sp.cols.rt <- grep("otu_", names(live.bac))
# Make Bray Curtis dissimilarity matrix 
# Square root transformation, Wisconsin double standardization,#this emphasizes the environmental variables
scam.root.mat<-vegdist((live.bac[,sp.cols.rt]), method="bray", binary=FALSE, metaMDS=TRUE, sqrt.dist=TRUE)

dbRDA_all.root.1 <- capscale(scam.root.mat ~ site+salinity+cohort+microbiome, 
                             data=live.bac, na.option=na.omit)
#Null Model
#generate one model with NOTHING to explain the braycurtis dm matrix
dbRDA_all.root.0 <- capscale(scam.root.mat~1,data=live.bac)
#use forward selection to choose which elements of the full model explain a significant amount of variaiton in the unifrac dm by comparing it to the null

dbRDA_all.root <- ordistep(dbRDA_all.root.0,scope=formula(dbRDA_all.root.1),
                           direction="forward",Pin=.1,trace=T,pstep=500000)
dbRDA_all.root
#test significance of analysis
anova(dbRDA_all.root)
#test 
anova(dbRDA_all.root, by="margin")
anova(dbRDA_all.root, by="terms", permu=200) # same above,test for sign. environ. variables
# test axes for significance
anova(dbRDA_all.root, by="axis", perm.max=500) 
#barplot inertia for each PCo axes
screeplot(dbRDA_all.root)
#retrieve coefficient
coef(dbRDA_all.root)
#get residual from the model
resid(dbRDA_all.root)
#PLOT
anyNA(live.bac)
#make model sumary
B <- summary(dbRDA_all.root)
#plot
A.1 <- scores(dbRDA_all.root)
A.2 <- A.1$sites
A.3 <- cbind(A.2, live.bac)
#scores for arows
A.4 <- data.frame(scores(dbRDA_all.root, display = "bp"))
#subset A4 for labeling
A.4 <- A.4[sort(rownames(A.4)),]
A4.sub1 <- A.4[2,]
A4.sub2 <- A.4[1,]
#make plot
p <- ggplot(data = A.3, aes(x = CAP1, y = CAP2))

p.dbrda.root.bac <- p +
  geom_point(data = A.3, alpha = 2/5, size=4,
             aes(shape = Site, color = Salinity, fill=Salinity, stroke = 1)) +
  theme_bw() + scale_color_manual(values=c("black", "red"))+
  scale_fill_manual(values=c("black", "red"))+
  xlab(label = paste("CAP1 (", round(B$concont$importance[2,1]*100, digits = 1), "%)", sep="")) +
  ylab(label = paste("CAP2 (", round(B$concont$importance[2,2]*100, digits = 1), "%)", sep="")) +
  scale_shape_manual(name = "Site", labels = c("Corn", "Sellman"), 
                     values = c(1,17))
p.dbrda.root.bac
ggsave(filename="Fig_dbrda_root_bacteria_live.tiff", plot=p.dbrda.root.bac, dpi=400, units=c("mm"), width=120, height=100)
#-------------------------#
#indicator species_using indicspecies
#-------------------------#
endo.bac<-read.csv(file="FinalData_scam_root_bacteria.csv")
#select columns with the otu and convert to mvabubd df
bac.otu.ind <-endo.bac[,33:7240]
group<-as.factor(interaction(endo.bac$Salinity,endo.bac$Mbiome.trt,endo.bac$Cohort))

#conduct analysis; restrict to 3 groupings
indc.bac<- multipatt(bac.otu.ind, group, max.order=1,
                     func = "r.g", control = how(nperm=9999))
data.frame(unclass(summary(indc.bac)))
write.table(indic.sum, file="root_bacteria_indic_spp.csv", sep=',',quote = FALSE, row.names = F)
#identify otus 
otu<-read.csv(file="otu_id.csv")
taxa<-read.csv(file="bacteria_root_tax_ID.csv")
#keep order of the rows in otu file
merge.d<-merge(otu, taxa, by="otu_id", sort=FALSE)
write.table(merge.d, file="top_10_indicator_spp.csv", sep=',', row.names=FALSE)
endo.fung<-read.csv(file="FinalData_scam_root_fungi.csv")
#FUNGI
endo.fung<-read.csv(file="FinalData_scam_root_fungi.csv")
#select columns with the otu and convert to mvabubd df
fungi.otu.ind <-endo.fung[,33:578]
group<-as.factor(interaction(endo.fung$Salinity,endo.fung$Mbiome.trt,endo.fung$Cohort))

#conduct analysis; restrict to 3 groupings
indc.fungi<- multipatt(fungi.otu.ind, group, max.order=1,
                       func = "r.g", control = how(nperm=9999))
summary(indc.fungi)

#Figure 2
library(reshape2)
library(ggplot2)
library(gridExtra)
library(cowplot)
library(VennDiagram)
library(plyr)

###############################
#Plotting Figure 2 combined
###############################
#================================#
#Fig 2a.height
#No Cohort since not significant
#================================#
#combine in one graph both ecotype, no cohort
trait.d<-read.csv(file="All_data_plant_measurement_less_samp_corrected.csv")
levels(trait.d$Site)[levels(trait.d$Site)=="Corn"]<-"CI"
levels(trait.d$Site)[levels(trait.d$Site)=="Sellman"]<-"SM"
levels(trait.d$Cohort)[levels(trait.d$Cohort)=="Modern"]<-"Descendant"

#graph labels for microbiome treatment
# (+) - live soil microbiome
#(-) sterile

#replace 0 with NAs
trait.d[trait.d==0]<-NA
height.all<-ddply(trait.d, c("Mbiome.trt","Salinity","Site"), summarise, N = length (final_aveheight), mean = mean(final_aveheight, na.rm=TRUE), 
                  sd= sd(final_aveheight, na.rm=TRUE), se = var(final_aveheight/sqrt(N), na.rm=TRUE))

#barplot
height.plot.all<-ggplot(subset(height.all, Mbiome.trt %in% c("Live", "Sterile")), aes(x=as.factor(interaction(Salinity,Mbiome.trt)), y=mean, color=Site, 
                                                                                      fill=Site, group=as.factor(interaction(Site, Mbiome.trt)))) + 
  scale_fill_manual(values=c("gray96", "dodgerblue"))+scale_color_manual(values=c("black", "black"))+
  geom_bar(stat="identity", position=position_dodge(), width=0.8)+ 
  ggtitle("")+ 
  ylim(0,60)+
  ylab("stem height (cm) \n")+xlab("")+
  guides(fill = guide_legend(override.aes = list(linetype = 0)))+
  scale_x_discrete(limits=c("Low_0ppt.Live","Low_0ppt.Sterile","High_15ppt.Live","High_15ppt.Sterile"),
                   labels=c("Low\n(+) soil-\nmicrobe","Low\n(-) soil-\nmicrobe","High\n(+) soil-\nmicrobe","High\n(-) soil-\nmicrobe"))+
  #guides(fill = guide_legend(override.aes = list(linetype = 0, color="black")))+
  theme(axis.text=element_text(size=10), axis.line.y.left=element_line(color="black"),
        axis.line.x.bottom=element_line(color="black"),
        axis.title.y=element_text(size=10),
        legend.position="top",legend.direction="horizontal",
        legend.key.size =unit(5, "mm"),legend.title=element_blank(),legend.text = element_text(size=12), 
        panel.grid.minor=element_blank(),panel.grid.major=element_blank(), 
        panel.background = element_blank(), plot.title=element_blank())+
  geom_errorbar(aes(ymin=mean, ymax=mean+se), position=position_dodge(0.75), width=0.2)
#geom_text(aes(label = .group),position = position_dodge(0.9), vjust = -5)
height.plot.all

#leg1<-get_legend(height.plot.all)
#leg1
#legend.position=c(0.54,0.9)
#height.plot.all1 <- height.plot.all + theme(legend.position = "none")

#================================#
#Fig 2c.SN
#
#================================#
#combine in one graph both ecotype, no cohort
trait.d<-read.csv(file="All_data_plant_measurement_less_samp_corrected.csv")
levels(trait.d$Site)[levels(trait.d$Site)=="Corn"]<-"CI"
levels(trait.d$Site)[levels(trait.d$Site)=="Sellman"]<-"SM"
levels(trait.d$Cohort)[levels(trait.d$Cohort)=="Modern"]<-"Descendant"

#replace 0 with NAs
trait.d[trait.d==0]<-NA
sn.all<-ddply(trait.d, c("Mbiome.trt","Salinity","Site"), summarise, N = length (stem_num), mean = mean(stem_num, na.rm=TRUE), 
              sd= sd(stem_num, na.rm=TRUE), se = var(stem_num/sqrt(N), na.rm=TRUE))
sn.all.sqrt<-ddply(trait.d, c("Mbiome.trt","Salinity","Site"), summarise, N = length (sqrt(stem_num)), mean = mean(sqrt(stem_num), na.rm=TRUE), 
                   sd= sd(sqrt(stem_num), na.rm=TRUE), se = var(sqrt(stem_num)/sqrt(N), na.rm=TRUE))


#barplot
sn.plot.all<-ggplot(subset(sn.all, Mbiome.trt %in% c("Live", "Sterile")), aes(x=as.factor(interaction(Salinity,Mbiome.trt)), y=mean, color=Site, 
                                                                              fill=Site, group=as.factor(interaction(Site, Mbiome.trt)))) + 
  scale_fill_manual(values=c("gray96", "dodgerblue"))+scale_color_manual(values=c("black", "black"))+
  geom_bar(stat="identity", position=position_dodge(), width=0.8)+ 
  ggtitle("")+ 
  ylim(0,7)+
  ylab("stem number \n")+xlab("")+
  guides(fill = guide_legend(override.aes = list(linetype = 0)))+
  scale_x_discrete(limits=c("Low_0ppt.Live","Low_0ppt.Sterile","High_15ppt.Live","High_15ppt.Sterile"),
                   labels=c("Low\n(+) soil-\nmicrobe","Low\n(-) soil-\nmicrobe","High\n(+) soil-\nmicrobe","High\n(-) soil-\nmicrobe"))+
  #guides(fill = guide_legend(override.aes = list(linetype = 0, color="black")))+
  theme(axis.text=element_text(size=10), axis.line.y.left=element_line(color="black"),
        axis.line.x.bottom=element_line(color="black"),
        axis.title.y=element_text(size=10),
        legend.position="top",legend.direction="horizontal",
        legend.key.size =unit(5, "mm"),legend.title=element_blank(),legend.text = element_text(size=12), 
        panel.grid.minor=element_blank(),panel.grid.major=element_blank(), 
        panel.background = element_blank(), plot.title=element_blank())+
  geom_errorbar(aes(ymin=mean, ymax=mean+se), position=position_dodge(0.75), width=0.2)
#geom_text(aes(label = .group),position = position_dodge(0.9), vjust = -5)
sn.plot.all


#================================#
#Fig 2e.SD
#
#================================#
#combine in one graph both ecotype, no cohort
trait.d<-read.csv(file="All_data_plant_measurement_less_samp_corrected.csv")
levels(trait.d$Site)[levels(trait.d$Site)=="Corn"]<-"CI"
levels(trait.d$Site)[levels(trait.d$Site)=="Sellman"]<-"SM"
levels(trait.d$Cohort)[levels(trait.d$Cohort)=="Modern"]<-"Descendant"

#replace 0 with NAs
trait.d[trait.d==0]<-NA
trait.d$stemden<-trait.d$stem_num/36

sd.all<-ddply(trait.d, c("Mbiome.trt","Salinity","Site"), summarise, N = length (stemden), mean = mean(stemden, na.rm=TRUE), 
              sd= sd(stemden, na.rm=TRUE), se = var(stemden/sqrt(N), na.rm=TRUE))
#barplot
sd.plot.all<-ggplot(subset(sd.all, Mbiome.trt %in% c("Live", "Sterile")), aes(x=as.factor(interaction(Salinity,Mbiome.trt)), y=mean, color=Site, 
                                                                              fill=Site, group=as.factor(interaction(Site, Mbiome.trt)))) + 
  scale_fill_manual(values=c("gray96", "dodgerblue"))+scale_color_manual(values=c("black", "black"))+
  geom_bar(stat="identity", position=position_dodge(), width=0.8)+ 
  ggtitle("")+ 
  ylim(0,0.2)+
  ylab(expression(paste("stem density ",(SN/cm),"\n ")))+xlab("\nTreatments")+
  guides(fill = guide_legend(override.aes = list(linetype = 0)))+
  scale_x_discrete(limits=c("Low_0ppt.Live","Low_0ppt.Sterile","High_15ppt.Live","High_15ppt.Sterile"),
                   labels=c("Low\n(+) soil-\nmicrobe","Low\n(-) soil-\nmicrobe","High\n(+) soil-\nmicrobe","High\n(-) soil-\nmicrobe"))+
  #guides(fill = guide_legend(override.aes = list(linetype = 0, color="black")))+
  theme(axis.text=element_text(size=10), axis.line.y.left=element_line(color="black"),
        axis.line.x.bottom=element_line(color="black"),
        axis.title.y=element_text(size=10),
        axis.title.x=element_text(size=12),
        legend.position="top",legend.direction="horizontal",
        legend.key.size =unit(5, "mm"),legend.title=element_blank(),legend.text = element_text(size=12), 
        panel.grid.minor=element_blank(),panel.grid.major=element_blank(), 
        panel.background = element_blank(), plot.title=element_blank())+
  geom_errorbar(aes(ymin=mean, ymax=mean+se), position=position_dodge(0.75), width=0.5, size=0.8)
#geom_text(aes(label = .group),position = position_dodge(0.9), vjust = -5)
sd.plot.all

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#Fig. 1 plot microbe effects height
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
miceff.t<-read.csv(file="plant_trait_microbe_eff_traits_lsmeans.csv")
levels(miceff.t$Provenance)[levels(miceff.t$Provenance)=="Corn"]<-"CI"
levels(miceff.t$Provenance)[levels(miceff.t$Provenance)=="Sellman"]<-"SM"
levels(miceff.t$Cohort)[levels(miceff.t$Cohort)=="Modern"]<-"Descendant"

ph<-miceff.t[miceff.t$Trait=="height",]
ph <- na.omit(ph) 
mic.eff.ph<-ggplot(ph, aes(x = Salinity, y = lsmean, shape= Cohort, color=Cohort)) + 
  geom_line() + #ylim(-25,25)+
  geom_hline(yintercept=0, size=1, lty=2, color='grey50')+
  scale_color_manual(values=c("black","black"))+
  scale_shape_manual(values=c(17,1))+
  scale_x_discrete(limits=c( "Low Salinity","High Salinity"), 
                   labels=c("Low", "High"))+
  geom_errorbar(aes(ymin=lsmean-SE, ymax=lsmean+SE), width=0.01, colour="gray50") +
  geom_point(size=4) +
  #geom_text(aes(label=Cohort), size=4,hjust=-0.5, vjust=0.1)+
  #ggtitle("") +
  facet_grid(~Provenance, scales = "free_x", space = "free_x", shrink=TRUE) +
  #facet_wrap(~Provenance, ncol=2)+
  labs(x="", y= "microbe effect SH") +
  theme(panel.spacing = unit(0.5, "lines"),  #panel.background = element_blank(),
        strip.background = element_blank(), legend.position='none',
        strip.placement = "inside", plot.title = element_text(hjust = 0.2),
        axis.line.y.left=element_line(color="black"),axis.line.x.bottom=element_line(color="black"),
        strip.text = element_text(face="bold", vjust=1.0, size = 14),
        axis.title.x = element_text(vjust=1.0, size = 12),
        axis.title.y = element_text(vjust=1.0, size = 12),
        axis.text.x = element_text(size=12, colour ="black"),
        axis.text.y = element_text(size=14, colour ="black"),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank()) 
mic.eff.ph

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#Fig. 1d plot microbe effects stem number
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
miceff.t<-read.csv(file="plant_trait_microbe_eff_traits_lsmeans.csv")
levels(miceff.t$Provenance)[levels(miceff.t$Provenance)=="Corn"]<-"CI"
levels(miceff.t$Provenance)[levels(miceff.t$Provenance)=="Sellman"]<-"SM"
levels(miceff.t$Cohort)[levels(miceff.t$Cohort)=="Modern"]<-"Descendant"

sn<-miceff.t[miceff.t$Trait=="stem_num",]
sn <- na.omit(sn) 
mic.eff.sn<-ggplot(sn, aes(x = Salinity, y = lsmean, shape= Cohort, color=Cohort)) + 
  geom_line() + #ylim(-25,25)+
  geom_hline(yintercept=0, size=1, lty=2, color='grey50')+
  scale_color_manual(values=c("black","black"))+
  scale_shape_manual(values=c(17,1))+
  scale_x_discrete(limits=c( "Low Salinity","High Salinity"), 
                   labels=c("Low", "High"))+
  geom_errorbar(aes(ymin=lsmean-SE, ymax=lsmean+SE), width=0.01, colour="gray50") +
  geom_point(size=4) +
  #geom_text(aes(label=Cohort), size=4,hjust=-0.5, vjust=0.1)+
  #ggtitle("") +
  facet_grid(~Provenance, scales = "free_x", space = "free_x", shrink=TRUE) +
  #facet_wrap(~Provenance, ncol=2)+
  labs(x="", y= "microbe effect SN\n") +
  theme(panel.spacing = unit(0.5, "lines"),  #panel.background = element_blank(),
        strip.background = element_blank(), legend.position='none',
        strip.placement = "inside", plot.title = element_text(hjust = 0.2),
        axis.line.y.left=element_line(color="black"),axis.line.x.bottom=element_line(color="black"),
        strip.text = element_text(face="bold", vjust=1.0, size = 14),
        axis.title.x = element_text(vjust=1.0, size = 12),
        axis.title.y = element_text(vjust=1.0, size = 12),
        axis.text.x = element_text(size=12, colour ="black"),
        axis.text.y = element_text(size=14, colour ="black"),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank()) 
mic.eff.sn

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#Fig. 1f plot microbe effects stem density
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
miceff.t<-read.csv(file="plant_trait_microbe_eff_traits_lsmeans.csv")
#miceff.t<-read.csv(file="plant_trait_microbe_eff_traits_lsmeans-exp.csv")
levels(miceff.t$Provenance)[levels(miceff.t$Provenance)=="Corn"]<-"CI"
levels(miceff.t$Provenance)[levels(miceff.t$Provenance)=="Sellman"]<-"SM"
levels(miceff.t$Cohort)[levels(miceff.t$Cohort)=="Modern"]<-"Descendant"

den<-miceff.t[miceff.t$Trait=="density",]
den <- na.omit(den) 

mic.eff.den<-ggplot(den, aes(x = Salinity, y = lsmean, shape= Cohort, color=Cohort)) + 
  geom_line() + #ylim(-25,25)+
  geom_hline(yintercept=0, size=1, lty=2, color='grey50')+
  scale_color_manual(values=c("black","black"))+
  scale_shape_manual(values=c(17,1))+
  scale_x_discrete(limits=c( "Low Salinity","High Salinity"), 
                   labels=c("Low", "High"))+
  geom_errorbar(aes(ymin=lsmean-SE, ymax=lsmean+SE), width=0.01, colour="gray50") +
  geom_point(size=4) +
  #geom_text(aes(label=Cohort), size=4,hjust=-0.5, vjust=0.1)+
  #ggtitle("") +
  facet_grid(~Provenance, scales = "free_x", space = "free_x", shrink=TRUE) +
  #facet_wrap(~Provenance, ncol=2)+
  labs(x="\nSalinity Level", y= "microbe effect SD\n") +
  theme(panel.spacing = unit(0.5, "lines"),  #panel.background = element_blank(),
        strip.background = element_blank(), legend.position='none',
        strip.placement = "inside", plot.title = element_text(hjust = 0.2),
        axis.line.y.left=element_line(color="black"),axis.line.x.bottom=element_line(color="black"),
        strip.text = element_text(face="bold", vjust=1.0, size = 14),
        axis.title.x = element_text(vjust=1.0, size = 12),
        axis.title.y = element_text(vjust=1.0, size = 11),
        axis.text.x = element_text(size=12, colour ="black"),
        axis.text.y = element_text(size=14, colour ="black"),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank()) 
mic.eff.den
#====================#
#Combine
library(ggpubr)
tiff("Fig2.tiff", width = 146, height = 169, units = 'mm', res = 500)
ggarrange(height.plot.all, 
          mic.eff.ph,
          sn.plot.all,
          mic.eff.sn,
          sd.plot.all,
          mic.eff.den,
          labels=c("a","b","c","d","e","f"),
          vjust=1,
          nrow=3, ncol=2)
dev.off()

###############################
#Figure 3 PLOTTING BIOMASS#
#plot by barplot
trait.d<-read.csv(file="All_data_plant_measurement_less_samp_corrected.csv")
levels(trait.d$Site)[levels(trait.d$Site)=="Corn"]<-"CI"
levels(trait.d$Site)[levels(trait.d$Site)=="Sellman"]<-"SM"
levels(trait.d$Cohort)[levels(trait.d$Cohort)=="Modern"]<-"Descendant"
#replace 0 with NAs
trait.d[trait.d==0]<-NA
#Fig 3b. AG BIOMASS CORN
corn.d<-trait.d[trait.d$Site=="CI",]
ag.biom.ave<-ddply(corn.d, c("Mbiome.trt","Salinity", "Cohort"), summarise, N = length (AG_biom), mean = mean(AG_biom, na.rm=TRUE), 
                   sd= sd(AG_biom, na.rm=TRUE), se = var((AG_biom)/sqrt(N), na.rm=TRUE))
ag.biom.corn.plot<-ggplot(subset(ag.biom.ave, Mbiome.trt %in% c("Live", "Sterile")), aes(x=as.factor(interaction(Salinity,Mbiome.trt)), y=mean, color=Cohort, 
                                                                                         fill=Cohort, group=as.factor(interaction(Cohort, Mbiome.trt)))) + 
  scale_fill_manual(values=c("gray10", "gray90"))+
  scale_color_manual(values=c("black", "black"))+ylim(0.0, 0.8)+
  geom_bar(stat="identity", position=position_dodge(), width=0.8)+ 
  ggtitle("CI")+ 
  ylab(bquote(atop('AG Biomass' ~(g~m^2))))+xlab("")+
  scale_x_discrete(limits=c("Low_0ppt.Live","Low_0ppt.Sterile","High_15ppt.Live","High_15ppt.Sterile"),
                   labels=c("Low\n(+) soil-\nmicrobe","Low\n(-) soil-\nmicrobe","High\n(+) soil-\nmicrobe","High\n(-) soil-\nmicrobe"))+
  guides(fill = guide_legend(override.aes = list(linetype = 0)))+
  theme(axis.text=element_text(size=12),
        axis.text.y=element_text(size=11),
        axis.line.y.left=element_line(color="black"),
        axis.line.x.bottom=element_line(color="black"),
        axis.title.y=element_text(size=12),
        legend.position=c(0.54,0.9),legend.direction="horizontal",
        legend.key.size =unit(5, "mm"),legend.title=element_blank(), legend.text = element_text(size=10),
        panel.grid.minor=element_blank(),panel.grid.major=element_blank(), 
        panel.background = element_blank(), plot.title=element_blank())+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(0.75), width=0.2)
ag.biom.corn.plot

#Fig 3c. BG BIOMASS SELLMAN
sell.d<-trait.d[trait.d$Site=="SM",]

bg.biom.ave<-ddply(sell.d, c("Mbiome.trt","Salinity", "Cohort"), summarise, N = length (BG_biom), mean = mean(BG_biom, na.rm=TRUE), 
                   sd= sd(BG_biom, na.rm=TRUE), se = var((BG_biom)/sqrt(N), na.rm=TRUE))

bg.biom.sell.plot<-ggplot(subset(bg.biom.ave, Mbiome.trt %in% c("Live", "Sterile")), aes(x=as.factor(interaction(Salinity,Mbiome.trt)), y=mean, color=Cohort, 
                                                                                         fill=Cohort, group=as.factor(interaction(Cohort, Mbiome.trt)))) + 
  scale_fill_manual(values=c("gray10", "gray90"))+scale_color_manual(values=c("black", "black"))+
  ylim(0.0, 0.4)+
  geom_bar(stat="identity", position=position_dodge(), width=0.8)+ 
  ggtitle("SM")+ 
  ylab(bquote(atop('BG Biomass' ~(g~m^2))))+xlab("")+
  scale_x_discrete(limits=c("Low_0ppt.Live","Low_0ppt.Sterile","High_15ppt.Live","High_15ppt.Sterile"),
                   labels=c("Low\n(+) soil-\nmicrobe","Low\n(-) soil-\nmicrobe","High\n(+) soil-\nmicrobe","High\n(-) soil-\nmicrobe"))+
  guides(fill = guide_legend(override.aes = list(linetype = 0)))+
  theme(axis.text.x=element_text(size=12),
        axis.text.y=element_text(size=11),
        axis.line.y.left=element_line(color="black"),
        axis.line.x.bottom=element_line(color="black"),
        axis.title.y=element_text(size=12),
        legend.position=c(0.54,0.9),legend.direction="horizontal",
        legend.key.size =unit(5, "mm"),legend.title=element_blank(),legend.text = element_text(size=10), 
        panel.grid.minor=element_blank(),panel.grid.major=element_blank(), 
        panel.background = element_blank(), plot.title=element_blank())+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(0.75), width=0.2)
#plot.title=element_text(hjust=0.5, size=14)
bg.biom.sell.plot
#Combine
library(ggpubr)
tiff("Fig3_biomass_25June2022.tiff", width = 195, height = 100, units = 'mm', res = 500)
ggarrange(ag.biom.corn.plot, 
          bg.biom.sell.plot, 
          labels=c("a","b"),
          vjust=1,
          nrow=1, ncol=2)
dev.off()

library(ggplot2)
library(ggrepel)

library(ggplot2)
library(ggrepel)

##########################################################################
#PLOTTING BLUPs
##########################################################################

#===================#
#Fig4a. BG sellman
#===================#

bg_sell<-read.csv(file="BG_blups_sellman_28Oct.csv")

levels(bg_sell$term)[levels(bg_sell$term)=="Live"]<-"Live"

#rename genotype to prevent overlap, note: as of November 15, 2021 somehow, the renaming doesn't work!!!
#it used to work before, I don't get it. It then worked December 15, 2021
#tried again June 24, 2022, didn't work with the updated R
#sellman
levels(bg_sell$Genotype)[levels(bg_sell$Genotype)=="AS1"]<-"A2"
levels(bg_sell$Genotype)[levels(bg_sell$Genotype)=="AS2"]<-"A3"
levels(bg_sell$Genotype)[levels(bg_sell$Genotype)=="MSR1"]<-"M1"
levels(bg_sell$Genotype)[levels(bg_sell$Genotype)=="MS2"]<-"M2"
levels(bg_sell$Genotype)[levels(bg_sell$Genotype)=="MS1"]<-"M3"


bg.sell.plot<-ggplot(bg_sell, aes(x = level, y = estimate, group=Genotype, shape= Cohort, color=Cohort)) + 
  geom_line() + 
  scale_color_manual(values=c("black","black"))+
  #scale_fill_manual(values=c("black","blue"))+
  scale_shape_manual(values=c(17,1))+
  scale_x_discrete(limits=c("Low Salinity", "High Salinity"), 
                   labels=c("Low", "High"))+
  #geom_errorbar(aes(ymin=estimate-std.error, ymax=estimate+std.error), width=0.01, colour="gray50") +
  geom_point(size=3,position = position_jitter(w = 0, h = 0.01)) +
  #geom_text(aes(label=genotype), size=4,hjust=-0.5, vjust=0.05)+
  geom_text_repel(aes(label = Genotype),size=4,hjust=-0.5, vjust=0.05)+
  ggtitle("Belowground") +
  facet_grid(~term, switch = "x", scales = "free_x", space = "free_x") +
  labs(x="", y= "estimate \n") +
  theme(panel.spacing = unit(0.5, "lines"),  #panel.background = element_blank(),
        strip.background = element_blank(), legend.position='none',
        strip.placement = "outside", plot.title = element_text(hjust = 0.5),
        axis.line.y.left=element_line(color="black"),axis.line.x.bottom=element_line(color="black"),
        #strip.text = element_text(face="bold", vjust=-1.5, size = 14),
        strip.text.x = element_blank(),
        axis.title.x = element_text(face="bold", vjust=1.0, size = 16),
        axis.title.y = element_text(face="bold", vjust=1.0, size = 16),
        axis.text.x = element_text(size=14, colour ="black"),
        axis.text.y = element_text(size=16, colour ="black"),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank()) 
bg.sell.plot

#ggsave(filename="blups_bg_sellman_colored_cohort.pdf", plot=bg.sell.plot, dpi=600, units=c("mm"), width=110, height=110)
#ggsave(filename="blups_bg_sellman_colored_cohort.tiff", plot=bg.sell.plot, dpi=600, units=c("mm"), width=110, height=110)

#===================#
#Fig4b. AG sellman
#===================#
ag_sell<-read.csv(file="AG_blups_sellman_28Oct.csv")
levels(ag_sell$term)[levels(ag_sell$term)=="Live"]<-"Live"
#sellman
levels(ag_sell$genotype)[levels(ag_sell$genotype)=="AS1"]<-"A2"
levels(ag_sell$genotype)[levels(ag_sell$genotype)=="AS2"]<-"A3"
levels(ag_sell$genotype)[levels(ag_sell$genotype)=="MSR1"]<-"M1"
levels(ag_sell$genotype)[levels(ag_sell$genotype)=="MS2"]<-"M2"
levels(ag_sell$genotype)[levels(ag_sell$genotype)=="MS1"]<-"M3"

ag.sell.plot<-ggplot(ag_sell, aes(x = level, y = estimate, group=genotype, shape= Cohort, color=Cohort)) + 
  geom_line() + 
  scale_color_manual(values=c("black","black"))+
  #scale_color_manual(values=c("red","blue","red","purple4","orange1"))+
  #scale_fill_manual(values=c("black","blue"))+
  scale_shape_manual(values=c(17,1))+
  scale_x_discrete(limits=c("Low Salinity", "High Salinity"), 
                   labels=c("Low", "High"))+
  #geom_errorbar(aes(ymin=estimate-std.error, ymax=estimate+std.error), width=0.01, colour="gray50") +
  geom_point(size=3) +
  #geom_errorbar(aes(ymin=estimate-std.error, ymax=estimate+std.error), width=0.01, colour="gray50") +
  #geom_text_repel(aes(label = genotype), size=3, point.padding= unit(0.2, "cm"))+
  #geom_text(aes(label=genotype), size=4,hjust=-0.5, vjust=0.05)+
  geom_text_repel(aes(label = genotype),size=4,hjust=-0.5, vjust=0.05)+
  ggtitle("Aboveground") +
  facet_grid(~term, switch = "x", scales = "free_x", space = "free_x") +
  labs(x="", y= "estimate \n") +
  theme(panel.spacing = unit(0.5, "lines"),  #panel.background = element_blank(),
        strip.background = element_blank(), legend.position='none',
        strip.placement = "outside", plot.title = element_text(hjust = 0.5),
        axis.line.y.left=element_line(color="black"),axis.line.x.bottom=element_line(color="black"),
        #strip.text = element_text(face="bold", vjust=-1.5, size = 14),
        strip.text.x = element_blank(),
        axis.title.x = element_text(face="bold", vjust=1.0, size = 16),
        axis.title.y = element_text(face="bold", vjust=1.0, size = 16),
        axis.text.x = element_text(size=14, colour ="black"),
        axis.text.y = element_text(size=16, colour ="black"),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank()) 
ag.sell.plot

#===================#
#Fig4c. Total sellman
#===================#
total_sell<-read.csv(file="total_biom_blups_sellman_new_28Oct.csv")
levels(total_sell$term)[levels(total_sell$term)=="Live"]<-"Live"
#sellman
levels(total_sell$genotype)[levels(total_sell$genotype)=="AS1"]<-"A2"
levels(total_sell$genotype)[levels(total_sell$genotype)=="AS2"]<-"A3"
levels(total_sell$genotype)[levels(total_sell$genotype)=="MSR1"]<-"M1"
levels(total_sell$genotype)[levels(total_sell$genotype)=="MS2"]<-"M2"
levels(total_sell$genotype)[levels(total_sell$genotype)=="MS1"]<-"M3"

new_labels<-c("Live"= "(+) soil-microbe", "Sterile"="(-) soil-microbe")

total.sell.plot<-ggplot(total_sell, aes(x = level, y = estimate, group=genotype, shape= Cohort, color=Cohort)) + 
  geom_line() + 
  scale_color_manual(values=c("black","black"))+
  #scale_fill_manual(values=c("black","blue"))+
  scale_shape_manual(values=c(17,1))+
  scale_x_discrete(limits=c("Low Salinity", "High Salinity"), 
                   labels=c("Low", "High"))+
  #geom_errorbar(aes(ymin=estimate-std.error, ymax=estimate+std.error), width=0.01, colour="gray50") +
  geom_point(size=3) +
  #geom_text(aes(label=genotype), size=4,hjust=-0.5, vjust=0.1)+
  geom_text_repel(aes(label = genotype),size=4,hjust=-0.5, vjust=0.05)+
  ggtitle("Total Biomass") +
  facet_grid(~term, switch = "x", scales = "free_x", space = "free_x", labeller = as_labeller(new_labels)) +
  labs(x="", y= "estimate \n") +
  theme(panel.spacing = unit(0.5, "lines"),#panel.background = element_blank(),  
        strip.background = element_blank(), legend.position='none',
        strip.placement = "outside", plot.title = element_text(hjust = 0.5),
        axis.line.y.left=element_line(color="black"),axis.line.x.bottom=element_line(color="black"),
        #strip.text.x = element_blank(),
        strip.text = element_text(face="bold", vjust=-1.5, size = 14),
        axis.title.x = element_text(face="bold", vjust=1.0, size = 16),
        axis.title.y = element_text(face="bold", vjust=1.0, size = 16),
        axis.text.x = element_text(size=14, colour ="black"),
        axis.text.y = element_text(size=16, colour ="black"),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank()) 
total.sell.plot

#=============================================#
#Fig4 Plant size both site/provenance combined
#==============================================#
size_all<-read.csv(file="plantsize_blups_comb_provenance.csv")

#rename just the sellman
levels(size_all$genotype)[levels(size_all$genotype)=="AS1"]<-"A2"
levels(size_all$genotype)[levels(size_all$genotype)=="AS2"]<-"A3"
levels(size_all$genotype)[levels(size_all$genotype)=="MSR1"]<-"M1"
levels(size_all$genotype)[levels(size_all$genotype)=="MS2"]<-"M2"
levels(size_all$genotype)[levels(size_all$genotype)=="MS1"]<-"M3"
#ci


size_all.plot<-ggplot(size_all, aes(x = level, y = estimate, group=genotype, shape= Cohort, color=Provenance)) + 
  geom_line() + 
  #scale_color_manual(values=c("black","black"))+
  scale_color_manual(values=c("black","purple"))+
  scale_fill_manual(values=c("black","purple"))+
  scale_shape_manual(values=c(17,1))+
  scale_x_discrete(limits=c("Low Salinity", "High Salinity"), 
                   labels=c("Low", "High"))+
  #geom_errorbar(aes(ymin=estimate-std.error, ymax=estimate+std.error), width=0.01, colour="gray50") +
  geom_point(size=3) +
  geom_text(aes(label=genotype), size=4,hjust=-0.5, vjust=0.1)+
  #geom_text_repel(aes(label = genotype),hjust=-0.5, vjust=0.05)+
  ggtitle("Plant Size") +
  facet_grid(~term, switch = "x", scales = "free_x", space = "free_x") +
  labs(x="", y= "") +
  theme(panel.spacing = unit(0.5, "lines"),  #panel.background = element_blank(),
        strip.background = element_blank(), legend.position='none',
        strip.placement = "outside", plot.title = element_text(hjust = 0.5),
        axis.line.y.left=element_line(color="black"),axis.line.x.bottom=element_line(color="black"),
        #strip.text = element_text(face="bold", vjust=-1.5, size = 14),
        strip.text.x = element_blank(),
        axis.title.x = element_text(face="bold", vjust=1.0, size = 16),
        axis.title.y = element_text(face="bold", vjust=1.0, size = 16),
        axis.text.x = element_text(size=14, colour ="black"),
        axis.text.y = element_text(size=16, colour ="black"),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank()) 
size_all.plot

#++++++++++++++++++++++++++++#
#CORN
#++++++++++++++++++++++++++++#
#===================#
#Fig4d. height corn
#===================#
height.corn<-read.csv(file="height_blups_corn_28Oct.csv")
levels(height.corn$term)[levels(height.corn$term)=="Live"]<-"Live"

height.corn.plot<-ggplot(height.corn, aes(x = level, y = estimate, group=genotype, shape= Cohort, color=Cohort)) + 
  geom_line() + 
  scale_color_manual(values=c("black","black"))+
  #scale_color_manual(values=c("red","blue","red","purple4","orange1"))+
  #scale_fill_manual(values=c("black","blue"))+
  scale_shape_manual(values=c(17,1))+
  scale_x_discrete(limits=c("Low Salinity", "High Salinity"), 
                   labels=c("Low", "High"))+
  #geom_errorbar(aes(ymin=estimate-std.error, ymax=estimate+std.error), width=0.01, colour="gray50") +
  geom_point(size=3) +
  #geom_text(aes(label=genotype), size=4,hjust=-0.5, vjust=0.1)+
  geom_text_repel(aes(label = genotype),size=4,hjust=-0.5, vjust=0.05)+
  #geom_text(data=bg_sell.high, label=c("AS1"), size = 3, show.legend = FALSE, hjust = 1,nudge_y = 0.15) +
  ggtitle("Stem height") +
  facet_grid(~term, switch = "x", scales = "free_x", space = "free_x") +
  labs(x="", y= "") +
  theme(panel.spacing = unit(0.5, "lines"),  #panel.background = element_blank(),
        strip.background = element_blank(), legend.position='none',
        strip.placement = "outside", plot.title = element_text(hjust = 0.5),
        axis.line.y.left=element_line(color="black"),axis.line.x.bottom=element_line(color="black"),
        #strip.text = element_text(face="bold", vjust=-1.5, size = 14),
        strip.text.x = element_blank(),
        axis.title.x = element_text(face="bold", vjust=1.0, size = 16),
        axis.title.y = element_text(face="bold", vjust=1.0, size = 16),
        axis.text.x = element_text(size=14, colour ="black"),
        axis.text.y = element_text(size=16, colour ="black"),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank()) 
height.corn.plot


#===================#
#Fig4f. stem diameter corn
#===================#
diam.corn<-read.csv(file="diam_blups_corn_28Oct.csv")
levels(diam.corn$term)[levels(diam.corn$term)=="Inoculated"]<-"Live"

new_labels<-c("Live"= "(+) soil-microbe", "Sterile"="(-) soil-microbe")

diam.corn.plot<-ggplot(diam.corn, aes(x = level, y = estimate, group=genotype, shape= Cohort, color=Cohort)) + 
  geom_line() + 
  scale_color_manual(values=c("black","black"))+
  #scale_color_manual(values=c("red","blue","red","purple4","orange1"))+
  #scale_fill_manual(values=c("black","blue"))+
  scale_shape_manual(values=c(17,1))+
  scale_x_discrete(limits=c("Low Salinity", "High Salinity"), 
                   labels=c("Low", "High"))+
  #geom_errorbar(aes(ymin=estimate-std.error, ymax=estimate+std.error), width=0.01, colour="gray50") +
  geom_point(size=3) +
  #geom_text(aes(label=genotype), size=4,hjust=-0.5, vjust=0.1)+
  geom_text_repel(aes(label = genotype),size=4, hjust=-0.5, vjust=0.05)+
  #geom_text(data=bg_sell.high, label=c("AS1"), size = 3, show.legend = FALSE, hjust = 1,nudge_y = 0.15) +
  ggtitle("Stem diameter") +
  facet_grid(~term, switch = "x", scales = "free_x", space = "free_x", labeller = as_labeller(new_labels)) +
  labs(x="", y= "") +
  theme(panel.spacing = unit(0.5, "lines"),  #panel.background = element_blank(),
        strip.background = element_blank(), legend.position='none',
        strip.placement = "outside", plot.title = element_text(hjust = 0.5),
        axis.line.y.left=element_line(color="black"),axis.line.x.bottom=element_line(color="black"),
        strip.text = element_text(face="bold", vjust=-1.5, size = 14),
        axis.title.x = element_text(face="bold", vjust=1.0, size = 16),
        axis.title.y = element_text(face="bold", vjust=1.0, size = 16),
        axis.text.x = element_text(size=14, colour ="black"),
        axis.text.y = element_text(size=16, colour ="black"),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank()) 
diam.corn.plot

#******************************#
#COMBINE ALL
#******************************#
library(ggpubr)
library(gridExtra)
#vertical
tiff("Fig4_blups_24une2022.tiff", width = 230, height = 330, units = 'mm', res = 300)
ggarrange(ag.sell.plot,
          size_all.plot,
          bg.sell.plot,
          height.corn.plot,
          total.sell.plot,
          diam.corn.plot,
          labels=c("a","d","b","e","c","f"),
          font.label = list(size = 16, face = "bold"),
          vjust=0.95,
          hjust=-1,
          nrow=3, ncol=2)
dev.off()
#===================================================#
#Fig 5a. A panel; fungi by salinity, and ecotype root, live
#===================================================#
endo.fun<-read.csv(file="FinalData_scam_root_fungi.csv")
live.fun<-endo.fun[endo.fun$Mbiome.trt=="Live"  ,]
levels(live.fun$Site)[levels(live.fun$Site)=="Corn"]<-"CI"
levels(live.fun$Site)[levels(live.fun$Site)=="Sellman"]<-"SM"
div.fun.root.live<-ddply(live.fun, c("Salinity", "Site"), summarise, N = length (enspie), mean = mean(enspie, na.rm=TRUE), 
                         sd= sd(enspie, na.rm=TRUE), se = var((enspie)/sqrt(N), na.rm=TRUE))
root.fungi.live.p<-ggplot(div.fun.root.live, aes(x=Salinity, y=mean, shape=Site, group=as.factor(interaction(Site, Salinity)))) + 
  scale_shape_manual(values = c(1,17))+ 
  geom_line()+ geom_point(cex=4)+ 
  labs(y=expression(paste("Fungal diversity  ",italic('ENS'[PIE]))))+xlab("")+scale_x_discrete(limits=c("Low_0ppt","High_15ppt"),
                                                                                               labels=c("Low","High"))+
  theme(axis.text=element_text(size=14), axis.line.y.left=element_line(color="black"),
        axis.line.x.bottom=element_line(color="black"),
        axis.title.y=element_text(size=14),legend.position=c(0.85,0.85),
        legend.key.size =unit(6, "mm"),legend.title=element_blank(), legend.text = element_text(size=10),
        panel.grid.minor=element_blank(),panel.grid.major=element_blank(), 
        panel.background = element_blank(), plot.title=element_text(hjust=0.5, size=16))+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.02)
root.fungi.live.p

#================================#
#Fig 5b PCo
#================================#
endo.fung<-read.csv(file="FinalData_scam_root_fungi.csv")
live.fung<-endo.fung[endo.fung$Mbiome.trt=="Live"  ,]
#make as dummy variable
# make dummy variable for Salinity
live.fung$salinity<-NA
live.fung$salinity[live.fung$Salinity=="Low_0ppt"] <-0
live.fung$salinity[live.fung$Salinity=="High_15ppt"] <-1
with(live.fung, table(salinity))
# make dummy variable for microbiome treatment
live.fung$cohort<-NA
live.fung$cohort[live.fung$Cohort=="Ancestral"] <-0
live.fung$cohort[live.fung$Cohort=="Modern"] <-1
with(live.fung, table(cohort))
# make dummy variable for microbiome treatment
live.fung$site<-NA
live.fung$site[live.fung$Site=="Corn"] <-0
live.fung$site[live.fung$Site=="Sellman"] <-1
with(live.fung, table(site))
#make species matrix
sp.cols.rt <- grep("denovo", names(live.fung))
# Make Bray Curtis dissimilarity matrix 
# Square root transformation, Wisconsin double standardization,#this emphasizes the environmental variables
scam.root.mat<-vegdist((live.fung[,sp.cols.rt]), method="bray", binary=FALSE, 
                       metaMDS=TRUE, sqrt.dist=TRUE)
dbRDA_all.root.1 <- capscale(scam.root.mat ~ Site+Salinity+Cohort+Mbiome.trt, 
                             data=live.fung, na.option=na.omit)
#Null Model
#generate one model with NOTHING to explain the braycurtis dm matrix
dbRDA_all.root.0 <- capscale(scam.root.mat~1,data=live.fung)

#use forward selection to choose which elements of the full model explain a significant amount of variaiton in the unifrac dm by comparing it to the null

dbRDA_all.root <- ordistep(dbRDA_all.root.0,scope=formula(dbRDA_all.root.1),
                           direction="forward",Pin=.1,trace=T,pstep=500000)
dbRDA_all.root

#test significance of analysis
anova(dbRDA_all.root)

#test 
#anova(dbRDA_all.root, by="margin")
#anova(dbRDA_all.root, by="terms", permu=200) # same above,test for sign. environ. variables
# test axes for significance
#anova(dbRDA_all.root, by="axis", perm.max=500) 

#any rid of NA values?
anyNA(live.fung)
#make model sumary
B <- summary(dbRDA_all.root)
#plot
A.1 <- scores(dbRDA_all.root)
A.2 <- A.1$sites
A.3 <- cbind(A.2, live.fung)

#scores for arows
A.4 <- data.frame(scores(dbRDA_all.root, display = "bp"))

#subset A4 for labeling
A.4 <- A.4[sort(rownames(A.4)),]
A4.sub1 <- A.4[2,]
A4.sub2 <- A.4[1,]

#make plot
p <- ggplot(data = A.3, aes(x = CAP1, y = CAP2), bty="n")
p.dbrda.fungi.bac <- p +
  geom_point(data = A.3, alpha = 2/5, size=4,
             aes(shape = Site, color=Salinity, stroke = 1)) +
  theme_bw() + 
  scale_color_manual(name = "", labels = c("High Salinity", "Low Salinity"),values=c("black", "red"))+
  #scale_shape_manual(values=c(15,17,0,2))+
  #scale_fill_manual(name = "", labels = c("High Salinity", "Low Salinity"),values=c("black", "red"))+
  xlab(label = paste("CAP1 (", round(B$concont$importance[2,1]*100, digits = 1), "%)", sep="")) +
  ylab(label = paste("CAP2 (", round(B$concont$importance[2,2]*100, digits = 1), "%)", sep="")) +
  scale_shape_manual(name = "Provenance", labels = c("CI", "SM"), 
                     values = c(1,17))
p.dbrda.fungi.bac


#================================#
#Fig 5c panel bacteria endosphere
#================================#
endo.bac<-read.csv(file="FinalData_scam_root_bacteria.csv")
live.bac<-endo.bac[endo.bac$Mbiome.trt=="Live"  ,]
levels(live.bac$Site)[levels(live.bac$Site)=="Corn"]<-"CI"
levels(live.bac$Site)[levels(live.bac$Site)=="Sellman"]<-"SM"

div.root.bac<-ddply(live.bac, c("Salinity","Site"), summarise, N = length (enspie), mean = mean(enspie, na.rm=TRUE), 
                    sd= sd(enspie, na.rm=TRUE), se = var((enspie)/sqrt(N), na.rm=TRUE))
#
div.root.bac.plot<-ggplot(div.root.bac, aes(x=Salinity, y=mean, shape=Site,group=as.factor(interaction(Site, Salinity)))) + 
  scale_shape_manual(values = c(1,17))+ 
  geom_line()+ geom_point(cex=4)+ 
  labs(y=expression(paste("Bacterial diversity  ",italic(  'ENS'[PIE]))))+xlab("")+scale_x_discrete(limits=c("Low_0ppt","High_15ppt"),
                                                                                                    labels=c("Low","High"))+
  theme(axis.text=element_text(size=14), axis.line.y.left=element_line(color="black"),
        axis.line.x.bottom=element_line(color="black"),
        axis.title.y=element_text(size=14),legend.position=c(0.85,0.85),
        legend.key.size =unit(6, "mm"),legend.title=element_blank(), legend.text = element_text(size=10),
        panel.grid.minor=element_blank(),panel.grid.major=element_blank(), 
        panel.background = element_blank(), plot.title=element_text(hjust=0.5, size=16))+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.02)

div.root.bac.plot
#ggsave(filename="bacteria_root_all.tiff", plot=div.root.bac.plot, dpi=400, units=c("mm"), width=100, height=100)


#================================#
#Fig 5d PCo
#================================#
endo.bac<-read.csv(file="FinalData_scam_root_bacteria.csv")
live.bac<-endo.bac[endo.bac$Mbiome.trt=="Live"  ,]

#make as dummy variable
# make dummy variable for Salinity
live.bac$salinity<-NA
live.bac$salinity[live.bac$Salinity=="Low_0ppt"] <-0
live.bac$salinity[live.bac$Salinity=="High_15ppt"] <-1
with(live.bac, table(salinity))
# make dummy variable for microbiome treatment
live.bac$cohort<-NA
live.bac$cohort[live.bac$Cohort=="Ancestral"] <-0
live.bac$cohort[live.bac$Cohort=="Modern"] <-1
with(live.bac, table(cohort))
# make dummy variable for microbiome treatment
live.bac$site<-NA
live.bac$site[live.bac$Site=="Corn"] <-0
live.bac$site[live.bac$Site=="Sellman"] <-1
with(live.bac, table(site))
#make species matrix
sp.cols.rt.bac <- grep("otu_", names(live.bac))
# Make Bray Curtis dissimilarity matrix 
# Square root transformation, Wisconsin double standardization,#this emphasizes the environmental variables
scam.root.mat.bac<-vegdist((live.bac[,sp.cols.rt.bac]), method="bray", binary=FALSE, 
                           metaMDS=TRUE, sqrt.dist=TRUE)
dbRDA_all.root.bac1 <- capscale(scam.root.mat.bac ~ Site+Salinity+Cohort+Mbiome.trt, 
                                data=live.bac, na.option=na.omit)
#Null Model
#generate one model with NOTHING to explain the braycurtis dm matrix
dbRDA_all.root.bac0 <- capscale(scam.root.mat.bac~1,data=live.bac)

#use forward selection to choose which elements of the full model explain a significant amount of variaiton in the unifrac dm by comparing it to the null

dbRDA_all.root.bac <- ordistep(dbRDA_all.root.bac0,scope=formula(dbRDA_all.root.bac1),
                               direction="forward",Pin=.1,trace=T,pstep=500000)
dbRDA_all.root.bac

#test significance of analysis
anova(dbRDA_all.root.bac)

#test 
#anova(dbRDA_all.root.bac, by="margin")
#anova(dbRDA_all.root.bac, by="terms", permu=200) # same above,test for sign. environ. variables
# test axes for significance
#anova(dbRDA_all.root.bac, by="axis", perm.max=500) 

#any rid of NA values?
anyNA(live.bac)
#make model sumary
B.bac <- summary(dbRDA_all.root.bac)
#plot
Abac.1 <- scores(dbRDA_all.root.bac)
Abac.2 <- Abac.1$sites
Abac.3 <- cbind(Abac.2, live.bac)

#scores for arows
Abac.4 <- data.frame(scores(dbRDA_all.root.bac, display = "bp"))

#subset A4 for labeling
Abac.4 <- Abac.4[sort(rownames(Abac.4)),]
Abac.4.sub1 <- Abac.4[2,]
A.4.sub2 <- Abac.4[1,]

#make plot
p.bac <- ggplot(data = Abac.3, aes(x = CAP1, y = CAP2), bty="n")
p.dbrda.bac.plot <- p.bac +
  geom_point(data = Abac.3, alpha = 2/5, size=4,
             aes(shape = Site, color=Salinity, stroke = 1)) +
  theme_bw() + 
  scale_color_manual(name = "", labels = c("High Salinity", "Low Salinity"),values=c("black", "red"))+
  #scale_shape_manual(values=c(15,17,0,2))+
  #scale_fill_manual(name = "", labels = c("High Salinity", "Low Salinity"),values=c("black", "red"))+
  xlab(label = paste("CAP1 (", round(B.bac$concont$importance[2,1]*100, digits = 1), "%)", sep="")) +
  ylab(label = paste("CAP2 (", round(B.bac$concont$importance[2,2]*100, digits = 1), "%)", sep="")) +
  scale_shape_manual(name = "Provenance", labels = c("CI", "SM"), 
                     values = c(1,17))
p.dbrda.bac.plot

#******************************#
#COMBINE ALL
#******************************#
library(ggpubr)
library(gridExtra)


#exporting to pdf or png
#below is much easier
# Export to pdf
ggarrange(root.fungi.live.p, 
          p.dbrda.fungi.bac,
          div.root.bac.plot,
          p.dbrda.bac.plot,
          labels=c("a","b","c","d"),
          nrow=2, ncol=2, widths=c(1,1.5))%>%
  ggexport(filename = "Figure5_mec.pdf")

#---------------------------------------------#
#Fig 6b. Fungi combined height and SD residual
#---------------------------------------------#
#read height data
df.fun<-read.csv(file="residual_score_plant_height_fungi.csv")
#read SD
df.fun.gb<-read.csv(file="residual_score_gb_fungi.csv")
#read plant size
df.fun.size<-read.csv(file="residual_score_size_fungi.csv")
#read gb
#df.fun.gb<-read.csv(file="residual_score_gb_fungi.csv")

#combine both in singe plot
ht.fung.p<-ggplot() + 
  geom_point(data = df.fun, aes(x=CAP1, y=resid.ht.mod), pch=2, cex=1.8) +
  geom_point(data = df.fun.size, aes(x=CAP1, y=resid.size.mod0), pch=15, cex=2)+
  geom_point(data = df.fun.gb, aes(x=CAP1, y=resid.lm.mod.grn2), pch=10, cex=2)+
  geom_smooth(data=df.fun.gb, aes(x=CAP1, y=resid.lm.mod.grn2), method = "lm", se=TRUE, color="black", formula= y ~ x) + 
  ylab("Residual\n")+xlab("\nPCo1")+ ylim(-2.5,2.0)+
  geom_smooth(data=df.fun.gb, aes(x=CAP1, y=resid.lm.mod.grn2), method = "lm", se=TRUE, color="black", formula= y ~ x, lty=1) + 
  theme(axis.text=element_text(size=16), axis.line.y.left=element_line(color="black"),
        axis.line.x.bottom=element_line(color="black"),
        axis.title=element_text(size=16),legend.position=c(0.85,0.85),
        legend.key.size =unit(6, "mm"),legend.title=element_blank(), legend.text = element_text(size=16),
        panel.grid.minor=element_blank(),panel.grid.major=element_blank(), 
        panel.background = element_blank(), plot.title=element_text(hjust=0.5, size=16))
ht.fung.p
ggsave(filename="fungi_residPco2_comb_new.tiff", plot=ht.fung.p, width = 90, height = 90, units = 'mm',dpi=500)

#---------------------------------------------#
#Fig 6d. Bacteria combined height and SD residual
#---------------------------------------------#
#height
df.bac<-read.csv(file="residual_score_plant_height_bacteria.csv")
#diam
df.bac.size<-read.csv(file="residual_score_size_bacteria.csv")
#read plant size
df.bac.num<-read.csv(file="residual_score_number_bacteria.csv")

#combine both in single plot
ht.bac.p<-ggplot() + 
  geom_point(data = df.bac, aes(x=CAP2, y=resid.ht.mod), pch=2, cex=1.5) +
  geom_point(data = df.bac.size, aes(x=CAP2, y=resid.size.mod0), pch=15, cex=2)+
  geom_point(data = df.bac.num, aes(x=CAP2, y=resid.sn.mod), pch=4, cex=2)+
  ylab("Residual\n")+xlab("\nPCo2")+ ylim(-1.5,1)+
  geom_smooth(data=df.bac, aes(x=CAP2, y=resid.ht.mod), method = "lm", se=TRUE, color="black", formula= y ~ x) + 
  geom_smooth(data=df.bac.size, aes(x=CAP2, y=resid.size.mod0), method = "lm", se=TRUE, color="red", formula= y ~ x, lty=2) + 
  geom_smooth(data=df.bac.num, aes(x=CAP2, y=resid.sn.mod), method = "lm", se=TRUE, color="blue", formula= y ~ x, lty=3) + 
  theme(axis.text=element_text(size=16), axis.line.y.left=element_line(color="black"),
        axis.line.x.bottom=element_line(color="black"),
        axis.title=element_text(size=16),legend.position=c(0.85,0.85),
        legend.key.size =unit(6, "mm"),legend.title=element_blank(), legend.text = element_text(size=16),
        panel.grid.minor=element_blank(),panel.grid.major=element_blank(), 
        panel.background = element_blank(), plot.title=element_text(hjust=0.5, size=16))
ht.bac.p
ggsave(filename="bacteria_residPco2_comb_new_trait.tiff", plot=ht.bac.p, width = 90, height = 90, units = 'mm',dpi=500)

dev.off()