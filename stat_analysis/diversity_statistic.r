library(lme4)
library(lmerTest)
library(emmeans)
library(rsq)
data_all <- read.delim("E:/Hic_vcf/re_analysis/ld_diversity/analysis_data.txt")

diversity <- data_all[,c(2,3,4)]
scale_div <- scale(diversity)
pca_result <- princomp(scale_div )
data_all$PC1 <- pca_result$scores[,1]


mod_rand_PC1 <- lmer(data=data_all,PC1~(1|clade))
data_all$res_PC1 <- residuals(mod_rand_PC1)
mod_rand_Pi <- lmer(data=data_all,Pi~(1|clade))
data_all$res_Pi <- residuals(mod_rand_Pi)
mod_rand_He <- lmer(data=data_all,He~(1|clade))
data_all$res_He <- residuals(mod_rand_He)
mod_rand_theta_w <- lmer(data=data_all,theta_w~(1|clade))
data_all$res_theta_w <- residuals(mod_rand_theta_w)


#############
#   ANOVA   #
#############


oneway.test(res_PC1 ~ bio_geo, data = data_all,var.equal = TRUE)

#       One-way analysis of means
#
#       data:  res_PC1 and bio_geo
#       F = 7.6893, num df = 3, denom df = 26, p-value = 0.0007738

########### posthoc
with(data_all,pairwise.t.test(res_PC1, bio_geo, pool.sd=FALSE, p.adjust.method="none"))

#       Pairwise comparisons using t tests with non-pooled SD 
#
#       data:  res_PC1 and bio_geo 
#
#             continu_area core_refug leading_edge
#       core_refug   0.5572       -          -           
#       leading_edge 0.0012       0.0014     -           
#       rear_edge    0.0230       0.0174     0.1915      
#
#       P value adjustment method: none 

oneway.test(res_He ~ bio_geo, data = data_all,var.equal = TRUE)

#       One-way analysis of means
#
#       data:  res_He and bio_geo
#       F = 8.4972, num df = 3, denom df = 26, p-value = 0.0004228

########### posthoc
with(data_all,pairwise.t.test(res_He, bio_geo, pool.sd=FALSE, p.adjust.method="none"))

#       Pairwise comparisons using t tests with non-pooled SD 
#
#       data:  res_He and bio_geo 
#
#             continu_area core_refug leading_edge
#       core_refug   0.858        -          -           
#       leading_edge 0.022        0.037      -           
#       rear_edge    0.012        0.010      0.058       
#
#       P value adjustment method: none 


oneway.test(res_Pi ~ bio_geo, data = data_all,var.equal = TRUE)

#       One-way analysis of means
#
#       data:  res_Pi and bio_geo
#       F = 12.495, num df = 3, denom df = 26, p-value = 3.003e-05

########### posthoc
with(data_all,pairwise.t.test(res_Pi, bio_geo, pool.sd=FALSE, p.adjust.method="none"))

#       Pairwise comparisons using t tests with non-pooled SD 
#
#       data:  res_Pi and bio_geo 
#
#             continu_area core_refug leading_edge
#       core_refug   0.27398      -          -           
#       leading_edge 0.00068      0.00023    -           
#       rear_edge    0.00866      0.00628    0.06273     
#
#       P value adjustment method: none



oneway.test(res_theta_w ~ bio_geo, data = data_all,var.equal = TRUE)

#       One-way analysis of means
#
#       data:  res_theta_w and bio_geo
#       F = 3.5079, num df = 3, denom df = 26, p-value = 0.02931

########### posthoc
with(data_all,pairwise.t.test(res_theta_w, bio_geo, pool.sd=FALSE, p.adjust.method="none"))

#       Pairwise comparisons using t tests with non-pooled SD 
#
#       data:  res_theta_w and bio_geo 
#
#             continu_area core_refug leading_edge
#       core_refug   0.72366      -          -           
#       leading_edge 0.00079      0.00107    -           
#       rear_edge    0.13521      0.11640    0.86058     
#
#       P value adjustment method: none



########################
#    founder effect    #
########################


data_lead <- data_all[which(data_all$group != "RE"),]

mod_hist_pc1 <- lmer(data=data_lead,PC1~hist_periphel*group+(1|clade))
summary(mod_hist_pc1)
rsq.glmm(mod_hist_pc1)
post_trends <- emtrends(mod_hist_pc1,pairwise~group,var="hist_periphel")
test(post_trends)

mod_hist_he <- lmer(data=data_lead,He~hist_periphel*group+(1|clade))
summary(mod_hist_he)
rsq.glmm(mod_hist_he)
post_trends <- emtrends(mod_hist_he,pairwise~group,var="hist_periphel")
test(post_trends)

mod_hist_pi <- lmer(data=data_lead,Pi~hist_periphel*group+(1|clade))
summary(mod_hist_pi)
rsq.glmm(mod_hist_pi)
post_trends <- emtrends(mod_hist_pi,pairwise~group,var="hist_periphel")
test(post_trends)

mod_hist_w <- lmer(data=data_lead,theta_w~hist_periphel*group+(1|clade))
summary(mod_hist_w)
rsq.glmm(mod_hist_w)
post_trends <- emtrends(mod_hist_w,pairwise~group,var="hist_periphel")
test(post_trends)



####################################
#    habitat suitability effect    #
####################################


mod_pc1_all <- lm(data=data_all,PC1~suitability)
summary(mod_pc1_all)
mod_he_all <- lm(data=data_all,He~suitability)
summary(mod_he_all)
mod_pi_all <- lm(data=data_all,Pi~suitability)
summary(mod_pi_all)
mod_theta_all <- lm(data=data_all,theta_w~suitability)
summary(mod_theta_all)

mod_pc1_all <- lm(data=data_all,PC1~brood_size)
summary(mod_pc1_all)
mod_he_all <- lm(data=data_all,He~brood_size)
summary(mod_he_all)
mod_pi_all <- lm(data=data_all,Pi~brood_size)
summary(mod_pi_all)
mod_theta_all <- lm(data=data_all,theta_w~brood_size)
summary(mod_theta_all)


mod_suit_pc1 <- lmer(data=data_all,PC1~suitability*group+(1|clade))
summary(mod_suit_pc1)
rsq.glmm(mod_suit_pc1)
post_trends <- emtrends(mod_suit_pc1,pairwise~group,var="suitability")
test(post_trends)

mod_suit_he <- lmer(data=data_all,He~suitability*group+(1|clade))
summary(mod_suit_he)
rsq.glmm(mod_suit_he)
post_trends <- emtrends(mod_suit_he,pairwise~group,var="suitability")
test(post_trends)

mod_suit_pi <- lmer(data=data_all,Pi~suitability*group+(1|clade))
summary(mod_suit_pi)
rsq.glmm(mod_suit_pi)
post_trends <- emtrends(mod_suit_pi,pairwise~group,var="suitability")
test(post_trends)

mod_suit_w <- lmer(data=data_all,theta_w~suitability*group+(1|clade))
summary(mod_suit_w)
rsq.glmm(mod_suit_w)
post_trends <- emtrends(mod_suit_w,pairwise~group,var="suitability")
test(post_trends)


mod_brod_pc1 <- lmer(data=data_all,PC1~brood_size*group+(1|clade))
summary(mod_brod_pc1)
rsq.glmm(mod_brod_pc1)
post_trends <- emtrends(mod_brod_pc1,pairwise~group,var="brood_size")
test(post_trends)

mod_brod_he <- lmer(data=data_all,He~brood_size*group+(1|clade))
summary(mod_brod_he)
rsq.glmm(mod_brod_he)
post_trends <- emtrends(mod_brod_he,pairwise~group,var="brood_size")
test(post_trends)

mod_brod_pi <- lmer(data=data_all,Pi~brood_size*group+(1|clade))
summary(mod_brod_pi)
rsq.glmm(mod_brod_pi)
post_trends <- emtrends(mod_brod_pi,pairwise~group,var="brood_size")
test(post_trends)

mod_brod_w <- lmer(data=data_all,theta_w~brood_size*group+(1|clade))
summary(mod_brod_w)
rsq.glmm(mod_brod_w)
post_trends <- emtrends(mod_brod_w,pairwise~group,var="brood_size")
test(post_trends)


#############################
#   link selection effect   #
#############################

mod <- lm(data=analysis_data,fec_pi~suitability+genome_pi)
summary(mod)
