

setwd("/Users/jerry/Documents/CSBQ/hijri/Acadian_seaplants/")
library(nlme)
source("/Library/Frameworks/R.framework/Versions/3.4/Resources/library/Legendre_numericalecology/anova.2way.R")

####setting things up -----
productivity = read.table("results/plant_productivity_data.tsv", header = T, stringsAsFactors = F,sep = "\t")
productivity.norm = productivity


#we do a sqrt transformation on all variables to help in normalizing (it helps in all cases, but it's not perfect. Some residuals are still not normally distributed. But the effect is so strong that we don't care too much.)
#####
#NO NORMALITY
productivity.norm$fruit.number = sqrt(productivity$fruit.number)
lmm.fruit.number <- lme(fruit.number~treatment*plant,data = productivity.norm,random = ~1|block/replicate, method = "REML")
anova(lmm.fruit.number)
shapiro.test(lmm.fruit.number$residuals[,1]) #test the normality of the residualds for the fixed effects...

    #So we do an ANOVA with premutation... results are the same.
anova.2way(fruit.number~treatment*plant,data = productivity.norm,nperm = 999,model = 1)


#NO NORMALITY
productivity.norm$average.fruit.weight = sqrt(productivity$average.fruit.weight)
lmm.average.fruit.weight <- lme(average.fruit.weight~treatment*plant,data = productivity.norm,random = ~1|block/replicate, method = "REML")
anova(lmm.average.fruit.weight)
shapiro.test(lmm.average.fruit.weight$residuals[,1]) #test the normality of the residualds for the fixed effects...

          productivity.norm.p = productivity.norm[productivity.norm[,1]== "Pepper",] #test only peppers
          lmm.average.fruit.weight.p <- lme(average.fruit.weight~treatment,data = productivity.norm.p,random = ~1|block/replicate, method = "REML")
          anova(lmm.average.fruit.weight.p)

          #So we do an ANOVA with premutation... results are the same.
          anova.2way(average.fruit.weight~treatment*plant,data = productivity.norm,nperm = 999,model = 1)
          
          
#NORMALITY LOG
productivity.norm$shoots.fresh.weight = log(productivity$shoots.fresh.weight)
lmm.shoots.fresh.weight <- lme(shoots.fresh.weight~treatment*plant,data = productivity.norm,random = ~1|block/replicate, method = "REML")
anova(lmm.shoots.fresh.weight)
shapiro.test(lmm.shoots.fresh.weight$residuals[,1]) #test the normality of the residualds for the fixed effects...

#NORMALITY SQRT
productivity.norm$shoots.dry.weight = sqrt(productivity$shoots.dry.weight)
lmm.shoots.dry.weight <- lme(shoots.dry.weight~treatment*plant,data = productivity.norm,random = ~1|block/replicate, method = "REML")
anova(lmm.shoots.dry.weight)
shapiro.test(lmm.shoots.dry.weight$residuals[,1]) #test the normality of the residualds for the fixed effects...

#NORMALITY LOG
productivity.norm$roots.fresh.weight = sqrt(productivity$roots.fresh.weight)
lmm.roots.fresh.weight <- lme(roots.fresh.weight~treatment*plant,data = productivity.norm,random = ~1|block/replicate, method = "REML")
anova(lmm.roots.fresh.weight)
shapiro.test(lmm.roots.fresh.weight$residuals[,1]) #test the normality of the residualds for the fixed effects...

#NORMALITY SQRT
productivity.norm$roots.dry.weight = sqrt(productivity$roots.dry.weight)
lmm.roots.dry.weight <- lme(roots.dry.weight~treatment*plant,data = productivity.norm,random = ~1|block/replicate, method = "REML")
anova(lmm.roots.dry.weight)
shapiro.test(lmm.roots.dry.weight$residuals[,1]) #test the normality of the residualds for the fixed effects...

#boxplot
#boxplot it?
dev.new(width=14, height=10,units = "cm",noRStudioGD = TRUE)
par(mfrow = c(2,3),mar = c(6,5,4,2))
par(mgp = c(3.5,2,1),cex.lab = 1.5)

names = expression(italic(control),italic(fertilized),italic(control),italic(fertilized))

x=boxplot(fruit.number~treatment*plant,data = productivity,ylim = c(0,max(productivity$fruit.number)*1.0),names = names,font = 2,cex.main = 1.5, main = "Fruit number")
abline(v = 2.5,lty = 2)
mtext(text = c("Pepper","Tomato"),at = c(1.5,3.5),side = 1,line = -25,cex = 1.5)
text(labels = c("a","b","a","b"),x = c(1.2,2.2,3.2,4.2),y = x$stats[5,]+max(x$stats[5,])/20,cex = 1.5,font =3)

x=boxplot(average.fruit.weight~treatment*plant,data = productivity,ylim = c(0,max(productivity$average.fruit.weight)*1.0),names = names,font = 2,cex.main = 1.5, main = "Avg fruit fresh weight",ylab = "grams")
abline(v = 2.5,lty = 2)
mtext(text = c("Pepper","Tomato"),at = c(1.5,3.5),side = 1,line = -25,cex = 1.5)
#text(labels = c("a","b","a","b"),x = c(1.2,2.2,3.2,4.2),y = x$stats[5,]+max(x$stats[5,])/20,cex = 1.5,font =3)

x=boxplot(shoots.fresh.weight~treatment*plant,data = productivity,ylim = c(0,max(productivity$shoots.fresh.weight)*1.0),names = names,font = 2,cex.main = 1.5, main = "Shoot fresh weight",ylab = "grams")
abline(v = 2.5,lty = 2)
mtext(text = c("Pepper","Tomato"),at = c(1.5,3.5),side = 1,line = -25,cex = 1.5)
text(labels = c("a","b","a","b"),x = c(1.2,2.2,3.2,4.2),y = x$stats[5,]+max(x$stats[5,])/20,cex = 1.5,font =3)

x=boxplot(shoots.dry.weight~treatment*plant,data = productivity,ylim = c(0,max(productivity$shoots.dry.weight)*1.0),names = names,font = 2,cex.main = 1.5, main = "Shoot dry weight",ylab = "grams")
abline(v = 2.5,lty = 2)
mtext(text = c("Pepper","Tomato"),at = c(1.5,3.5),side = 1,line = -25,cex = 1.5)
text(labels = c("a","b","a","b"),x = c(1.2,2.2,3.2,4.2),y = x$stats[5,]+max(x$stats[5,])/20,cex = 1.5,font =3)

x=boxplot(roots.fresh.weight~treatment*plant,data = productivity,ylim = c(0,max(productivity$roots.fresh.weight)*1.0),names = names,font = 2,cex.main = 1.5, main = "Root fresh weight",ylab = "grams")
abline(v = 2.5,lty = 2)
mtext(text = c("Pepper","Tomato"),at = c(1.5,3.5),side = 1,line = -25,cex = 1.5)
text(labels = c("a","b","a","b"),x = c(1.2,2.2,3.2,4.2),y = x$stats[5,]+max(x$stats[5,])/20,cex = 1.5,font =3)

x=boxplot(roots.dry.weight~treatment*plant,data = productivity,ylim = c(0,max(productivity$roots.dry.weight)*1.05),names = names,font = 2,cex.main = 1.5, main = "Root dry weight",ylab = "grams")
abline(v = 2.5,lty = 2)
mtext(text = c("Pepper","Tomato"),at = c(1.5,3.5),side = 1,line = -25,cex = 1.5)
text(labels = c("a","b","a","b"),x = c(1.2,2.2,3.2,4.2),y = x$stats[5,]+max(x$stats[5,])/20,cex = 1.5,font =3)

dev.print(device=pdf, "figures/Figure_3_productivity.pdf", onefile=FALSE)
dev.off()

#legend
system("touch figures/legend") 
system("echo 'Figure_2_productivity.pdf : all 6 traits are different (higher) in fertilized plants (Linear mixed effect model (block and replicate are random factors) p<0.0001). Species are different too (p<0.001) except for the fruit fresh weight' >>figures/legend")
