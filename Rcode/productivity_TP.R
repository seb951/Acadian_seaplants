

setwd("/Users/jerry/Documents/CSBQ/hijri/Acadian_seaplants/")
library(nlme)
source("/Library/Frameworks/R.framework/Versions/3.4/Resources/library/Legendre_numericalecology/anova.2way.R")

####setting things up -----
productivity = read.table("results/plant_productivity_data.tsv", header = T, stringsAsFactors = F,sep = "\t")
productivity.norm = productivity
productivity_all = productivity

#we do a sqrt transformation on all variables to help in normalizing (it helps in all cases, but it's not perfect. Some residuals are still not normally distributed. But the effect is so strong that we don't care too much.)
#####
#NO NORMALITY
for(species in c("Tomato","Pepper"))
{
#per species
productivity.norm = productivity_all[productivity_all[,1] == species,]
productivity = productivity_all[productivity_all[,1] == species,]


productivity.norm$fruit.number = sqrt(productivity$fruit.number)
lmm.fruit.number <- lme(fruit.number~treatment,data = productivity.norm,random = ~1|block/replicate, method = "REML")
print(anova(lmm.fruit.number))
print(shapiro.test(lmm.fruit.number$residuals[,1])) #test the normality of the residualds for the fixed effects...

    #So we do an ANOVA with premutation... results are the same.
#anova.2way(fruit.number~treatment,data = productivity.norm,nperm = 999,model = 1)


#NO NORMALITY
productivity.norm$average.fruit.weight = sqrt(productivity$average.fruit.weight)
lmm.average.fruit.weight <- lme(average.fruit.weight~treatment,data = productivity.norm,random = ~1|block/replicate, method = "REML")
print(anova(lmm.average.fruit.weight))
print(shapiro.test(lmm.average.fruit.weight$residuals[,1])) #test the normality of the residualds for the fixed effects...

          #So we do an ANOVA with premutation... results are the same.
#         anova.2way(average.fruit.weight~treatment,data = productivity.norm,nperm = 999,model = 1)
          
          
#NORMALITY LOG
productivity.norm$shoots.fresh.weight = log(productivity$shoots.fresh.weight)
lmm.shoots.fresh.weight <- lme(shoots.fresh.weight~treatment,data = productivity.norm,random = ~1|block/replicate, method = "REML")
print(anova(lmm.shoots.fresh.weight))
print(shapiro.test(lmm.shoots.fresh.weight$residuals[,1])) #test the normality of the residualds for the fixed effects...

#NORMALITY SQRT
productivity.norm$shoots.dry.weight = sqrt(productivity$shoots.dry.weight)
lmm.shoots.dry.weight <- lme(shoots.dry.weight~treatment,data = productivity.norm,random = ~1|block/replicate, method = "REML")
print(anova(lmm.shoots.dry.weight))
print(shapiro.test(lmm.shoots.dry.weight$residuals[,1])) #test the normality of the residualds for the fixed effects...

#NORMALITY LOG
productivity.norm$roots.fresh.weight = sqrt(productivity$roots.fresh.weight)
lmm.roots.fresh.weight <- lme(roots.fresh.weight~treatment,data = productivity.norm,random = ~1|block/replicate, method = "REML")
print(anova(lmm.roots.fresh.weight))
print(shapiro.test(lmm.roots.fresh.weight$residuals[,1])) #test the normality of the residualds for the fixed effects...

#NORMALITY SQRT
productivity.norm$roots.dry.weight = sqrt(productivity$roots.dry.weight)
lmm.roots.dry.weight <- lme(roots.dry.weight~treatment,data = productivity.norm,random = ~1|block/replicate, method = "REML")
print(anova(lmm.roots.dry.weight))
print(shapiro.test(lmm.roots.dry.weight$residuals[,1])) #test the normality of the residualds for the fixed effects...

}

#get productivity fold changes as well...
p_p = productivity_all[productivity_all[,1]=="Pepper",]
p_t = productivity_all[productivity_all[,1]=="Tomato",]

#boxplot
#boxplot it?
dev.new(width=14, height=10,units = "cm",noRStudioGD = TRUE)
par(mfrow = c(2,3),mar = c(6,5,4,2))
par(mgp = c(3.5,2,1),cex.lab = 1.5)

names = expression(italic(control),italic(amendment),italic(control),italic(amendment))

x=boxplot(fruit.number~treatment*plant,data = productivity_all,ylim = c(0,max(productivity_all$fruit.number)*1.0),show.names=F, font = 2,cex.main = 1.5, main = "Fruit number")
abline(v = 2.5,lty = 2)
#Fold changes
p_fc = signif(mean(p_p$fruit.number[p_p[,2]== "F+"])/mean(p_p$fruit.number[p_p[,2]== "F-"]),3)
t_fc = signif(mean(p_t$fruit.number[p_t[,2]== "F+"])/mean(p_t$fruit.number[p_t[,2]== "F-"]),3)
mtext(text = c(paste("Pepper\n(",p_fc,"X)",sep=""),paste("Tomato\n(",t_fc,"X)",sep="")),at = c(1.5,3.5),side = 1,line = -23,cex = 1.5)
text(labels = c("a","b","a","b"),x = c(1.2,2.2,3.2,4.2),y = x$stats[5,]+max(x$stats[5,])/20,cex = 1.5,font =3)
text(labels = names,x = c(1,2,3,4),c(-2,-2,-2,-2),srt = 60,xpd = T,adj = 1)

x=boxplot(average.fruit.weight~treatment*plant,data = productivity_all,ylim = c(0,max(productivity_all$average.fruit.weight)*1.0),names = names,font = 2,cex.main = 1.5, main = "Avg fruit fresh weight",ylab = "grams")
abline(v = 2.5,lty = 2)
p_fc = signif(mean(p_p$average.fruit.weight[p_p[,2]== "F+"])/mean(p_p$average.fruit.weight[p_p[,2]== "F-"]),3)
t_fc = signif(mean(p_t$average.fruit.weight[p_t[,2]== "F+"])/mean(p_t$average.fruit.weight[p_t[,2]== "F-"]),3)
mtext(text = c(paste("Pepper\n(",p_fc,"X)",sep=""),"Tomato"),at = c(1.5,3.5),side = 1,line = c(-23,-25),cex = 1.5)
#text(labels = c("a","b","a","b"),x = c(1.2,2.2,3.2,4.2),y = x$stats[5,]+max(x$stats[5,])/20,cex = 1.5,font =3)

x=boxplot(shoots.fresh.weight~treatment*plant,data = productivity_all,ylim = c(0,max(productivity_all$shoots.fresh.weight)*1.0),names = names,font = 2,cex.main = 1.5, main = "Shoot fresh weight",ylab = "grams")
abline(v = 2.5,lty = 2)
#Fold changes
p_fc = signif(mean(p_p$shoots.fresh.weight[p_p[,2]== "F+"])/mean(p_p$shoots.fresh.weight[p_p[,2]== "F-"]),3)
t_fc = signif(mean(p_t$shoots.fresh.weight[p_t[,2]== "F+"])/mean(p_t$shoots.fresh.weight[p_t[,2]== "F-"]),3)
mtext(text = c(paste("Pepper\n(",p_fc,"X)",sep=""),paste("Tomato\n(",t_fc,"X)",sep="")),at = c(1.5,3.5),side = 1,line = -23,cex = 1.5)
text(labels = c("a","b","a","b"),x = c(1.2,2.2,3.2,4.2),y = x$stats[5,]+max(x$stats[5,])/20,cex = 1.5,font =3)

x=boxplot(shoots.dry.weight~treatment*plant,data = productivity_all,ylim = c(0,max(productivity_all$shoots.dry.weight)*1.0),names = names,font = 2,cex.main = 1.5, main = "Shoot dry weight",ylab = "grams")
abline(v = 2.5,lty = 2)
#Fold changes
p_fc = signif(mean(p_p$shoots.dry.weight[p_p[,2]== "F+"])/mean(p_p$shoots.dry.weight[p_p[,2]== "F-"]),3)
t_fc = signif(mean(p_t$shoots.dry.weight[p_t[,2]== "F+"])/mean(p_t$shoots.dry.weight[p_t[,2]== "F-"]),3)
mtext(text = c(paste("Pepper\n(",p_fc,"X)",sep=""),paste("Tomato\n(",t_fc,"X)",sep="")),at = c(1.5,3.5),side = 1,line = -23,cex = 1.5)
text(labels = c("a","b","a","b"),x = c(1.2,2.2,3.2,4.2),y = x$stats[5,]+max(x$stats[5,])/20,cex = 1.5,font =3)

x=boxplot(roots.fresh.weight~treatment*plant,data = productivity_all,ylim = c(0,max(productivity_all$roots.fresh.weight)*1.0),names = names,font = 2,cex.main = 1.5, main = "Root fresh weight",ylab = "grams")
abline(v = 2.5,lty = 2)
#Fold changes
p_fc = signif(mean(p_p$roots.fresh.weight[p_p[,2]== "F+"])/mean(p_p$roots.fresh.weight[p_p[,2]== "F-"]),3)
t_fc = signif(mean(p_t$roots.fresh.weight[p_t[,2]== "F+"])/mean(p_t$roots.fresh.weight[p_t[,2]== "F-"]),3)
mtext(text = c(paste("Pepper\n(",p_fc,"X)",sep=""),paste("Tomato\n(",t_fc,"X)",sep="")),at = c(1.5,3.5),side = 1,line = -23,cex = 1.5)
text(labels = c("a","b","a","b"),x = c(1.2,2.2,3.2,4.2),y = x$stats[5,]+max(x$stats[5,])/20,cex = 1.5,font =3)

x=boxplot(roots.dry.weight~treatment*plant,data = productivity_all,ylim = c(0,max(productivity_all$roots.dry.weight)*1.05),names = names,font = 2,cex.main = 1.5, main = "Root dry weight",ylab = "grams")
abline(v = 2.5,lty = 2)
#Fold changes
p_fc = signif(mean(p_p$roots.dry.weight[p_p[,2]== "F+"])/mean(p_p$roots.dry.weight[p_p[,2]== "F-"]),3)
t_fc = signif(mean(p_t$roots.dry.weight[p_t[,2]== "F+"])/mean(p_t$roots.dry.weight[p_t[,2]== "F-"]),3)
mtext(text = c(paste("Pepper\n(",p_fc,"X)",sep=""),paste("Tomato\n(",t_fc,"X)",sep="")),at = c(1.5,3.5),side = 1,line = -23,cex = 1.5)
text(labels = c("a","b","a","b"),x = c(1.2,2.2,3.2,4.2),y = x$stats[5,]+max(x$stats[5,])/20,cex = 1.5,font =3)

dev.print(device=pdf, "figures/Figure_3_productivity.pdf", onefile=FALSE)
dev.off()

#legend
system("touch figures/legend") 
system("echo 'Figure_2_productivity.pdf : all 6 traits are different (higher) in fertilized plants (Linear mixed effect model (block and replicate are random factors) p<0.0001). Species are different too (p<0.001) except for the fruit fresh weight' >>figures/legend")




####sandbox


mean(p_p[p_p[,2]== "F+",6])/mean(p_p[p_p[,2]== "F-",6])
