# install.packages("pROC")
library(ggplot2)
library(pROC)

pheno = read.table('EA23k_plink.txt', header=T)

staph = read.table('StaphAureusBloodCulture_GRIDs.txt', header=T)

pheno.1 = pheno[,colnames(pheno) %in% c("FID", "IID", "AGE_in_days", "ARRAY", "GENDER", "PC1", "PC2", "PC3", "PC4", "PC5", "X41.1")]

pheno.1$STAPH = pheno.1$FID %in% staph[,1]
pheno.1 = subset(pheno.1, pheno.1$X41.1 >= 0)

# logistic model
x = lm(STAPH ~ AGE_in_days + GENDER + PC1 + PC2 + PC3 + PC4 + PC5 , data=pheno.1)

# phecode classifier
plot(roc(pheno.1$STAPH, pheno.1$X41.1, direction="<"),  auc.polygon.col="lightblue",
     col="red", lwd=3, main="", cex.lab=1.8)

auc(roc(pheno.1$STAPH, pheno.1$X41.1, direction="<"))
# Setting levels: control = FALSE, case = TRUE
# Area under the curve: 0.9379

# PC1
par(new=T)
plot(roc(pheno.1$STAPH, pheno.1$PC1, direction="<"), 
     col="orange", lwd=3, main="", cex.lab=1.8)
auc(roc(pheno.1$STAPH, pheno.1$PC1, direction="<"))
# Setting levels: control = FALSE, case = TRUE
# Area under the curve: 0.514

# Fitted values
par(new=T)
plot(roc(pheno.1$STAPH, x$fitted.values, direction="<"),  auc.polygon.col="lightblue",
     col="blue", lwd=3, main="", cex.lab=1.8)

auc(roc(pheno.1$STAPH, x$fitted.values, direction="<"))
# Setting levels: control = FALSE, case = TRUE
# Area under the curve: 0.5683

legend("bottomright", legend=c("Phecode (041.1)", "PC1", "Age+Sex+5PCs"),
       col=c("red", "orange", "blue"), lwd=2)

# Bootstrap
for (i in 1:100)
{

	pheno.2 = pheno.1[sample(nrow(pheno.1), nrow(pheno.1)),]
	tmp = auc(roc(pheno.2$STAPH, pheno.2$X41.1, direction="<"))
	print(tmp)
}
