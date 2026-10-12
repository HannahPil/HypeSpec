library(lme4)

# Import Phenotypic Data that was filtered to just have genotype with with multiple replicates
laData = read.table(file  ="LVs_Reflectance_Data_MAIZE_GENOID_10_20_2021__21.csv",header = TRUE, na.strings = "NA",skipNul = TRUE, sep=",")

# View the head of phenotypic file
head(laData)

#Filter out unrepicated data
Tally = table(laData[,1])
repTally = Tally[!Tally == 1]
repGeno = names(repTally)
laData = laData[laData[,1] %in% repGeno,]

# Set genotypes as a character type
genotypes = as.character(laData[,1])


# Count number of columns and set the columns with numbers as numeric, excluding genotype column that is indexed at 1
colCount = ncol(laData)
laData[,2:colCount] = sapply(laData[,2:colCount],as.numeric)


############################################################################################
#For loop to calculate heritability of all columns that has a phenotypic value for genotype#
############################################################################################

# Range of all columns containing phenotypic data
my_range = 2:colCount

# Create empty 'y' list to store heritability values and empty 'pheno' list to store column names that is being used as the phenotypic value for each genotype
y = NULL
pheno = NULL

# Get column names for dataframe
names = colnames(laData)

for (i in my_range)
{
  Ang = laData[,i]
  
 # Linear Model
 model = lmer(Ang~(1|genotypes),data = laData)
 summary(model)


# Extracting variance components
 Variances = as.data.frame(VarCorr(model))
 Variances

 GenoVariance = Variances[1,4]
 GenoVariance
 ResidVariance = Variances[2,4]
 ResidVariance


# Calculating Heritability

 H2 = GenoVariance/(GenoVariance+(ResidVariance/2))
 H2 

 
 y =  rbind(y, H2) # List of all heritability measurements for specific phenotype, list appends heritability measurements to list as it is being measured
 pheno = rbind(pheno,names[i]) # List of phenotypes being measured, list appends the phentoype to this list as it is being measured


 heriData = data.frame(pheno,y) # Dataframe with all different phenotypes and their heritability measurements


}

# Save dataframe as a csv file
write.csv(heriData, file = "HeriData_Maize2020.csv")
