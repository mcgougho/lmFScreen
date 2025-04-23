# load in libraries
library(haven)
library(dplyr)
library(lmFScreen)

# load in data 2005 - 2006
BMD_2005_2006 <- read_xpt("DXXFEM_2005_2006.xpt")
DEMO_2005_2006 <- read_xpt("DEMO_2005_2006.xpt")
merged_2005_2006 <- merge(BMD_2005_2006, DEMO_2005_2006, by = "SEQN")
# load in data 2007 - 2008
BMD_2007_2008 <- read_xpt("DXXFEM_2007_2008.xpt")
DEMO_2007_2008 <- read_xpt("DEMO_2007_2008.xpt")
merged_2007_2008 <- merge(BMD_2007_2008, DEMO_2007_2008, by = "SEQN")
# load in data 2009 - 2010
BMD_2009_2010 <- read_xpt("DXXFEM_2009_2010.xpt")
DEMO_2009_2010 <- read_xpt("DEMO_2009_2010.xpt")
merged_2009_2010 <- merge(BMD_2009_2010, DEMO_2009_2010, by = "SEQN")
# load in data 2013 - 2014
BMD_2013_2014 <- read_xpt("DXXFEM_2013_2014.xpt")
DEMO_2013_2014 <- read_xpt("DEMO_2013_2014.xpt")
merged_2013_2014 <- merge(BMD_2013_2014, DEMO_2013_2014, by = "SEQN")

# Prepare data for 2005-2006
bmd_2005_2006 <- merged_2005_2006[, c("RIDAGEYR", "DXXNKBMD", "RIAGENDR")]
colnames(bmd_2005_2006) <- c("Age", "Femur_Neck_BMD", "Gender")
bmd_2005_2006 <- na.omit(bmd_2005_2006)
bmd_2005_2006 <- bmd_2005_2006 %>% filter(Age > 30)
# Prepare data for 2007-2008
bmd_2007_2008 <- merged_2007_2008[, c("RIDAGEYR", "DXXNKBMD", "RIAGENDR")]
colnames(bmd_2007_2008) <- c("Age", "Femur_Neck_BMD", "Gender")
bmd_2007_2008 <- na.omit(bmd_2007_2008)
bmd_2007_2008 <- bmd_2007_2008 %>% filter(Age > 30)
# Prepare data for 2009-2010
bmd_2009_2010 <- merged_2009_2010[, c("RIDAGEYR", "DXXNKBMD", "RIAGENDR")]
colnames(bmd_2009_2010) <- c("Age", "Femur_Neck_BMD", "Gender")
bmd_2009_2010 <- na.omit(bmd_2009_2010)
bmd_2009_2010 <- bmd_2009_2010 %>% filter(Age > 30)
# Prepare data for 2013-2014
bmd_2013_2014 <- merged_2013_2014[, c("RIDAGEYR", "DXXNKBMD", "RIAGENDR")]
colnames(bmd_2013_2014) <- c("Age", "Femur_Neck_BMD", "Gender")
bmd_2013_2014 <- na.omit(bmd_2013_2014)
bmd_2013_2014 <- bmd_2013_2014 %>% filter(Age > 30)


# Split data by gender for 2005-2006
bmd_2005_2006_men <- bmd_2005_2006 %>% filter(Gender == 1)
bmd_2005_2006_women <- bmd_2005_2006 %>% filter(Gender == 2)
# Split data by gender for 2007-2008
bmd_2007_2008_men <- bmd_2007_2008 %>% filter(Gender == 1)
bmd_2007_2008_women <- bmd_2007_2008 %>% filter(Gender == 2)
# Split data by gender for 2009-2010
bmd_2009_2010_men <- bmd_2009_2010 %>% filter(Gender == 1)
bmd_2009_2010_women <- bmd_2009_2010 %>% filter(Gender == 2)
# Split data by gender for 2013-2014
bmd_2013_2014_men <- bmd_2013_2014 %>% filter(Gender == 1)
bmd_2013_2014_women <- bmd_2013_2014 %>% filter(Gender == 2)

# Combine the BMD values for men and create a factor for years (cycles)
bmd_men_all <- data.frame(
  BMD = c(bmd_2005_2006_men$Femur_Neck_BMD, bmd_2007_2008_men$Femur_Neck_BMD,
          bmd_2009_2010_men$Femur_Neck_BMD, bmd_2013_2014_men$Femur_Neck_BMD),
  Age = c(bmd_2005_2006_men$Age, bmd_2007_2008_men$Age,
          bmd_2009_2010_men$Age, bmd_2013_2014_men$Age),
  Year = factor(rep(c("2005-2006", "2007-2008", "2009-2010", "2013-2014"),
                    times = c(nrow(bmd_2005_2006_men), nrow(bmd_2007_2008_men),
                              nrow(bmd_2009_2010_men), nrow(bmd_2013_2014_men))),
                levels = c("2005-2006", "2007-2008", "2009-2010", "2013-2014"))
)
# Make into indicator variables
bmd_men_all$Year_2005_2006 <- ifelse(bmd_men_all$Year == "2005-2006", 1, 0)
bmd_men_all$Year_2007_2008 <- ifelse(bmd_men_all$Year == "2007-2008", 1, 0)
bmd_men_all$Year_2009_2010 <- ifelse(bmd_men_all$Year == "2009-2010", 1, 0)
bmd_men_all$Year_2013_2014 <- ifelse(bmd_men_all$Year == "2013-2014", 1, 0)

# testing 2005_2006 vs 2007_2008
set.seed(5678)
# only use 3 out of 4 year columns to avoid rank deficiency of X matrix
# add first two columns together to test coefficient for 2005_2006 vs 2007_2008
X <- cbind(bmd_men_all$Year_2005_2006 + bmd_men_all$Year_2007_2008, bmd_men_all$Year_2007_2008, bmd_men_all$Year_2009_2010, bmd_men_all$Year_2013_2014)
n <- dim(X)[1]
p <- dim(X)[2]
project_out <- cbind(bmd_men_all$Age)
svdP <- svd(project_out, nu = nrow(project_out))
tol <- max(dim(project_out)) * max(svdP$d) * .Machine$double.eps
r <- sum(svdP$d > tol)
U_full <- svdP$u
U_perp <- U_full[, (r+1):ncol(U_full)]
X <- t(U_perp) %*% X
Y <- bmd_men_all$BMD
Y <- t(U_perp) %*% Y
# run lmFScreen on X and Y after projecting out  age column
# testing second column to extract coefficient for 2005_2006 vs 2007_2008
result <- lmFScreen.fit(X, Y, alpha = 0.05, alpha_ov = 0.05, test_cols = 2)
summary(result)


# testing 2005_2006 vs 2009_2010
set.seed(56910)
# add first and third columns together to test coefficient for 2005_2006 vs 2009_2010
X <- cbind(bmd_men_all$Year_2005_2006 + bmd_men_all$Year_2009_2010, bmd_men_all$Year_2007_2008, bmd_men_all$Year_2009_2010, bmd_men_all$Year_2013_2014)
n <- dim(X)[1]
p <- dim(X)[2]
project_out <- cbind(bmd_men_all$Age)
svdP <- svd(project_out, nu = nrow(project_out))
tol <- max(dim(project_out)) * max(svdP$d) * .Machine$double.eps
r <- sum(svdP$d > tol)
U_full <- svdP$u
U_perp <- U_full[, (r+1):ncol(U_full)]
X <- t(U_perp) %*% X
Y <- bmd_men_all$BMD
Y <- t(U_perp) %*% Y
# run lmFScreen on X and Y after projecting out age column
# testing third column to extract coefficient for 2005_2006 vs 2009_2010
result <- lmFScreen.fit(X, Y, alpha = 0.05, alpha_ov = 0.05, test_cols = 3)
summary(result)


# testing 2005_2006 vs 2013_2014
set.seed(561314)
# add first and third columns together to test coefficient for 2005_2006 vs 2013_2014
X <- cbind(bmd_men_all$Year_2005_2006 + bmd_men_all$Year_2013_2014, bmd_men_all$Year_2007_2008, bmd_men_all$Year_2009_2010, bmd_men_all$Year_2013_2014)
n <- dim(X)[1]
p <- dim(X)[2]
project_out <- cbind(bmd_men_all$Age)
svdP <- svd(project_out, nu = nrow(project_out))
tol <- max(dim(project_out)) * max(svdP$d) * .Machine$double.eps
r <- sum(svdP$d > tol)
U_full <- svdP$u
U_perp <- U_full[, (r+1):ncol(U_full)]
X <- t(U_perp) %*% X
Y <- bmd_men_all$BMD
Y <- t(U_perp) %*% Y
# run lmFScreen on X and Y after projecting out age column
# testing fourth column to extract coefficient for 2005_2006 vs 2013_2014
result <- lmFScreen.fit(X, Y, alpha = 0.05, alpha_ov = 0.05, test_cols = 4)
summary(result)


# testing 2007_2008 vs 2009_2010
set.seed(78910)

# add first and second columns together to test coefficient for 2007_2008 vs 2009_2010
X <- cbind(bmd_men_all$Year_2005_2006, bmd_men_all$Year_2007_2008+bmd_men_all$Year_2009_2010, bmd_men_all$Year_2009_2010, bmd_men_all$Year_2013_2014)
n <- dim(X)[1]
p <- dim(X)[2]
project_out <- cbind(bmd_men_all$Age)
svdP <- svd(project_out, nu = nrow(project_out))
tol <- max(dim(project_out)) * max(svdP$d) * .Machine$double.eps
r <- sum(svdP$d > tol)
U_full <- svdP$u
U_perp <- U_full[, (r+1):ncol(U_full)]
X <- t(U_perp) %*% X
Y <- bmd_men_all$BMD
Y <- t(U_perp) %*% Y
# run lmFScreen on X and Y after projecting out age column
# testing third column to extract coefficient for 2007_2008 vs 2009_2010
result <- lmFScreen.fit(X, Y, alpha = 0.05, alpha_ov = 0.05, test_cols = 3)
summary(result)



# testing 2007_2008 vs 2013_2014
set.seed(781314)
# only use 3 out of 4 year columns to avoid rank deficiency of X matrix
# add first and third columns together to test coefficient for 2007_2008 vs 2013_2014
X <- cbind(bmd_men_all$Year_2005_2006, bmd_men_all$Year_2007_2008+bmd_men_all$Year_2013_2014, bmd_men_all$Year_2009_2010, bmd_men_all$Year_2013_2014)
n <- dim(X)[1]
p <- dim(X)[2]
project_out <- cbind(bmd_men_all$Age)
svdP <- svd(project_out, nu = nrow(project_out))
tol <- max(dim(project_out)) * max(svdP$d) * .Machine$double.eps
r <- sum(svdP$d > tol)
U_full <- svdP$u
U_perp <- U_full[, (r+1):ncol(U_full)]
X <- t(U_perp) %*% X
Y <- bmd_men_all$BMD
Y <- t(U_perp) %*% Y
# run lmFScreen on X and Y after projecting out age column
# testing fourth column to extract coefficient for 2007_2008 vs 2013_2014
result <- lmFScreen.fit(X, Y, alpha = 0.05, alpha_ov = 0.05, test_cols = 4)
summary(result)


# testing 2009_2010 vs 2013_2014
set.seed(9101314)
# add second and third columns together to test coefficient for 2009_2014 vs 2013_2014
X <- cbind(bmd_men_all$Year_2005_2006, bmd_men_all$Year_2007_2008, bmd_men_all$Year_2009_2010+bmd_men_all$Year_2013_2014, bmd_men_all$Year_2013_2014)
n <- dim(X)[1]
p <- dim(X)[2]
project_out <- cbind(bmd_men_all$Age)
svdP <- svd(project_out, nu = nrow(project_out))
tol <- max(dim(project_out)) * max(svdP$d) * .Machine$double.eps
r <- sum(svdP$d > tol)
U_full <- svdP$u
U_perp <- U_full[, (r+1):ncol(U_full)]
X <- t(U_perp) %*% X
Y <- bmd_men_all$BMD
Y <- t(U_perp) %*% Y
# run lmFScreen on X and Y after projecting out age column
# testing fourth column to extract coefficient for 2009_2010 vs 2013_2014
result <- lmFScreen.fit(X, Y, alpha = 0.05, alpha_ov = 0.05, test_cols = 4)
summary(result)
