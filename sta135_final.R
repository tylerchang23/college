# PART 1: MULTIPLE LINEAR REGRESSION

setwd('~/Downloads')
data = read.delim('auto-mpg.data', header = FALSE, sep = "")
data = data[ , c(1,3:6)]
colnames(data) = c('Y', 'z1', 'z2', 'z3', 'z4')
data$z2 = as.numeric(data$z2)

Y = data$Y
n = length(Y)
Z = cbind(rep(1,n), data[ ,2:5])
Z = data.matrix(Z)
r = dim(Z)[2]-1

# Summary of the Data
summary(Y)
summary(Z)

library(ggplot2)
# Histogram of MPG (Y)
p = ggplot(data = data, aes(x = Y))+
  geom_histogram(binwidth = 2, color="black", fill="white")+
  geom_vline(aes(xintercept=mean(Y)), color="blue", linetype="dashed", size=1)+
  labs(title="Distribution of Miles Per Gallon (mpg)", x="mpg", y = "Count")

# Histogram of Displacement (z1)
p = ggplot(data = data, aes(x = z1))+
  geom_histogram(binwidth = 75,color="black", fill="white")+
  geom_vline(aes(xintercept=mean(z1)), color="blue", linetype="dashed", size=1)+
  labs(title="Distribution of Displacement", x = "Displacement", y = "Count")

# Histogram of Horsepower (z2)
p = ggplot(data = data, aes(x = z2))+
  geom_histogram(binwidth = 20, color="black", fill="white")+
  geom_vline(aes(xintercept=mean(z2)), color="blue", linetype="dashed", size=1)+
  labs(title="Distribution of Horsepower", x = "Horsepower", y = "Count")

# Histogram of Weight (z3)
p = ggplot(data = data, aes(x = z3))+
  geom_histogram(binwidth = 300, color="black", fill="white")+
  geom_vline(aes(xintercept=mean(z3)), color="blue", linetype="dashed", size=1)+
  labs(title="Distribution of Weight", x = "Weight", y = "Count")

# Histogram of Acceleration (z4)
p = ggplot(data = data, aes(x = z4))+
  geom_histogram(binwidth = 2.5,color="black", fill="white")+
  geom_vline(aes(xintercept=mean(z4)), color="blue", linetype="dashed", size=1)+
  labs(title="Distribution of Acceleration", x = "Acceleration", y = "Count")


# Least Squares Estimates, Beta Hats
beta_hat <- solve(t(Z)%*%Z)%*%t(Z)%*%Y
beta_hat

# R^2 statistic
R_squared <- 1 - sum((Y - Z%*%beta_hat)^2)/sum((Y-mean(Y))^2)
R_squared

# Sigma_Hat Squared 
sigma_hat_square <- sum((Y - Z%*%beta_hat)^2)/(n-r-1)
sigma_hat_square  

# Cov(Beta_hat)
estimated_cov = sigma_hat_square * solve(t(Z)%*%Z)
estimated_cov

alpha <- 0.05
# 95% One-At-A-Time CI for Beta_0 
j <- 0
cat('[',
    beta_hat[j+1] - qt(1-alpha/2, n-r-1)*sqrt(sigma_hat_square * solve(t(Z)%*%Z)[j+1,j+1]),
    ',',
    beta_hat[j+1] + qt(1-alpha/2, n-r-1)*sqrt(sigma_hat_square * solve(t(Z)%*%Z)[j+1,j+1]),
    ']')

# 95% Simultaneous CI for Beta_0 based on CONFIDENCE REGION
j <- 0
cat('[',
    beta_hat[j+1] - sqrt((r+1)*qf(1-alpha, r+1, n-r-1))*sqrt(sigma_hat_square * solve(t(Z)%*%Z)[j+1,j+1]),
    ',',
    beta_hat[j+1] + sqrt((r+1)*qf(1-alpha, r+1, n-r-1))*sqrt(sigma_hat_square * solve(t(Z)%*%Z)[j+1,j+1]),
    ']')

# 95% Simultaneous CI for Beta_0 based on BONFERRONI CORRECTION
j <- 0
cat('[',
    beta_hat[j+1] - qt(1-alpha/(2*(r+1)), n-r-1)*sqrt(sigma_hat_square * solve(t(Z)%*%Z)[j+1,j+1]),
    ',',
    beta_hat[j+1] + qt(1-alpha/(2*(r+1)), n-r-1)*sqrt(sigma_hat_square * solve(t(Z)%*%Z)[j+1,j+1]),
    ']')

# 95% One-At-A-Time CI for Beta_1
j <- 1
cat('[',
    beta_hat[j+1] - qt(1-alpha/2, n-r-1)*sqrt(sigma_hat_square * solve(t(Z)%*%Z)[j+1,j+1]),
    ',',
    beta_hat[j+1] + qt(1-alpha/2, n-r-1)*sqrt(sigma_hat_square * solve(t(Z)%*%Z)[j+1,j+1]),
    ']')

# 95% Simultaneous CI for Beta_1 based on CONFIDENCE REGION
j <- 1
cat('[',
    beta_hat[j+1] - sqrt((r+1)*qf(1-alpha, r+1, n-r-1))*sqrt(sigma_hat_square * solve(t(Z)%*%Z)[j+1,j+1]),
    ',',
    beta_hat[j+1] + sqrt((r+1)*qf(1-alpha, r+1, n-r-1))*sqrt(sigma_hat_square * solve(t(Z)%*%Z)[j+1,j+1]),
    ']')

# 95% Simultaneous CI for Beta_1 based on BONFERRONI CORRECTION
j <- 1
cat('[',
    beta_hat[j+1] - qt(1-alpha/(2*(r+1)), n-r-1)*sqrt(sigma_hat_square * solve(t(Z)%*%Z)[j+1,j+1]),
    ',',
    beta_hat[j+1] + qt(1-alpha/(2*(r+1)), n-r-1)*sqrt(sigma_hat_square * solve(t(Z)%*%Z)[j+1,j+1]),
    ']')


# 95% One-At-A-Time CI for Beta_2 
j <- 2
cat('[',
    beta_hat[j+1] - qt(1-alpha/2, n-r-1)*sqrt(sigma_hat_square * solve(t(Z)%*%Z)[j+1,j+1]),
    ',',
    beta_hat[j+1] + qt(1-alpha/2, n-r-1)*sqrt(sigma_hat_square * solve(t(Z)%*%Z)[j+1,j+1]),
    ']')

# 95% Simultaneous CI for Beta_2 based on CONFIDENCE REGION
j <- 2
cat('[',
    beta_hat[j+1] - sqrt((r+1)*qf(1-alpha, r+1, n-r-1))*sqrt(sigma_hat_square * solve(t(Z)%*%Z)[j+1,j+1]),
    ',',
    beta_hat[j+1] + sqrt((r+1)*qf(1-alpha, r+1, n-r-1))*sqrt(sigma_hat_square * solve(t(Z)%*%Z)[j+1,j+1]),
    ']')

# 95% Simultaneous CI for Beta_2 based on BONFERRONI CORRECTION
j <- 2
cat('[',
    beta_hat[j+1] - qt(1-alpha/(2*(r+1)), n-r-1)*sqrt(sigma_hat_square * solve(t(Z)%*%Z)[j+1,j+1]),
    ',',
    beta_hat[j+1] + qt(1-alpha/(2*(r+1)), n-r-1)*sqrt(sigma_hat_square * solve(t(Z)%*%Z)[j+1,j+1]),
    ']')

# 95% One-At-A-Time CI for Beta_3 
j <- 3
cat('[',
    beta_hat[j+1] - qt(1-alpha/2, n-r-1)*sqrt(sigma_hat_square * solve(t(Z)%*%Z)[j+1,j+1]),
    ',',
    beta_hat[j+1] + qt(1-alpha/2, n-r-1)*sqrt(sigma_hat_square * solve(t(Z)%*%Z)[j+1,j+1]),
    ']')

# 95% Simultaneous CI for Beta_3 based on CONFIDENCE REGION
j <- 3
cat('[',
    beta_hat[j+1] - sqrt((r+1)*qf(1-alpha, r+1, n-r-1))*sqrt(sigma_hat_square * solve(t(Z)%*%Z)[j+1,j+1]),
    ',',
    beta_hat[j+1] + sqrt((r+1)*qf(1-alpha, r+1, n-r-1))*sqrt(sigma_hat_square * solve(t(Z)%*%Z)[j+1,j+1]),
    ']')

# 95% Simultaneous CI for Beta_3 based on BONFERRONI CORRECTION
j <- 3
cat('[',
    beta_hat[j+1] - qt(1-alpha/(2*(r+1)), n-r-1)*sqrt(sigma_hat_square * solve(t(Z)%*%Z)[j+1,j+1]),
    ',',
    beta_hat[j+1] + qt(1-alpha/(2*(r+1)), n-r-1)*sqrt(sigma_hat_square * solve(t(Z)%*%Z)[j+1,j+1]),
    ']')

# 95% One-At-A-Time CI for Beta_4 
j <- 4
cat('[',
    beta_hat[j+1] - qt(1-alpha/2, n-r-1)*sqrt(sigma_hat_square * solve(t(Z)%*%Z)[j+1,j+1]),
    ',',
    beta_hat[j+1] + qt(1-alpha/2, n-r-1)*sqrt(sigma_hat_square * solve(t(Z)%*%Z)[j+1,j+1]),
    ']')

# 95% Simultaneous CI for Beta_4 based on CONFIDENCE REGION
j <- 4
cat('[',
    beta_hat[j+1] - sqrt((r+1)*qf(1-alpha, r+1, n-r-1))*sqrt(sigma_hat_square * solve(t(Z)%*%Z)[j+1,j+1]),
    ',',
    beta_hat[j+1] + sqrt((r+1)*qf(1-alpha, r+1, n-r-1))*sqrt(sigma_hat_square * solve(t(Z)%*%Z)[j+1,j+1]),
    ']')

# 95% Simultaneous CI for Beta_4 based on BONFERRONI CORRECTION
j <- 4
cat('[',
    beta_hat[j+1] - qt(1-alpha/(2*(r+1)), n-r-1)*sqrt(sigma_hat_square * solve(t(Z)%*%Z)[j+1,j+1]),
    ',',
    beta_hat[j+1] + qt(1-alpha/(2*(r+1)), n-r-1)*sqrt(sigma_hat_square * solve(t(Z)%*%Z)[j+1,j+1]),
    ']')
# 95% CI for new Mean Response

# mustang mpg is 20
z_0 = c(1,302,460, 3705,4)

# P5% Prediction Interval for New Response
cat('[',
    z_0%*%beta_hat - sqrt(sigma_hat_square)*sqrt(1+t(z_0)%*%solve(t(Z)%*%Z)%*%z_0)*qt(1-alpha/2, n-r-1),
    ',',
    z_0%*%beta_hat + sqrt(sigma_hat_square)*sqrt(1+t(z_0)%*%solve(t(Z)%*%Z)%*%z_0)*qt(1-alpha/2, n-r-1),
    ']')

# Successfully predicts 2019 Mustang GT

# T Test for each variable

j <- 1
t_stat1 <- (beta_hat[j+1] - 0)/sqrt(sigma_hat_square * solve(t(Z)%*%Z)[j+1,j+1])

j <- 2
t_stat2 <- (beta_hat[j+1] - 0)/sqrt(sigma_hat_square * solve(t(Z)%*%Z)[j+1,j+1])

j <- 3
t_stat3 <- (beta_hat[j+1] - 0)/sqrt(sigma_hat_square * solve(t(Z)%*%Z)[j+1,j+1])

j <- 4
t_stat4 <- (beta_hat[j+1] - 0)/sqrt(sigma_hat_square * solve(t(Z)%*%Z)[j+1,j+1])


alpha <- 0.05
cval_t <- qt(1-alpha/2, n-r-1)
cval_t

# F-test
# H_0: beta_1 = beta_2 = 0

C = rbind(c(0,1,0,0,0),
          c(0,0,1,0,0),
          c(0,0,0,1,0),
          c(0,0,0,0,1))

df_1 <- qr(C)$rank # df_1: rank of matrix R

f_stat <- (t(C%*%beta_hat)%*%solve(C%*%solve(t(Z)%*%Z)%*%t(C))%*%(C%*%beta_hat)/df_1)/sigma_hat_square
f_stat

cval_f <- qf(1-alpha, df_1, n-r-1)
cval_f * (df_1)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# PART 2: TWO SAMPLE TEST, LINEAR DISCRIMINANT ANALYSIS

setwd('~/Downloads')
data = read.csv('StudentsPerformance.csv', header = TRUE, sep = ",")

# TWO SAMPLE TEST

# Trying to test if there's a difference in test scores between test preparation course and no prep
df = data[data$test.preparation.course == "completed" | data$test.preparation.course == 'none', ]
prep = df[df$test.preparation.course == "completed", 6:8]
noPrep = df[df$test.preparation.course == "none", 6:8]

# Summary
table(df$test.preparation.course)
summary(prep)
summary(noPrep)

# Plot the distribution of scores (Math, Writing, Reading)
ggplot(df, aes(x = math.score)) +
  geom_histogram()

ggplot(df, aes(x = writing.score)) +
  geom_histogram()

ggplot(df, aes(x = reading.score)) +
  geom_histogram()

# Plot the distribution of scores by group

p = ggplot(df, aes(x = math.score, fill = test.preparation.course, color = test.preparation.course)) +
  geom_histogram(position = 'identity',binwidth = 5, alpha = 0.3)+
  geom_vline(aes(xintercept=mean(df$math.score[df$test.preparation.course == 'completed'])), color="red", linetype="dashed", size=1)+
  geom_vline(aes(xintercept=mean(df$math.score[df$test.preparation.course == 'none'])), color="blue", linetype="dashed", size=1)+
  labs(title="Math Scores for Test Prep vs. No Prep",x="Math Score", y = "Count")
p + scale_color_brewer(palette="Accent") + 
  theme_minimal()+theme(legend.position="top")

p = ggplot(df, aes(x = writing.score, fill = test.preparation.course, color = test.preparation.course)) +
  geom_histogram(position = 'identity',binwidth = 5, alpha = 0.3)+
  geom_vline(aes(xintercept=mean(df$math.score[df$test.preparation.course == 'completed'])), color="red", linetype="dashed", size=1)+
  geom_vline(aes(xintercept=mean(df$math.score[df$test.preparation.course == 'none'])), color="blue", linetype="dashed", size=1)+
  labs(title="Writing Scores for Test Prep vs. No Prep",x="Writing Score", y = "Count")
p + scale_color_brewer(palette="Accent") + 
  theme_minimal()+theme(legend.position="top")

p = ggplot(df, aes(x = reading.score, fill = test.preparation.course, color = test.preparation.course)) +
  geom_histogram(position = 'identity',binwidth = 5, alpha = 0.3)+
  geom_vline(aes(xintercept=mean(df$math.score[df$test.preparation.course == 'completed'])), color="red", linetype="dashed", size=1)+
  geom_vline(aes(xintercept=mean(df$math.score[df$test.preparation.course == 'none'])), color="blue", linetype="dashed", size=1)+
  labs(title="Reading Scores for Test Prep vs. No Prep",x="Reading Score", y = "Count")
p + scale_color_brewer(palette="Accent") + 
  theme_minimal()+theme(legend.position="top")

# Perform Hotelling T^2 Test
n<-c(358,642)
p<-3

xmean1<-colMeans(prep)
xmean2<-colMeans(noPrep)

d<-xmean1-xmean2

S1<-var(prep)
S2<-var(noPrep)
Sp<-((n[1]-1)*S1+(n[2]-1)*S2)/(sum(n)-2)

t2 <- t(d)%*%solve(sum(1/n)*Sp)%*%d
t2

alpha<-0.05
cval <- (sum(n)-2)*p/(sum(n)-p-1)*qf(1-alpha,p,sum(n)-p-1)
cval

# Simultaneous Confidence Intervals (Difference in Means)
wd<-sqrt(cval*diag(Sp)*sum(1/n))
Cis<-cbind(d-wd,d+wd)
Cis

wd.b<- qt(1-alpha/(2*p),n[1]+n[2]-2) *sqrt(diag(Sp)*sum(1/n))
Cis.b<-cbind(d-wd.b,d+wd.b)
Cis.b

# LINEAR DISCRIMINANT ANALYSIS (LDA)

library(rrcov)
library(MASS)

# 3-D Plot of Test Scores
colors <- c("#E69F00", "#56B4E9")
colors <- colors[as.numeric(df$test.preparation.course)]
s3d = scatterplot3d(df[,6:8], pch = 16, color=colors,grid=TRUE, box=FALSE)
legend("right", legend = levels(df$test.preparation.course),
       col =  c("#E69F00", "#56B4E9"), pch = 16)
title(main = "3-D Plot of Student Test Scores")

# Obtain sample mean vectors for each group
xmean1<-colMeans(prep)
xmean2<-colMeans(noPrep)

# Fisher's Rule:
S.u <- 357*(var(prep)+var(noPrep))/998
w <- solve(S.u)%*%(xmean1-xmean2)
w0 <- -(xmean1+xmean2)%*%w/2

# Use built in LDA function to classify based on training data
lda.obj <- lda(test.preparation.course~math.score+reading.score+writing.score,data=df,prior=c(1,1)/2)
plda <- predict(object=lda.obj,newdata=df)

# Figure out how many students were correctly classified
predicted = data.frame(plda$class)
actual = data.frame(df$test.preparation.course)
combined = cbind.data.frame(actual,predicted)
sum(combined$df.test.preparation.course == combined$plda.class)


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# PART 3: PRINCIPAL COMPONENT ANALYSIS

library(ggplot2)
# Remove labels from built-in Iris
iris.df = data.frame(iris)
iris.df = iris.df[,1:4]

# Plot of Sepal Width vs. Sepal Length
ggplot(iris, aes(x=Sepal.Width, y=Sepal.Length, color=Species)) +
  geom_point()+
  labs(title="Sepal Width vs. Sepal Length", x = "Sepal Width", y = "Sepal Length")

# Obtain principal components
pca_results = data.frame(prcomp(iris.df)$x)
pca_results$species = iris$Species

# Plot the first 2 principal components
ggplot(pca_results, aes(x=PC1, y=PC2, color=species))+
  geom_point()+
  labs(title="First Two Principal Components", x = "First Principal Component", y = "Second Principal Compo")

# Obtain proportion of variance/loadings
iris.pc <- princomp(iris.df,cor=TRUE)
summary(iris.pc, loadings = TRUE)

# Scree Plot
plot(1:(length(iris.pc$sdev)),  (iris.pc$sdev)^2, type='b', 
     main="Scree Plot", xlab="Number of Components", ylab="Eigenvalue Size")

# Plotting the PC scores for the sample data in the space of the first two principal components:
par(pty="s")
plot(iris.pc$scores[,1], iris.pc$scores[,2], 
     xlab="PC 1", ylab="PC 2", type ='n', lwd=2)
text(iris.pc$scores[,1], iris.pc$scores[,2], labels = colnames(iris.df), cex=0.7, lwd=2)
biplot(iris.pc)

