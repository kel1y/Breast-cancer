options(digits = 3)
library(matrixStats)
library(tidyverse)
library(caret)
library(dslabs)
data(brca)

# how many sample are in the dataset

NROW(brca$y)
NROW(brca$x)

# how many predictors are in the matrix

NCOL(brca$x)

# proportion of malgninant

prop <- sum(brca$y == "M") / length(brca$y)

# column with high mean

col_means <- colMeans(brca$x)
max_mean <- which.max(col_means)

# which column with low standard deviation

col_sd <- colSds(brca$x)
min <- which.min(col_sd)

# sweep the column two times, subtract column means then divide column standard deviation

x_centerd <- sweep(brca$x, 2, colMeans(brca$x))
x_scaled <- sweep(x_centerd, 2, colSds(brca$x), FUN = "/")

first_coumn_sd <- sd(x_scaled[,1])
first_coumn_sd <- median(x_scaled[,1])

# proportion of variance explained by first component 
# and number of component that can represent 90% of variance

pca_result <- prcomp(x_scaled)
pca_var <- cumsum(pca_result$sdev^2) / sum(pca_result$sdev^2)
number_of_component <- sum(pca_var <= 0.9) + 1

# Plot the first two principal components with color representing tumor type
# Malignant tumors tend to have larger values of PC1 than benign tumors.

data.frame(pca_result$x[,1:2], type =  brca$y) %>% ggplot(aes(PC1,type)) + geom_point()
data.frame(pca_result$x[,1:2], type =  brca$y) %>% ggplot(aes(PC2,type)) + geom_point()

# Make a boxplot of the first 10 PCs grouped by tumor type
# Which PCs are significantly different enough by tumor type that there is no overlap in the interquartile ranges (IQRs) for benign and malignant samples?

data.frame(type = brca$y, pca_result[,1:10]) %>% gather(key="PC", value="value", -type) %>% ggplot(aes(PC, value, fill=type))

# scale brca data into test set and trainning set

set.seed(1) 
test_index <- createDataPartition(brca$y, times = 1, p = 0.2, list = FALSE)
test_x <- x_scaled[test_index,]
test_y <- brca$y[test_index]
train_x <- x_scaled[-test_index,]
train_y <- brca$y[-test_index]

# proportion of benign in training and test set

prop_b_test <- sum(test_y == 'B') / length(test_y)
prop_b_train <- sum(train_y == 'B') / length(train_y)

# train with logistic regression using caret 

train_glm <- train(train_x, train_y, method="glm")
y_hat_glm <- predict(train_glm, test_x)
mean(y_hat_glm == test_y) #accuracy

# train with loess using caret 

set.seed(5)
train_loess <- train(train_x, train_y, method = "gamLoess")
y_hat_loess <- predict(train_loess, test_x)
mean(y_hat_loess == test_y) #accuracy

# train with knn using caret

set.seed(7)
train_knn <- train(train_x, train_y, method = 'knn', tuneGrid = data.frame(k=seq(3,21)))
train_knn$bestTune #utilized K
y_hat_knn <- predict(train_knn, test_x)
mean(y_hat_knn == test_y)

# train with random forest using caret while keeping importances

set.seed(9)
grid <- data.frame(mtry = c(3,5,7,9))
train_randomforest <- train(train_x, train_y, method = 'rf', tuneGrid = grid, importance = TRUE)
y_hat_randomforest <- predict(train_randomforest, test_x)
mean(y_hat_randomforest == test_y)
varImp(train_randomforest)

# creating ensemble of the above data starting from logistic regression to the end

ensemble <- cbind(glm = y_hat_glm == "B", loess = y_hat_loess == 'B', rf = y_hat_randomforest == 'B', knn = y_hat_knn == 'B')
ensemble_hats <- ifelse(rowMeans(ensemble) > 0.5, "B", "M")
mean(ensemble_hats == test_y)