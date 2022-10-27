##################### DATA STRUCTURES #########################################

print("Hello world!")

# create a vector named a with 5 components
a <- c(12, 13, 14, 15) # creating object
a

class(a) # check the class of object a
class("Hello world!")


x <- 1:10 # create a vector of integers
x + 2 # scalar addition
5 * x # scalar multiplication
x^2 # the second power

x
x <- x * 2
x # it is now changed


vec <- rep("a", 12) # create a vector of as, length 12
length(vec) # get the length of the vector
class(vec) # check the class of object vec


x <- c(1, 2, 3, 4)
y <- c(4, 5, 6, 7)

mat <- cbind(x, y) # bind the columns
mat
mat <- rbind(x, y) # bind the rows
mat
t(mat) # transpose of mat
dim(mat) # 4 by 2 matrix
dim(t(mat)) # 2 by 4 matrix


mat2 <- matrix(c(1, 3, 7, 5, -1, 4, 5, 3, 9), nrow = 3)
mat2
mat2[2, ] # get the row 2
mat2[, 2] # get the column 2

# data.frames
a <- c("chr1", "chr2", "chr3", "chr4")
b <- c(250, 410, 200, 450)
c <- c(2000, 4000, 1000, 4000)
d <- c("-", "-", "+", "+")
data <- data.frame(a, b, c, d)
data
colnames(data)

# change column names
colnames(data) <- c("chr", "start", "end", "strand")
data

# another way to create data.frames
data2 <- data.frame(chr = a, start = b, end = c, strand = d)
data2


data2[, 1:3] # columns 2,3,4 of data frame
data2$chr # get the start variablein the data.frame
data2[c(1, 4), ] # get 1st and 4th rows
data2[data2$start > 400, ] # get all rows where start>400


features <- c("promoter", "exon", "intron")
factor(features) # create a factor from features
# Factors are used to store categorical data

# example of a list with 5 components
# a string, a numeric vector, a matrix, and a scalar
l1 <- list(
  name = "Irem",
  numbers = c(1, 2, 3),
  matrix = matrix(1:4, ncol = 2),
  age = 23
)
l1

l1[[2]] # get the second component from the list
l1$age # get the age component from the list
l1[["numbers"]] # component named numbers in list


##############################################################################

############################## FUNCTIONS ######################################
# Functions take an input and perform an operation
numberGenerator <- function(number) {
  print("Your number is:")
  print(sample(1:number, 1))
}
numberGenerator(100)

##############################################################################

########################## LOOPS #############################################

for (i in 1:10) {
  print(i)
}

i <- 1
while (i < 10) {
  print(i)
  i <- i + 1
}

if (!dir.exists("Week-1/data")) {
  dir.create("Week-1/data")
}else{
  message("Directory already exists!")
}
################################ APPLY FAMILY ################################

# The apply family

# sapply
my_function <- function(day) {
  if (day == "Saturday") {
    return("Party time! :)")
  } else if (day == "Monday") {
    return("I HATE mondays! :( :(")
  } else {
    return("No parties today... :(")
  }
}

my_function("Monday")
my_function("Saturday")

days_of_week <- c("Monday","Saturday")
res_vec <- sapply(days_of_week, FUN = function(day) {
  if (day == "Saturday") {
    return("Party time! :)")
  } else if (day == "Monday") {
    return("I HATE mondays! :( :(")
  } else {
    return("No parties today... :(")
  }
})
res_vec

# lapply
res_lst <- lapply(days_of_week, FUN = function(day) {
  if (day == "Saturday") {
    return("Party time! :)")
  } else if (day == "Monday") {
    return("I HATE mondays! :( :(")
  } else {
    return("No parties today... :(")
  }
})
names(res_lst) <- days_of_week
res_lst

# apply (for dataframes/matrices)
my_mat <- matrix(1:10, nrow = 2, ncol = 5)
my_mat

# Apply across all matrix elements
res_mat <- apply(my_mat, MARGIN = 1:2, FUN = function(number) {
  if (number > 5) {
    return(TRUE)
  } else {
    return(FALSE)
  }
})
res_mat

# Apply across columns
res_mat <- apply(my_mat, MARGIN = 2, function(column) {
  return(sum(column))
})
res_mat

# Apply across rows
res_mat <- apply(my_mat, MARGIN = 1, function(row) {
  return(sum(row))
})
res_mat
##############################################################################

############################ PACKAGES ######################################

# install packages
install.packages("dplyr") # insall the package from CRAN
devtools::install_github("BIMSBbioinfo/deconvR") # install package from github
BiocManager::install("GenomicRanges") # install package from Bioconductor

# remove packages
remove.packages("dplyr")

#Call packages to environment
library(tidyverse)

##############################################################################
library(GenomicRanges)

#GRanges
gr <- GRanges(Rle(c("chr2", "chr2", "chr1", "chr3"), c(1, 3, 2, 4)),
                             IRanges::IRanges(1:10, width=10:1))
gr
