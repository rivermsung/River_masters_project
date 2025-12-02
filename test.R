# Set the working directory
setwd("/Users/riversung/Library/Mobile Documents/com~apple~CloudDocs/River's documents/NTU/Vianney's Lab/PAM/For R")

# Load necessary libraries
library(tidyr)
library(readr)
library(plyr)
library(dplyr)

# Define a function to read in csv files and add a "disc" column with the filename
read_csv_filename1 <- function(filename){
  ret <- read_csv(filename)
  ret$disc <- filename
  ret
}

# Read in all csv files in the working directory
filenames <- list.files(full.names = TRUE)
data1 <- ldply(filenames, read_csv_filename1)

# Extract the ID from the start of the file name and add it as a new column
data1$disc1 <- data1 %>% separate(disc, c("a", "b", 'c', 'd')) %>% .$b

# Read in environmental factors data
#env.fact <- read.table(file="https://raw.githubusercontent.com/gerard-ricardo/data/master/postset%20treat%20mil", 
                       #header = TRUE,
                       #dec = ",", 
                       #na.strings = c("",".","NA"))

# Join data1 and env.fact by the disc1 column
#data2 <- left_join(data1, env.fact, by = 'disc1')

# Select only necessary columns from data1
data.s <- dplyr::select(data1, c('disc1', 'PAR', 'Y (II)'))

# Create long format and pivot data
data1.long <- data.s %>% pivot_longer(-c(disc1, PAR), 
                                      names_to = "rep" ,
                                      values_to = "meas")

# Create individual ID for each RLC
data1.long$id <- paste(data1.long$disc1, data1.long$rep, sep = "")


# Calculate rETR 
data1.long$rETR <- data1.long$PAR * data1.long$meas

# Load ggplot2 library
library(ggplot2)

# source theme 
source("https://raw.githubusercontent.com/gerard-ricardo/data/master/theme_sleek1")

# create a scatter plot with ggplot2
p0 <- ggplot() + 
  geom_point(data = data1.long, 
             mapping = aes(x = PAR, y = rETR), 
             size = 1) + 
  facet_wrap(~id)

# print the plot
print(p0)

data1.l = data1.long[which(data1.long$rETR<100),] #remove rTER anomalies
data1.l$PAR <- ifelse(data1.l$PAR <= 0, 0.1, data1.l$PAR) #
data1.l$rETR <- ifelse(data1.l$rETR <= 0, 0.01, data1.l$rETR)


####test 1####
data1.s = split(data1.l, data1.l$id)
data1.s$`LQ53Y (II)`
source("https://raw.githubusercontent.com/gerard-ricardo/data/master/ssplattmy") #this is the starter values. Full #credit to the author of the Platt package
library(minpack.lm)
start = unname(getInitial(rETR ~ SSPlatt.my(PAR, alpha, beta, Pmax), data1.s$`LQ53Y (II)`)) #this finds the starting #values
md1 = nlsLM(rETR ~ Pmax*(1-exp(-alpha*PAR/Pmax))*(exp(-beta*PAR/Pmax)), start=list(Pmax=start[3],alpha=start[1], beta=start[2]), data = data1.s$LQ53Y (II)) #notice I added the starting values in the equation
df.x <- data.frame(PAR = seq(0.1, 926, length = 100)) #setting up new data frame (df) defining log.x values to run
vec.x =df.x[,1]
plot(data1.s$LQ531:Y (II)$PAR, data1.s$LQ531:Y (II)$rETR, col = 'red')
lines(vec.x, predict(md1, df.x)) #looks good for m001Y(II)1

summary(data1.s$`LQ53Y (II)`)
plot(`LQ53Y (II)`$PAR, `LQ53Y (II)`$rETR)


####test 2####
data1.s$`TI016Y (II)`
source("https://raw.githubusercontent.com/gerard-ricardo/data/master/ssplattmy") #this is the starter values. Full #credit to the author of the Platt package
library(minpack.lm)
start = unname(getInitial(rETR ~ SSPlatt.my(PAR, alpha, beta, Pmax), data1.s$`TI016Y (II)`)) #this finds the starting #values
md1 = nlsLM(rETR ~ Pmax*(1-exp(-alpha*PAR/Pmax))*(exp(-beta*PAR/Pmax)), start=list(Pmax=start[3],alpha=start[1], beta=start[2]), data = data1.s$LQ53Y (II)) #notice I added the starting values in the equation
df.x <- data.frame(PAR = seq(0.1, 926, length = 100)) #setting up new data frame (df) defining log.x values to run
vec.x =df.x[,1]
plot(data1.s$LQ531:Y (II)$PAR, data1.s$LQ531:Y (II)$rETR, col = 'red')
lines(vec.x, predict(md1, df.x)) #looks good for m001Y(II)1

summary(data1.s$`LQ53Y (II)`)
plot(`LQ53Y (II)`$PAR, `LQ53Y (II)`$rETR)




