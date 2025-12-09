###############################################
### CLEAN RLC (RAPID LIGHT CURVE) PIPELINE  ###
###  For dataset with Unique ID, PAR, ETR   ###
###     Output: RLC parameters in csv.      ###
###############################################

setwd("/Users/riversung/Library/Mobile Documents/com~apple~CloudDocs/River's documents/NTU/Vianney's Lab/R/River_masters_project")
library(tidyr)
library(readr)
library(dplyr)
library(plyr)
library(ggplot2)
library(ggpubr)
library(readxl)
library(minpack.lm) #model for the curves
library(IDPmisc) #for data analyses 
library(nlme) #for the models nonlinear mixed-effects model
library(car)
library(FSA)
library(gridExtra)
library(patchwork)

### ----------------------------------------------------------
### 1. Import your data
### ----------------------------------------------------------

## load FIRST Excel sheet
raw <- read_excel("2025.04_TI_PAM_R.xlsx", sheet = 1)


### ----------------------------------------------------------
### 2. Keep only the columns needed (Unique ID, PAR, Y(II))
### ----------------------------------------------------------

mydata.long <- raw %>%
  select(`Unique ID`, PAR, 'Y(II)') %>%
  dplyr::rename(
    id   = `Unique ID`,
    'Y(II)' = 'Y(II)'
  )
str(mydata.long)

### Convert types
mydata.long$PAR  <- as.numeric(mydata.long$PAR) 
mydata.long$'Y(II)' <- as.numeric(mydata.long$'Y(II)') #convert PAR & Y(II) as numeric
mydata.long$id   <- as.factor(mydata.long$id)
mydata.long = mydata.long[complete.cases(mydata.long), ] #remove rows with no data 



### ----------------------------------------------------------
### 3. Calculate rETR + Plot raw RLC
### ----------------------------------------------------------

mydata.long$rETR = mydata.long$PAR * mydata.long$'Y(II)' #calculate rETR

{
  theme_sleek1 <- function(base_size = 11, base_family = "") {
    half_line <- base_size/2
    theme_light(base_size = base_size, base_family = base_family) +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.length = unit(half_line / 2.2, "pt"),
        strip.background = element_rect(fill = NA, colour = NA),
        strip.text.x = element_text(colour = "grey30"),
        strip.text.y = element_text(colour = "grey30"),
        #strip.placement	 = "inside", 
        strip.switch.pad.grid = unit(2, "cm"), 
        axis.text = element_text(colour = "grey30"),
        axis.title = element_text(colour = "grey30"),
        legend.title = element_text(colour = "grey30", size = rel(0.9)),
        panel.border = element_rect(fill = NA, colour = "grey30", size = 0.5),
        legend.key.size = unit(0.9, "lines"),
        legend.text = element_text(size = rel(0.7), colour = "grey30"),
        legend.key = element_rect(colour = NA, fill = NA),
        legend.background = element_rect(colour = NA, fill = NA),
        legend.position="none",
        plot.title = element_text(colour = "grey30", size = rel(1)),
        plot.subtitle = element_text(colour = "grey30", size = rel(.85))
      )
  }
  theme_set(theme_sleek1())}
p0 = ggplot()+geom_point(mydata.long, mapping = aes(x = PAR, y = rETR), size = 1 )+facet_wrap(~id)
p0

### ----------------------------------------------------------
### 3. Define Platt self-start model (Gerard Ricardo version)
### ----------------------------------------------------------

#adjust data to replace 0 with small positive #s
mydata.long$PAR  <- ifelse(mydata.long$PAR <= 0, 0.1, mydata.long$PAR)
mydata.long$rETR <- ifelse(mydata.long$rETR <= 0, 0.01, mydata.long$rETR)

#source("https://raw.githubusercontent.com/gerard-ricardo/data/master/ssplattmy")  
#for starting values using Platt package author (https://github.com/SWotherspoon/Platt)
{
  SSPlatt.my <- stats::selfStart(
    ~ Pmax*(1-exp(-alpha*I/Pmax))*exp(-beta*I/Pmax),
    function(mCall,data,LHS) {
      ## Extract x and y but do not average replicated x values
      x <- mCall[["I"]]
      y <- LHS
      if (is.language(x) || ((length(x) == 1L) && is.character(x)))
        x <- eval(asOneSidedFormula(x)[[2L]], data)
      x <- as.numeric(x)
      if (is.language(y) || ((length(y) == 1L) && is.character(y)))
        y <- eval(asOneSidedFormula(y)[[2L]], data)
      y <- as.numeric(y)
      keep <- !is.na(x) & !is.na(y)
      x <- x[keep]
      y <- y[keep]
      if(length(unique(x)) < 6) {
        stop("Too few distinct x values to fit a self starting Platt model")
      }
      ord <- order(x)
      x <- x[ord]
      y <- y[ord]
      ## Get initial estimate of smaller rate parameter
      ks <- (length(x)+1)-seq_len(max(4,which(rev(y) >= 0.9*max(y))[1]))
      beta <- -lsfit(x[ks],log(pmax(y[ks]-min(y),1.0E-8)))$coef[2]
      beta <- max(beta,1.0E-8)
      ## Fit partial linear model to estimate second rate parameter
      fit <- nls(y ~ exp(-x*exp(b))-exp(-x*exp(a)),
                 data = list(y=y,x=x,b=log(beta)),
                 start = list(a=log(beta)+3),
                 control=nls.control(maxiter=100,minFactor=1/4096),
                 algorithm = "plinear")
      cf <- coef(fit)
      ## Ensure the difference in exponetials is postive
      cf <- if(cf[2]>0) c(cf[1],log(beta),cf[2:3]) else c(log(beta),cf[1],-cf[2:3])
      cf0 <- cf
      ## Refit partial linear model to estimate both parameters
      fit <- tryCatch(nls(y ~ exp(-x*exp(b))-exp(-x*exp(a)),
                          data = list(y=y,x=x),
                          start = list(a=cf[1],b=cf[2]),
                          control=nls.control(maxiter=100,minFactor=1/4096),
                          algorithm = "plinear"),
                      error=function(e) NULL)
      if(!is.null(fit)) {
        cf1 <- coef(fit)
        ## Ensure the difference in exponetials is postive
        cf1 <- if(cf1[3]>0) cf1 else c(cf1[2:1],-cf1[3:4])
        cf <- if(cf1[1] < cf1[2]) cf else cf1
      }
      setNames(c((exp(cf[1])-exp(cf[2]))*cf[3],exp(cf[2])*cf[3],-cf[3:4]),
               c("alpha","beta","Pmax"))
    },
    c("alpha","beta","Pmax"))
}

#setting up  new  data frame (df) defining log.x values to run
df.x <- data.frame(PAR = seq(0.1, 1500, length = 100))  
vec.x =df.x[,1]

#Find starting values 
starts = mydata.long %>% group_by(id) %>% do(broom::tidy(stats::getInitial(rETR ~ SSPlatt.my(PAR, alpha, beta, Ys), data = . ))) %>% 
  pivot_wider(names_from = names, values_from = x, names_prefix = "") %>% dplyr::select (.,-c('NA'))
colnames(starts) <- c("id", "alpha.s", 'beta.s', 'Pmax.s')
starts = NaRV.omit(starts)


### ----------------------------------------------------------
### 4. Fit model for each ID
### ----------------------------------------------------------

test2 = mydata.long  %>% right_join(.,starts, by = 'id') %>% 
  group_by(id) %>%
  do(model = try(nlsLM(rETR ~ Pmax*(1-exp(-alpha*PAR/Pmax))*(exp(-beta*PAR/Pmax)), 
                       start = list(Pmax = mean(.$Pmax.s),
                                    alpha = mean(.$alpha.s),
                                    beta = mean(.$beta.s)),
                       data = .),silent = TRUE)) 
test2$model[[1]]

### ----------------------------------------------------------
### 5. Predict smooth RLC 
### ----------------------------------------------------------

usq=list()#
for(i in 1:nrow(test2)) {
  out <- try(predict(test2$model[[i]], df.x))
  usq=c(usq,list(out))
}
usq

#put all prediciton in data.frame
df3 = data.frame(t(do.call(rbind.data.frame, usq)), row.names = paste0("", 1:100))  
str(df3)
names = test2$id

#add col names
names <- as.factor(as.character(names)) 
colnames(df3) <- names

#convert all to numeric which adds NAs for error
df3[] <- lapply(df3, function(x) as.numeric(as.character(x)))   
df3$PAR = df.x$PAR

#keep vec.x, add all other columns to factors , add all their values to meas)
df3.long = df3 %>% pivot_longer(-PAR,  names_to = "id" ,values_to = "rETR") %>% data.frame() 
df3.long = dplyr::arrange(df3.long, id) 
str(df3.long)
df3.long$rETR <- as.numeric(as.character(df3.long$rETR)) 

#Building & Displaying plots
p0 = ggplot()+geom_point(mydata.long, mapping = aes(x = PAR, y = rETR), size = 1 )
p0 = p0 + geom_line(df3.long, mapping = aes(x = PAR, y = rETR))
p0 = p0 + facet_wrap(~id)
p0 


### ----------------------------------------------------------
### 6. Extract model parameters
### ----------------------------------------------------------

test = mydata.long  %>% right_join(.,starts, by = 'id') %>% 
  group_by(id) %>%
  do(model = try(broom::tidy(nlsLM(rETR ~ Pmax*(1-exp(-alpha*PAR/Pmax))*(exp(-beta*PAR/Pmax)), 
                                   start = list(Pmax = mean(.$Pmax.s),
                                                alpha = mean(.$alpha.s),
                                                beta = mean(.$beta.s)),
                                   data = .),silent = TRUE)) )  #this get parameters for all models

test$model[[1]]  #check for model 1

#Transform paramters into wide table
unest.test = test %>% unnest(model) 
df.param  = dplyr::select(unest.test, c(id, term, estimate))
dat_wide <- df.param %>% pivot_wider(names_from = term, values_from = estimate)  #%>% dplyr::select(.,-c("NA")) #year goes to columns, their areas go as the values, area is the prefix. Remove the NA cuz no column and give ERROR!

#Calculate derived photophysiological parameters
dat_wide$ETRm = dat_wide$Pmax*(dat_wide$alpha/(dat_wide$alpha+dat_wide$beta))*((dat_wide$beta/(dat_wide$alpha+dat_wide$beta)))^(dat_wide$beta/dat_wide$alpha)
dat_wide$Ek = dat_wide$ETRm/dat_wide$alpha
dat_wide$Em =(dat_wide$Pmax/dat_wide$alpha)*log((dat_wide$alpha+dat_wide$beta)/dat_wide$beta)

id_table <- mydata.long %>% distinct(id) #each sample with data per row 
final.df <- left_join(id_table, dat_wide, by = "id")
final.df <- final.df %>%
  dplyr::rename(`Unique ID` = id)
write.csv(final.df, 'RLC_parameters_TI_20251209.csv')


###long table with repeated data with calculated rETR
#final.df = left_join(dat_wide, mydata.long, by = 'id')
#final.df <- final.df %>% 
  #dplyr::rename(`Unique ID` = id)



