library(ggplot2)
library(ggbiplot)
library(dplyr)
library(lubridate)
library(caTools)
library("depmixS4")
library("fpp2")
library(corrplot)

# Set working directory to source file location 
setwd("F:/SFU/CMPT318/Term Project")

df <- read.table(file="TermProjectData.txt", header=T, sep=",")

df$Day <-  as.POSIXlt(df$Date)$wday
df <- subset(df, is.na(Global_active_power) == FALSE
                     & is.na(Global_reactive_power) == FALSE
                     & is.na(Voltage) == FALSE
                     & is.na(Global_intensity) == FALSE
                     & is.na(Sub_metering_1) == FALSE
                     & is.na(Sub_metering_2) == FALSE
                     & is.na(Sub_metering_3) == FALSE)

df <- subset(df, (df$Time >= "12:00:00") & (df$Time <= "15:59:59"))
#Part 1 Train and Test HMM
df <- df %>% 
  mutate(Date = lubridate::parse_date_time(Date,"%d/%m/%Y")) %>% 
  arrange( Date )

# Find corr
cor_df <- cor(df[,c(3,4,5,6,7,8,9)], use="pairwise.complete.obs")
corrplot(cor_df)

# Plot values
mean_GAP<- aggregate(list(Global_active_power = df$Global_active_power), by=list(Time = df$Time), mean)
mean_GAP$Time <- strptime(mean_GAP$Time, "%H:%M:%S")
plot(mean_GAP$Time, mean_GAP$Global_active_power,xlab="Time",ylab = 'Global active power' )
title('Average Global active power')

mean_GI<- aggregate(list(Global_intensity = df$Global_intensity), by=list(Time = df$Time), mean)
mean_GI$Time <- strptime(mean_GI$Time, "%H:%M:%S")
plot(mean_GI$Time, mean_GI$Global_intensity,xlab="Time",ylab = 'Global intensity' )
title('Average Global intensity')

mean_SM<- aggregate(list(Sub_metering_3 = df$Sub_metering_3), by=list(Time = df$Time), mean)
mean_SM$Time <- strptime(mean_SM$Time, "%H:%M:%S")
plot(mean_SM$Time, mean_SM$Sub_metering_3,xlab="Time",ylab = 'Sub metering 3' )
title('Average Sub metering 3')



# ----------------------------------------------------
# 1
df_initial <- df

df_weekly_avg <- df_initial %>%
  mutate(Date = floor_date(Date, "week")) %>%
  group_by(Date) %>%
  summarize(Global_active_power = mean(Global_active_power),
            Global_reactive_power = mean(Global_reactive_power),
            Voltage = mean(Voltage),
            Global_intensity = mean(Global_intensity),
            Sub_metering_1 = mean(Sub_metering_1),
            Sub_metering_2 = mean(Sub_metering_2),
            Sub_metering_3 = mean(Sub_metering_3))

pca <- prcomp(df_weekly_avg[,c(2,3,4,5,6,7,8)], scale=TRUE)
ggbiplot(pca,choices=c(1,2), obs.scale = 1, var.scale = 1,
         ellipse = TRUE, circle = TRUE) +
  scale_color_discrete(name = '') +
  theme(legend.direction = 'horizontal', legend.position = 'top')
summary(pca)
pca



df_weekly_avg$Global_active_power <- as.numeric(df_weekly_avg$Global_active_power)
df_weekly_avg$Sub_metering_3 <- as.numeric(df_weekly_avg$Sub_metering_3)
df_weekly_avg$Global_intensity <- as.numeric(df_weekly_avg$Global_intensity)

# Set training and test data
train_data <- subset(df_weekly_avg, (as.POSIXlt(df_weekly_avg$Date,format="%d/%m/%Y") < "2009-1-1"))
test_data <- subset(df_weekly_avg, (as.POSIXlt(df_weekly_avg$Date,format="%d/%m/%Y") >= "2009-1-1"))

# Train Uni and multi
set.seed(22)
model_uni_n5_train <- depmix(response = Global_active_power ~ 1, data = train_data, nstates = 5, ntimes = nrow(train_data))
fm_uni_n5 <- fit(model_uni_n5_train)
print(fm_uni_n5)


#'log Lik.' -8.060985 (df=34)
#'AIC:  84.12197 
#'BIC:  174.9981

model_multi_n5_train <- depmix(response = list(Global_intensity ~ 1, Global_active_power ~ 1, Sub_metering_3 ~ 1 ), data = train_data, nstates = 5, family=list(gaussian(),gaussian(), gaussian()),ntimes = nrow(train_data))
fm_multi_n5 <- fit(model_multi_n5_train)
print(fm_multi_n5)
#Best Result n=5 :
#'log Lik.' -324.3738 (df=54)
#AIC:  756.7475  
#BIC:  901.0803 

#Test Uni and multi

model_uni_n5_test <- depmix(response = Global_active_power ~ 1, data = test_data, nstates = 5, ntimes = nrow(test_data))
mod_uni_test <- setpars(model_uni_n5_test,getpars(fm_uni_n5))
fm_uni_test <- fit(mod_uni_test)
# converged at iteration 71 with logLik: 0.6713427 

logLik(fm_uni_n5) / nrow(train_data)
logLik(fm_uni_test) / nrow(test_data)


model_multi_n5_test <- depmix(response = list(Global_intensity ~ 1, Global_active_power ~ 1, Sub_metering_3 ~ 1 ), data = test_data, nstates = 5, family=list(gaussian(),gaussian(), gaussian()),ntimes = nrow(test_data))
mod_multi_test <- setpars(model_multi_n5_test,getpars(fm_multi_n5))
fm_multi_test <- fit(mod_multi_test)
# converged at iteration 44 with logLik: -142.4851 

logLik(fm_multi_n5) / nrow(train_data)
logLik(fm_multi_test) / nrow(test_data)


# =====================================================


# =====================================================
# No 2
# Data anom 1
df_anom_1 <- read.table(file="Data1(WithAnomalies).txt", header=T, sep=",")
df_anom_1$Day <-  as.POSIXlt(df_anom_1$Date)$wday
df_anom_1 <- subset(df_anom_1, is.na(Global_active_power) == FALSE
             & is.na(Global_reactive_power) == FALSE
             & is.na(Voltage) == FALSE
             & is.na(Global_intensity) == FALSE
             & is.na(Sub_metering_1) == FALSE
             & is.na(Sub_metering_2) == FALSE
             & is.na(Sub_metering_3) == FALSE)

df_anom_1 <- subset(df_anom_1, (df_anom_1$Time >= "12:00:00") & (df_anom_1$Time <= "15:59:59"))

df_anom_1 <- df_anom_1 %>% 
  mutate(Date = lubridate::parse_date_time(Date,"%d/%m/%Y")) %>% 
  arrange( Date )


df_weekly_avg <- df_anom_1 %>%
  mutate(Date = floor_date(Date, "week")) %>%
  group_by(Date) %>%
  summarize(Global_active_power = mean(Global_active_power),
            Global_reactive_power = mean(Global_reactive_power),
            Voltage = mean(Voltage),
            Global_intensity = mean(Global_intensity),
            Sub_metering_1 = mean(Sub_metering_1),
            Sub_metering_2 = mean(Sub_metering_2),
            Sub_metering_3 = mean(Sub_metering_3))

df_init_avg <- df_initial %>%
  mutate(Date = floor_date(Date, "week")) %>%
  group_by(Date) %>%
  summarize(Global_active_power = mean(Global_active_power),
            Global_reactive_power = mean(Global_reactive_power),
            Voltage = mean(Voltage),
            Global_intensity = mean(Global_intensity),
            Sub_metering_1 = mean(Sub_metering_1),
            Sub_metering_2 = mean(Sub_metering_2),
            Sub_metering_3 = mean(Sub_metering_3))


df_weekly_avg$Global_active_power <- as.numeric(df_weekly_avg$Global_active_power)
df_weekly_avg$Sub_metering_3 <- as.numeric(df_weekly_avg$Sub_metering_3)
df_weekly_avg$Global_intensity <- as.numeric(df_weekly_avg$Global_intensity)

# Test Uni
model_uni_n5_test <- depmix(response = Global_active_power ~ 1, data = df_weekly_avg, nstates = 5, ntimes = nrow(df_weekly_avg))
mod_uni_test <- setpars(model_uni_n5_test,getpars(fm_uni_n5))
fm_uni_test <- fit(mod_uni_test)

logLik(fm_uni_n5) / nrow(train_data)
logLik(fm_uni_test) / nrow(df_weekly_avg)

# Test Multi
model_multi_n5_test <- depmix(response = list(Global_intensity ~ 1, Global_active_power ~ 1, Sub_metering_3 ~ 1 ), data = df_weekly_avg, nstates = 5, family=list(gaussian(),gaussian(), gaussian()),ntimes = nrow(df_weekly_avg))
mod_multi_test <- setpars(model_multi_n5_test,getpars(fm_multi_n5))
fm_multi_test <- fit(mod_multi_test)

logLik(fm_multi_n5) / nrow(train_data)
logLik(fm_multi_test) / nrow(df_weekly_avg)

ggplot() + 
  geom_line(data=df_weekly_avg, aes(x=Date, y=Global_active_power), color='blue')
ggplot() +
  geom_line(data=train_data, aes(x=Date, y=Global_active_power), color='black')


# -------------------------------------------------------------

# Data Anom 2
df_anom_2 <- read.table(file="Data2(WithAnomalies).txt", header=T, sep=",")
df_anom_2$Day <-  as.POSIXlt(df_anom_2$Date)$wday
df_anom_2 <- subset(df_anom_2, is.na(Global_active_power) == FALSE
                    & is.na(Global_reactive_power) == FALSE
                    & is.na(Voltage) == FALSE
                    & is.na(Global_intensity) == FALSE
                    & is.na(Sub_metering_1) == FALSE
                    & is.na(Sub_metering_2) == FALSE
                    & is.na(Sub_metering_3) == FALSE)

df_anom_2 <- subset(df_anom_2, (df_anom_2$Time >= "12:00:00") & (df_anom_2$Time <= "15:59:59"))

df_anom_2 <- df_anom_2 %>% 
  mutate(Date = lubridate::parse_date_time(Date,"%d/%m/%Y")) %>% 
  arrange( Date )


df_weekly_avg <- df_anom_2 %>%
  mutate(Date = floor_date(Date, "week")) %>%
  group_by(Date) %>%
  summarize(Global_active_power = mean(Global_active_power),
            Global_reactive_power = mean(Global_reactive_power),
            Voltage = mean(Voltage),
            Global_intensity = mean(Global_intensity),
            Sub_metering_1 = mean(Sub_metering_1),
            Sub_metering_2 = mean(Sub_metering_2),
            Sub_metering_3 = mean(Sub_metering_3))

df_weekly_avg$Global_active_power <- as.numeric(df_weekly_avg$Global_active_power)
df_weekly_avg$Sub_metering_3 <- as.numeric(df_weekly_avg$Sub_metering_3)
df_weekly_avg$Global_intensity <- as.numeric(df_weekly_avg$Global_intensity)


# Test Uni
model_uni_n5_test <- depmix(response = Global_active_power ~ 1, data = df_weekly_avg, nstates = 5, ntimes = nrow(df_weekly_avg))
mod_uni_test <- setpars(model_uni_n5_test,getpars(fm_uni_n5))
fm_uni_test <- fit(mod_uni_test)

logLik(fm_uni_n5) / nrow(train_data)
logLik(fm_uni_test) / nrow(df_weekly_avg)

# Test Multi
model_multi_n5_test <- depmix(response = list(Global_intensity ~ 1, Global_active_power ~ 1, Sub_metering_3 ~ 1 ), data = df_weekly_avg, nstates = 5, family=list(gaussian(),gaussian(), gaussian()),ntimes = nrow(df_weekly_avg))
mod_multi_test <- setpars(model_multi_n5_test,getpars(fm_multi_n5))
fm_multi_test <- fit(mod_multi_test)

logLik(fm_multi_n5) / nrow(train_data)
logLik(fm_multi_test) / nrow(df_weekly_avg)

# -------------------------------------------------------------

# Data anom 3
df_anom_3 <- read.table(file="Data3(WithAnomalies).txt", header=T, sep=",")
df_anom_3$Day <-  as.POSIXlt(df_anom_3$Date)$wday
df_anom_3 <- subset(df_anom_3, is.na(Global_active_power) == FALSE
                    & is.na(Global_reactive_power) == FALSE
                    & is.na(Voltage) == FALSE
                    & is.na(Global_intensity) == FALSE
                    & is.na(Sub_metering_1) == FALSE
                    & is.na(Sub_metering_2) == FALSE
                    & is.na(Sub_metering_3) == FALSE)

df_anom_3 <- subset(df_anom_3, (df_anom_3$Time >= "12:00:00") & (df_anom_3$Time <= "15:59:59"))

df_anom_3 <- df_anom_3 %>% 
  mutate(Date = lubridate::parse_date_time(Date,"%d/%m/%Y")) %>% 
  arrange( Date )


df_weekly_avg <- df_anom_3 %>%
  mutate(Date = floor_date(Date, "week")) %>%
  group_by(Date) %>%
  summarize(Global_active_power = mean(Global_active_power),
            Global_reactive_power = mean(Global_reactive_power),
            Voltage = mean(Voltage),
            Global_intensity = mean(Global_intensity),
            Sub_metering_1 = mean(Sub_metering_1),
            Sub_metering_2 = mean(Sub_metering_2),
            Sub_metering_3 = mean(Sub_metering_3))



# Test Uni
model_uni_n5_test <- depmix(response = Global_active_power ~ 1, data = df_weekly_avg, nstates = 5, ntimes = nrow(df_weekly_avg))
mod_uni_test <- setpars(model_uni_n5_test,getpars(fm_uni_n5))
fm_uni_test <- fit(mod_uni_test)

logLik(fm_uni_n5) / nrow(train_data)
logLik(fm_uni_test) / nrow(df_weekly_avg)

# Test Multi
model_multi_n5_test <- depmix(response = list(Global_intensity ~ 1, Global_active_power ~ 1, Sub_metering_3 ~ 1 ), data = df_weekly_avg, nstates = 5, family=list(gaussian(),gaussian(), gaussian()),ntimes = nrow(df_weekly_avg))
mod_multi_test <- setpars(model_multi_n5_test,getpars(fm_multi_n5))
fm_multi_test <- fit(mod_multi_test)


logLik(fm_multi_n5) / nrow(train_data)
logLik(fm_multi_test) / nrow(df_weekly_avg)


