library(ggplot2)

CO2_df <- read.csv("CO2_raw_b2.csv")
ID_df <- read.csv("B2_ids.csv")
# N2_df <- read.csv("N2_raw_b2.csv")

CO2_df <- merge(ID_df, CO2_df, by = c("Row"), all = TRUE)

# N1 samples don't have carbon in them!
CO2_df <- CO2_df[!CO2_df$ID %in% c("N1.2A", "N2-1A", "N1-1A", "N2-2A"), ]
CO2_df <- merge(ID_df, CO2_df, by = c("Row", "ID"), all = TRUE)

## Put all row numbers and names back
n_samps <- data.frame(Row = 1:max(CO2_df$Row))
CO2_df <- merge(CO2_df, n_samps, by = c("Row"), all = TRUE)

## Count all the rows to see what runs had less than 3 peaks
CO2_df$count_Row <- 1 
counts <- aggregate(CO2_df$count_Row, by = list(Row = CO2_df$Row), FUN = "sum")

drops <- counts[counts$x < 3,]
CO2_df <- CO2_df[!CO2_df$Row %in% drops$Row,]

## Put all row numbers and names back
n_samps <- data.frame(Row = 1:max(CO2_df$Row))
CO2_df <- merge(ID_df, CO2_df, by = c("Row", "ID"), all = TRUE)


## Get peak for the sample!
CO2_df <- CO2_df[CO2_df$Peak.Nr == 4,]
CO2_df <- merge(ID_df, CO2_df, by = c("Row", "ID"), all = TRUE)


# "cond"       "2.5"        "N1.2A"      "USGS.65.2"  "IVA.Urea.4" "A1"         "A2"         "A3"         "A4"
# "A5"         "A6"         "A7"         "4"          "N2-1A"      "IVA.Urea.2" "A8"         "A9"         "A10"
# "A11"        "A12"        "B1"         "B2"         "B3"         "0.6"        "N1-1A"      "USGS.65.1"  "IVA.Urea.1"
# "B4"         "B5"         "B6"         "B7"         "B8"         "B9"         "B10"        "B11"        "1.5"
# "BR"         "B12"        "C1"         "C2"         "C3"         "C4"         "C5"         "C6"         "1"
# "0.2"        "N2-2A"


# ACA_nms <- c("2.5",  "4" ,   "0.6" , "1.5",  "1", "0.2")

# ACA_df <- CO2_df[CO2_df$ID %in% ACA_nms,]
# print(ACA_df)

## Check to make sure values are reasonably similar for the d13C
ACA_df <- CO2_df[CO2_df$ID %in%  c("2.5",  "4" ,   "0.6" , "1.5",  "1", "0.2"),]
ACA_df$UmolC <- as.numeric(ACA_df$ID)
ACA_df$UmolC <- ACA_df$UmolC*6



lm_fit_ACA <- summary(lm(d.13C.12C~log(Area.44), data = ACA_df ))
Inter_ACA <- lm_fit_ACA$coefficients[1,1]
slope_ACA <- lm_fit_ACA$coefficients[2,1]

plt_fit_ACA <- function(Inter_ACA, slope_ACA, x){
  y = log(x)*slope_ACA + Inter_ACA
  return(y)
}
ggplot() +
  geom_point(data = ACA_df, aes(x = Area.44, y = d.13C.12C)) +
  geom_line(aes(x = 0:26, y =plt_fit_ACA(Inter_ACA = Inter_ACA,
                                         slope_ACA = slope_ACA, 
                                         x = -0:26) ))  + 
  ylim(0,12)

ACA_corr <- function(Area.44,d.13C.12C ){
  lind13C = d.13C.12C -log(Area.44)*slope_ACA
  return(lind13C)
}

CO2_df$lind13C <- ACA_corr(Area.44 = CO2_df$Area.44, d.13C.12C = CO2_df$d.13C.12C)

####################
## Look at the known isotope values

df_stds_C <- data.frame(ID = c( "USGS.65.1", "USGS.65.2",  "IVA.Urea.1",  "IVA.Urea.2", "IVA.Urea.4"), 
                  Knownd13c =c(    -20.29,    -20.29,        -36.54,         -36.54,       -36.54) )


CO2_df <- merge(df_stds_C, CO2_df, by = "ID", all = TRUE)




lm_fit_STDs <- summary(lm(lind13C~Knownd13c, data = CO2_df ))
Inter_STDs <- lm_fit_STDs$coefficients[1,1]
slope_STDs <- lm_fit_STDs$coefficients[2,1]

plt_fit_STDs <- function(Inter_STDs, slope_STDs, x){
  y = slope_STDs*x + Inter_STDs
  return(y)
}

ggplot() +
  geom_point(data = CO2_df, aes(x = Knownd13c, y = lind13C)) +
  geom_line(aes(x = -40:-19, y =plt_fit_STDs(Inter_STDs = Inter_STDs,
                                           slope_STDs = slope_STDs, 
                                           x = -40:-19) ))

## Apply the funtion to correct with 
STD_corr <- function(Inter_STDs,lind13C, slope_STDs ){
  find13c = (lind13C - Inter_STDs)/slope_STDs
  return(find13c)
}

CO2_df$find13c <- STD_corr(Inter_STDs = Inter_STDs,
                          slope_STDs = slope_STDs, 
                          lind13C = CO2_df$lind13C)


## Look at the samples that were run! 
samples <- paste0(c(rep("A", 12), rep("B", 12), rep("C", 6)), 
                 c(1:12, 1:12, 1:6))
CO2_df_smps <- CO2_df[CO2_df$ID %in% samples,]

hist(CO2_df_smps$find13c)
