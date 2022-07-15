## Read in the files - make an ID file with the row names and actual sample names!
library(ggplot2)

ID_df <- read.csv("B2_ids.csv")
N2_df <- read.csv("N2_raw_b2.csv")

## Our sample was the thrid peak!
N2_df <- N2_df[N2_df$Peak.Nr == 3,]

## Fill in the blanks
n_samps <- data.frame(Row = 1:max(N2_df$Row))
N2_df <- merge(N2_df, n_samps, by = c("Row"), all = TRUE)

## Name things
N2_df <- merge(ID_df, N2_df, by = "Row", all = TRUE)



# "cond"       "2.5"        "N1.2A"      "USGS.65.2"  "IVA.Urea.4" "A1"         "A2"         "A3"         "A4"
# "A5"         "A6"         "A7"         "4"          "N2-1A"      "IVA.Urea.2" "A8"         "A9"         "A10"
# "A11"        "A12"        "B1"         "B2"         "B3"         "0.6"        "N1-1A"      "USGS.65.1"  "IVA.Urea.1"
# "B4"         "B5"         "B6"         "B7"         "B8"         "B9"         "B10"        "B11"        "1.5"
# "BR"         "B12"        "C1"         "C2"         "C3"         "C4"         "C5"         "C6"         "1"
# "0.2"        "N2-2A"


## Look to see if the values are consistant, or if something went weird with these samples 
ACA_nms <- c("2.5",  "4" ,   "0.6" , "1.5",  "1", "0.2")

ACA_df <- N2_df[N2_df$ID %in% ACA_nms,]
print(ACA_df)

ACA_df <- N2_df[N2_df$ID %in%  c("2.5",  "4" ,   "0.6" , "1.5",  "1"),]


## Fit a line to the standard:
lm_fit_ACA <- summary(lm(d.15N.14N~log(Area.28), data = ACA_df ))
Inter_ACA <- lm_fit_ACA$coefficients[1,1]
slope_ACA <- lm_fit_ACA$coefficients[2,1]

plt_fit_ACA <- function(Inter_ACA, slope_ACA, x){
  y = log(x)*slope_ACA + Inter_ACA
  return(y)
}
ggplot() +
  geom_point(data = ACA_df, aes(x = Area.28, y = d.15N.14N)) +
  geom_line(aes(x = 0:26, y =plt_fit_ACA(Inter_ACA = Inter_ACA,
                                         slope_ACA = slope_ACA, 
                                           x = -0:26) )) + 
  ylim(0,10)


## Apply the function to correct
ACA_corr <- function(Area.29,d.15N.14N ){
  lind15N = d.15N.14N -log(Area.29)*slope_ACA
  return(lind15N)
}

N2_df$lind15N <- ACA_corr(Area.29 = N2_df$Area.29, d.15N.14N = N2_df$d.15N.14N)
  


################
## These are the values we know! SO we can use these to do the next correction
df_stds_N <- data.frame(ID = c("N1-1A", "N1.2A", "N2-1A", "N2-2A", "USGS.65.1", "USGS.65.2",  "IVA.Urea.1",  "IVA.Urea.2", "IVA.Urea.4"), 
                       Knownd15n =c(   .4,     .4,      20.3,     20.3,      20.68,        20.68,    -2.35,           -2.35,       -2.35) )

# c("N1-1A", "N1.2A", "N2-1A", "N2-2A", "USGS.65.1", "USGS.65.2"  "IVA.Urea.1",  "IVA.Urea.2", "IVA.Urea.4")
# c(   .4,     .4,      20.3,     20.3,      20.68,        20.68,    -2.35,           -2.35,       -2.35)


N2_df <- merge(df_stds_N, N2_df, by = "ID", all = TRUE)

## Fit a line
lm_fit_STDs <- summary(lm(lind15N~Knownd15n, data = N2_df ))
Inter_STDs <- lm_fit_STDs$coefficients[1,1]
slope_STDs <- lm_fit_STDs$coefficients[2,1]

plt_fit_STDs <- function(Inter_STDs, slope_STDs, x){
  y = slope_STDs*x + Inter_STDs
  return(y)
}

ggplot() +
  geom_point(data = N2_df, aes(x = Knownd15n, y = lind15N)) +
  geom_line(aes(x = -5:21, y =plt_fit_STDs(Inter_STDs = Inter_STDs,
                                          slope_STDs = slope_STDs, 
                                          x = -5:21) ))

## Apply the correction
STD_corr <- function(Inter_STDs,lind15N, slope_STDs ){
  find15N = (lind15N - Inter_STDs)/slope_STDs
  return(find15N)
}

N2_df$find15N <- STD_corr(Inter_STDs = Inter_STDs,
                          slope_STDs = slope_STDs, 
                          lind15N = N2_df$lind15N)

## Check out the samples!

samles <- paste0(c(rep("A", 12), rep("B", 12), rep("C", 6)), 
       c(1:12, 1:12, 1:6))
N2_df_smps <- N2_df[N2_df$ID %in% samles,]

hist(N2_df_smps$find15N)

write.csv(N2_df, file = "Batch_2_N.csv")
