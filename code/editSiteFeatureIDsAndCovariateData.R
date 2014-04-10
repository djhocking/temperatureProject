#Replacing covariate data for a site after visual inspection of site location resulting in changes to assigned NHD FeatureID.
# Read in the .csv for the site changes which indicates a 1 for OK, 0 for unknown issue, or the new FeatureID.
rm(list=ls())

# Define agency to edit
agency <- 'MEDMR'

setwd(paste0('C:/KPONEIL/GitHub/projects/temperatureProject/dataIn/', agency))

# Sites to be changed
s <- read.csv(paste0('siteChanges_', agency, '.csv'))
s$site <- as.character(s$site)
# Specific sites to change
s1 <- s[which(s$correctFeatureID > 1),]


# Dataframe to be edited
load(paste0('covariateData_', agency, '.RData'))
d <- covariateData
d$site <- as.character(d$site)

# Master dataframe
load('C:/KPONEIL/GitHub/projects/temperatureProject/dataIn/NENY_CovariateData_2014-03-13.RData')

# Don't want to replace Lat/Lon of the site with Lat/Lon of catchment
u <- UpstreamStats[, - which(names(UpstreamStats) %in% c('Latitude', 'Longitude'))]

# Remove old site covariate data
d1 <- d[!(d$site %in% s1$site), ]

# Set up dataframe for replacements
n <- s1[, c('site', 'correctFeatureID')]
colnames(n) <- c('site', 'FEATUREID')
n1 <- merge(n, d[,c('site', 'Latitude', 'Longitude')], by = 'site', all.x = T, sort = F)

# Merge in new covariate data
n2 <- merge(n1, u, by = 'FEATUREID', all.x = T, sort = F)

# Remove unused columns
n3 <- n2[,which(names(n2) %in% names(d1))]

# Join new to existing keeper data
covariateData <- rbind(d1, n3)


# Save the new dataframe
save(covariateData, file = paste0('covariateData_', agency, '.RData'))

