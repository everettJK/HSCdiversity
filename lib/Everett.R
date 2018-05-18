library(stringr)
library(reshape2)
library(pander)
library(parallel)
library(reldist)
library(dplyr)
library(vegan)
source('abundanceCorrections.R')
options(stringsAsFactors = FALSE)

#
# Read in intSite data.
# The intSite positions and abundances are the result of processing the output of the INSPIIRED 
# intSite pipeline with C Berry's dRep software package.
#

intSites <- read.csv(gzfile('../data/intSites.csv.gz'), header = TRUE)
for(i in c('start', 'end', 'width', 'estAbund')) intSites[[i]] <- as.numeric(intSites[[i]]) 
intSites[which(is.na(intSites$nearestFeature)),]$nearestFeature <- 'NONE'



# Read in the list of samples that are to be exported to the SRA.
SRA_samples <- read.table('../data/SRA_samples.txt', sep = ';', header = TRUE)

# For plot clarity, change all m12_6 time points to m12.
# This only pertains to paient FR04.

intSites$timepoint <- as.character(intSites$timepoint)
intSites[which(intSites$timepoint == 'm12_6'),]$timepoint <- 'm12'


#
# Restrict data set to SRA samples.
#

all(SRA_samples$SpecimenAccNum %in% intSites$sampleName)
intSites <- subset(intSites, sampleName %in% SRA_samples$SpecimenAccNum)



#
# Rename cell types.
#

intSites$celltype <- as.character(intSites$celltype)
intSites$celltype <- gsub('NEUTROPHILS', 'GRANULOCYTES', intSites$celltype)  


#
# Convert the time point to days.

# Convert the timepoint character string to a numeric days column.
intSites$timePointDays <- toupper(intSites$timepoint)
intSites$timePointDays <- gsub('_', '.', intSites$timePointDays)
intSites$timePointDays <- sapply(intSites$timePointDays, function(x){
  if(grepl('D', x)){
    return(as.numeric(str_match(x, '[\\d\\.]+')))
  } else if(grepl('M', x)){
    return(as.numeric(str_match(x, '[\\d\\.]+'))*30.4167)
  } else if(grepl('Y', x)){
    return(as.numeric(str_match(x, '[\\d\\.]+'))*365)
  } else {
    message('Warning - could not determine date unit for: ', x)
    return(as.numeric(str_match(x, '[\\d\\.]+')))
  }
})



#
# Population estimates 

ise <- read.csv(gzfile('../data/berryDrep.csv.gz'))
ise$key <- toupper(paste(ise$patient, ise$sampleName, ise$seqnames, ise$strand, ise$start))
intSites$key <- toupper(paste(intSites$patient, intSites$sampleName, intSites$seqnames, intSites$strand, intSites$start))
table(intSites$key %in% ise$key)

ise.readsCounts <- 
  dplyr::select(ise, key, readCounts) %>%
  dplyr::group_by(key) %>%
  dplyr::mutate(totalReads = sum(readCounts)) %>%
  dplyr::ungroup()

intSites$reads <- ise.readsCounts[match(intSites$key, ise.readsCounts$key),]$totalReads

populationEstimates <-
  bind_rows(lapply(split(intSites, intSites$patient), function(x){
    x$TP <- ifelse(x$timePointDays <= 690, 'TP1','TP2')
    bind_rows(lapply(split(x, x$TP), function(x2){
      message(x2$patient[1], ':', x2$TP[1], ':', x2$timepoint[1])
      group_by(x2, posid) %>%
        mutate(nSamples = n_distinct(sampleName)) %>%
        top_n(1,  estAbund) %>%
        dplyr::slice(1) %>%  
      ungroup() %>%
      summarise(patient = x2$patient[1],
                TP      = x2$TP[1],
                percentSharedSites = round((sum(nSamples > 1)/n())*100, 2),
                uniqueSites        = n_distinct(x2$posid),
                chao1_frags = round(estimateR(x2$estAbund, index='chao')[2], 0),
                chao1_reads = round(estimateR(x2$reads, index='chao')[2], 0))       
    }))
})) %>% data.frame()




populationSizeEstimates <- 
  intSites %>% 
  mutate(TP = ifelse(timePointDays <= 690, 'TP1','TP2')) %>%
  mutate(patientPosid = paste(patient, posid)) %>%
  select(TP, patientPosid, estAbund) %>%
  group_by(TP, patientPosid) %>%
  top_n(1,  estAbund) %>%
  dplyr::slice(1) %>%  
  ungroup() %>%
  summarise(TP1_unique = n_distinct(patientPosid[TP=='TP1']),
            TP1_chao1  =  round(estimateR(estAbund[TP=='TP1'], index='chao')[2], 0),
            TP2_unique = n_distinct(patientPosid[TP=='TP2']),
            TP2_chao1  =  round(estimateR(estAbund[TP=='TP2'], index='chao')[2], 0))


#
# Create a second labeled nerest gene column than contains gene type flags.
#

intSites$nearestFeature2 <- 
  intSites %>%
  mutate(nearestFeature2 = paste0(nearestFeature, ' ')) %>% 
  mutate(nearestFeature2 = ifelse(!is.na(inFeature), paste0(nearestFeature2, '*'), nearestFeature2)) %>%
  mutate(nearestFeature2 = ifelse(abs(nearestOncoFeatureDist) <= 50000, paste0(nearestFeature2, '~'), nearestFeature2)) %>%
  mutate(nearestFeature2 = ifelse(abs(nearestlymphomaFeatureDist) <= 50000, paste0(nearestFeature2, '!'), nearestFeature2)) %>%
  select(nearestFeature2) %>%
  unlist() %>%
  unname()




#
# Read in the cell sorting cross over reports which are stored in a single file with records separted with '#%'.
# The counts table is identifiable by the key word COUNTS. 
#


crossOverReports <- readChar('../data/crossOverReports.tsv', file.info('../data/crossOverReports.tsv')$size)
crossOverReports <- unlist(strsplit(crossOverReports, '#%'))
crossOverReports <- unlist(lapply(crossOverReports, function(x){
  source         <- ifelse(grepl('#source:', x), sub('\\t+$', '', str_match(x, '#source:\\s+(.+)')[2]), '')
  patient        <- str_match(x, '#(\\S+)')[2]
  timePoint      <- str_match(x, '#\\S+\\s+(\\S+)')[2]
  t              <- gsub('\\t+\\n', '\n', substr(x, regexpr('COUNTS', x)+7, nchar(x)))
  m              <- read.table(tc <- textConnection(t), header = TRUE, fill = TRUE, check.names = FALSE, sep='\t'); close(tc); 
  
  # If the matrix is less than 5 x 5 then return NA
  if(nrow(m) < 5 | length(m) < 5) return(NA)
  
  # Some matrices will contain additional columns and rows - exclude these extra dimensions.
  m <- m[1:5, 1:5]
  
  # Return the data entry as a named list. 
  r <- list()
  key <- paste(patient, timePoint, source, sep='|')
  
  r[[key]] <- list(
     patient       = patient,
     timePoint     = timePoint,
     source        = source,
     table         = m,
     cellCount     = as.numeric(str_match_all(str_match(x, 'initialCellCounts,(.+)')[2], '([\\d\\.]+)')[[1]][,2]),
     pB_postsort   = as.numeric(str_match_all(str_match(x, 'pB_postsort,(.+)')[2], '([\\d\\.]+)')[[1]][,2]),
     postsortCount = as.numeric(str_match_all(str_match(x, 'postsortCount,(.+)')[2], '([\\d\\.]+)')[[1]][,2]),
     VCN           = as.numeric(str_match_all(str_match(x, 'VCN,(.+)')[2], '([\\d\\.]+)')[[1]][,2]))
  
  r
}), recursive = FALSE)

crossOverReports <- crossOverReports[sapply(crossOverReports, is.list)]



#
# Subset the resulting crossover reports to include only those reports which 
# will be used in the study and throw an error if any expected reports are missing.
#

samples <- c('FR01|m22|Blood',  'FR01|m48|Blood',       
             'FR03|m12|Blood',  'FR03|m36|Blood',       
             'FR04|m12|Blood',  'FR04|m36|Blood',      
             'p1202|m36|Blood', 'p1202|m11|Blood', 
             'p1204|m12|Blood', 'p1204|m24|Blood',
             'FR05|m30|Blood')

if(! all(samples %in% names(crossOverReports))) stop('All the requested crossover tables were not found.')
crossOverReports <- crossOverReports[samples]  



#
# Utilize the data in each crossover report to correct the intSite abundances.
# Abundance columns estAbund2 and estAbund3 will be used to store different abundance corrections.
# The two different corrections methods will store their results in intSites$estAbund2 and intSites$estAbund3.
#

intSites$estAbund2 <- NA
intSites$estAbund3 <- NA


#
# Create a key to be used to update corrected abundances.
#

intSites$key <- toupper(paste(intSites$posid, intSites$celltype, intSites$patient, intSites$timepoint))


invisible(lapply(crossOverReports, function(x){
  if(x$source == 'BoneMarrow') names(x$table) <- paste0('BM_', names(x$table))
  
  # Subset the intSite data to include only sites from a specific cell type and time point.
  i <- which(intSites$patient == x$patient & intSites$timepoint == x$timePoint & toupper(intSites$celltype) %in% toupper(names(x$table)))
  d <- intSites[i,]
  
  if(nrow(d) == 0) stop(paste0(x$patient, ' / ', x$timePoint, ' could not be found in the intSite data.'))
  
  message('Cell types in retrieved data subset(', x$patient, ' - ', x$timePoint, '): ', paste0(unique(d$celltype), collapse=', '))
  
  # Reorganize the data to create an intSite / cell count table.
  d2 <- reshape2::dcast(d, posid ~ celltype, value.var='estAbund', fun.aggregate=function(x){x[1]}, fill=0)
  
  # Add missing cell types.
  d2[names(x$table)[! toupper(names(x$table)) %in% toupper(names(d2))]] <- 0
  
  # Add nearest gene column.
  d2$gene <- intSites[match(d2$posid, intSites$posid),]$nearestFeature
  
  ###browser()
  
  cols <- c(1, grep('gene',   names(d2), ignore.case = TRUE),
            grep('Bcell',  names(d2), ignore.case = TRUE),
            grep('Mono',   names(d2), ignore.case = TRUE),
            grep('Granu', names(d2), ignore.case = TRUE),
            grep('NKcell', names(d2), ignore.case = TRUE),
            grep('Tcell',  names(d2), ignore.case = TRUE))
  

  if(any(duplicated(cols))) stop('Duplicate column headers found.')
  if(length(cols) != 7) stop('All the expected columns were not found.')
  
  # Reorganize the cloumn headers to match Correction_CutData_new() input structure.
  d2 <- d2[, cols]
  
  # Create a cell count column.
  d2$cellCount <-  apply(d2, 1, function(x){ sum(as.integer(x[3:7])) })
  
  # Rename the input columns for Correction_CutData_new().
  names(d2)=c("integrationSite", "gene", "B", "M", "N", "K", "T", "cellCount") 
  
  #
  # The intSite data has been reorganized into a count like table and formatted for Correction_CutData_new().
  # Now we replace all zeros in the crossover table with 0.1 and call Correction_CutData_new().
  #
  
  x$table[x$table == 0] <- 0.1
  
  message('Correction input header:')
  pander(head(d2, n=3))
  
  correction <- Correction_CutData_new(d2, 2, as.matrix(x$table), x$cellCount, x$pB_postsort)
  
  # Rename the returned data frame to match the cell types in the intSites data frame.
  n <- c("integrationSite", "gene", "Bcells", "Monocytes", "Granulocytes", "NKcells", "Tcells")
  if(x$source == 'BoneMarrow') n[3:7] <- paste0('BM_', n[3:7]) 
  names(correction$Cut_data_corrected2C) <- n
  names(correction$Cut_data_corrected3)  <- n
  
  # Update the intSites data frame with the corrected1 values (stored in estAbund2)
  m <- melt(correction$Cut_data_corrected2C)
  m$key <- toupper(paste(m$integrationSite, m$variable, x$patient, x$timePoint))
  m$i <- match(m$key, intSites$key)
  m <- m[!is.na(m$i),]
  intSites[m$i,]$estAbund2 <<- m$value
  
  # Update the intSites data frame with the corrected2C values (stored in estAbund3)
  m <- melt(correction$Cut_data_corrected3)
  m$key <- toupper(paste(m$integrationSite, m$variable, x$patient, x$timePoint))
  m$i <- match(m$key, intSites$key)
  m <- m[!is.na(m$i),]
  intSites[m$i,]$estAbund3 <<- m$value
}))




#
# Reduce the intSite data to only the 5 cell types of interest (both BM and non-BM)
# and the time points used in the study.
#

cellTypes <- toupper(c('Bcells', 'Monocytes', 'Granulocytes', 'NKcells', 'Tcells'))
intSites <- subset(intSites, celltype %in% cellTypes)

# Save a copy on intSites before trimming it down to only the patient time points for the manuscript.
intSites.allTimePoints <- intSites

intSites <- rbind(subset(intSites, patient=='FR01'  & timepoint %in% c('m22',   'm48')),
                  subset(intSites, patient=='FR03'  & timepoint %in% c('m12',   'm36')),
                  subset(intSites, patient=='FR04'  & timepoint %in% c('m12', 'm36')),
                  subset(intSites, patient=='FR05'  & timepoint %in% c('m12',   'm30')),
                  subset(intSites, patient=='p1202' & timepoint %in% c('m11',   'm36')),
                  subset(intSites, patient=='p1204' & timepoint %in% c('m12',   'm24')))

intSites$celltype  <- factor(as.character(intSites$celltype),  levels=c(unique(as.character(intSites$celltype))))
intSites$timepoint <- factor(as.character(intSites$timepoint), levels=c(unique(as.character(intSites$timepoint))))


# Create a column denoting source (required for plots)
intSites$sampleType <- 'MPB'
intSites[which(intSites$patient=='FR03' | intSites$patient=='FR04' | intSites$patient=='p1204' ),]$sampleType <- 'BM'

# Clean up unneeded columns.
intSites$key <- NULL

# Create wide view of the data.
intSites.wide <- reshape2::dcast(intSites, posid ~ celltype, value.var='estAbund', fun.aggregate=function(x){x[1]}, fill=0)
intSites.wide$gene <- intSites[match(intSites.wide$posid, intSites$posid),]$nearestFeature
intSites.wide <- intSites.wide[, c(1, grep('gene',   names(intSites.wide), ignore.case = TRUE),
                                   grep('Bcell',  names(intSites.wide), ignore.case = TRUE),
                                   grep('Mono',   names(intSites.wide), ignore.case = TRUE),
                                   grep('Gran', names(intSites.wide), ignore.case = TRUE),
                                   grep('NKcell', names(intSites.wide), ignore.case = TRUE),
                                   grep('Tcell',  names(intSites.wide), ignore.case = TRUE))]

intSites.wide$cellCount <-  apply(intSites.wide, 1, function(x){ sum(as.integer(x[3:7])) })
names(intSites.wide)=c("integrationSite", "gene", "B", "M", "N", "K", "T", "cellCount") 










# Select integration sites that are within 50KB of any gene and then determine
# which intSites where seen in the early timepoint and which were seen in the later timepoint.
# 690 days (23 months) is used as a dividing timepoint between all timepoints 1 and all timepoints 2
# for all subjects.

CPUs <- 25
dividingTimePointDays <- 690

cluster <- makeCluster(CPUs)
d <- intSites[which(abs(intSites$nearestFeatureDist) <= 50000),]
oncoGeneList <- scan(file='../data/humanOncoGenes.txt', what='character')
clusterExport(cl=cluster, list('oncoGeneList','dividingTimePointDays'), envir=.GlobalEnv)  

pvp.plotData <- do.call(rbind, parLapply(cluster, split(d, d$nearestFeature), function(x){
  pre  <- subset(x, timePointDays <= dividingTimePointDays)
  post <- subset(x, timePointDays >  dividingTimePointDays)
  
  oncoGene <- toupper(x$nearestFeature[1]) %in% toupper(oncoGeneList)
  
  data.frame(gene       = x$nearestFeature, 
             preIntNum  = ifelse(length(unique(pre$posid))  > 0, length(unique(pre$posid)),  0.001), 
             postIntNum = ifelse(length(unique(post$posid)) > 0, length(unique(post$posid)), 0.001), 
             oncoGene   = oncoGene, 
             maxAbund   = max(x$estAbund))
}))
stopCluster(cluster)

# Determine (sites near x) / (all sites) ratio for each gene.
pvp.plotData$preIntFreq  <- pvp.plotData$preIntNum/length(unique(subset(intSites,  timePointDays <= dividingTimePointDays)$posid))
pvp.plotData$postIntFreq <- pvp.plotData$postIntNum/length(unique(subset(intSites, timePointDays > dividingTimePointDays)$posid))

pvp.plotData <- pvp.plotData[!duplicated(pvp.plotData),]
row.names(pvp.plotData) <- NULL




#
# Create a simliar data frame for d0 vs last time point.
#

d0 <- data.frame(readRDS('../data/WAS_d0_intSites.rds'))
d0 <- d0[which(abs(d0$nearestFeatureDist) <= 50000),]
d0$timePointDays <- 0

d <- intSites[which(abs(intSites$nearestFeatureDist) <= 50000),]

d  <- rbind(d[,  c('nearestFeature', 'timePointDays', 'posid', 'estAbund')], 
            d0[, c('nearestFeature', 'timePointDays', 'posid', 'estAbund')])

d0vs2ndTimePoint <- subset(d, timePointDays < 1 | timePointDays > dividingTimePointDays)

cluster <- makeCluster(CPUs)
clusterExport(cl=cluster, list('oncoGeneList', 'dividingTimePointDays'), envir=.GlobalEnv)  

pvp.plotData2 <- do.call(rbind, parLapply(cluster, split(d0vs2ndTimePoint, d0vs2ndTimePoint$nearestFeature), function(x){
  pre  <- subset(x, timePointDays == 0)
  post <- subset(x, timePointDays >  0)
  
  oncoGene <- toupper(x$nearestFeature[1]) %in% toupper(oncoGeneList)
  
  data.frame(gene       = x$nearestFeature, 
             preIntNum  = ifelse(length(unique(pre$posid))  > 0, length(unique(pre$posid)),  0.001), 
             postIntNum = ifelse(length(unique(post$posid)) > 0, length(unique(post$posid)), 0.001), 
             oncoGene   = oncoGene, 
             maxAbund   = max(x$estAbund))
}))
stopCluster(cluster)

# Determine (sites near x) / (all sites) ratio for each gene.
pvp.plotData2$preIntFreq  <- pvp.plotData2$preIntNum/length(unique(subset(d0vs2ndTimePoint,  timePointDays == 0)$posid))
pvp.plotData2$postIntFreq <- pvp.plotData2$postIntNum/length(unique(subset(d0vs2ndTimePoint, timePointDays > 0)$posid))

pvp.plotData2 <- pvp.plotData2[!duplicated(pvp.plotData2),]
row.names(pvp.plotData2) <- NULL


#
# Cell name shorthands for figures.
#

### cellNameAbbrv <- list('TCELLS'='T', 'MONOCYTES'='M', 'GRANULOCYTES'='N', 'BCELLS'='B', 'NKCELLS'='K')



#
# Save all data objects for knittr.
#

save.image(file='Everett.RData')




#
# Create genomic heat maps
#

Rscript_path <- '/home/opt/R-3.4.0/bin/Rscript'

# Remove intSites in non-standard chromosomes since we already have to remove some sites
# in order not exceede the capacity of the calcualtion.

d <- intSites
d <- subset(d, seqnames %in% c(paste0('chr', 1:22), 'chrX', 'chrY'))

# Create a timepoint / cell type splitting factor.
d$s <- paste(d$timepoint, d$celltype, sep='.')

# Create the sample files for the genomic heatmap software.
# Rather than replicate samples, use the splitting factor as the sample otherwise
# the returned heatmap would have an unreasonable number of columns.

if(! dir.exists('../data/tmp')) dir.create('../data/tmp', showWarnings = FALSE)

invisible(lapply(split(d, d$patient), function(x){
  df <- do.call(rbind, lapply(split(x, x$s), function(x2){
    data.frame(sampleName=x2$s[1],
               GTSP=x2$s[1],
               patient=x2$patient[1],
               cellType=x2$celltype[1])
  }))
  
  df <- df[order(df$cellType),]
  df$cellType <- NULL
  fileName <- paste0('../data/tmp/', x$patient[1], '.samples.csv')
  write.table(df, quote = FALSE, row.names = FALSE, col.names = TRUE, file = fileName, sep=',')
}))


# For each of the generated sample files, create a script that will create the corresponding intSite data files.
# The genomic heatmap maker appears to be limted to 20,000 inSites.

#
# Create command script to create genomic heatmaps.
#

write('# Command script', file='../data/tmp/commandScript.sh', append = FALSE)

invisible(lapply(list.files(path='../data/tmp', pattern='\\.samples.csv'), function(file){
  sampleTable <- read.csv(file.path('../data/tmp', file))
  o <- subset(d, patient==sampleTable$patient[1] & s %in% sampleTable$GTSP)
  df <- data.frame(o$seqnames, o$strand, o$start, o$s, 'hg38') 
  names(df) <- c('seqnames', 'strand', 'position', 'sampleName', 'refGenome')
  
  if(nrow(df) > 20000) df <- df[sample(1:nrow(df), 20000),]
  
  sitesFileName <- file.path('../data/tmp', sub('samples', 'sites', file))
  write.table(df, quote = FALSE, row.names = FALSE, col.names = TRUE, file = sitesFileName, sep=',')
  dir.create(file.path('../data/tmp', sampleTable$patient[1]), showWarnings = FALSE)
  comm <- paste(Rscript_path, './software/genomicHeatmapMaker-from_input/genomic_heatmap_from_file.R', 
                file.path('../data/tmp', file),
                '-c ./software/genomicHeatmapMaker-from_input/INSPIIRED.yml',
                '-o', file.path('../data/tmp', sampleTable$patient[1]),
                '-f', sitesFileName,
                '-r hg38')
  write(paste('echo', sampleTable$patient[1]), file='../data/tmp/commandScript.sh', append = TRUE)
  write(comm, file='../data/tmp/commandScript.sh', append = TRUE)
}))


#
# Create command script to create epigenetic heatmaps.
#

write('# Command script', file='../data/tmp/epiCommandScript.sh', append = FALSE)

invisible(lapply(list.files(path='../data/tmp', pattern='\\.samples.csv'), function(file){
  sampleTable <- read.csv(file.path('../data/tmp', file))
  o <- subset(d, patient==sampleTable$patient[1] & s %in% sampleTable$GTSP)
  df <- data.frame(o$seqnames, o$strand, o$start, o$s, 'hg38') 
  names(df) <- c('seqnames', 'strand', 'position', 'sampleName', 'refGenome')
  
  sitesFileName <- file.path('../data/tmp', sub('samples', 'EpiSites', file))
  write.table(df, quote = FALSE, row.names = FALSE, col.names = TRUE, file = sitesFileName, sep=',')
  dir.create(file.path('../data/tmp', paste0('EPI_', sampleTable$patient[1])), showWarnings = FALSE)
  
  comm <- paste('/home/opt/R-3.4.0/bin/Rscript ./software/EpigeneticHeatmapMaker-from_input/epi_heatmap_from_file.R', 
                file.path('../data/tmp', file),
                '-c  ./software/EpigeneticHeatmapMaker-from_input/INSPIIRED.yml',
                '-o', file.path('../data/tmp', paste0('EPI_', sampleTable$patient[1])),
                '-t epiCellTypes',
                '-f', sitesFileName,
                '-r hg38')
  write(paste('echo', sampleTable$patient[1]), file='../data/tmp/epiCommandScript.sh', append = TRUE)
  write(comm, file='../data/tmp/epiCommandScript.sh', append = TRUE)
}))


#
# Create heat map files required for 'All subjects' heat maps.
#

if(! file.exists('../data/tmp/all.sites')) system(paste('cat', paste(list.files('../data/tmp', pattern='sites', full.names = TRUE), collapse = ' '), '> ../data/tmp/all.sites'))
if(! file.exists('../data/tmp/all.samples')) system(paste('cat', paste(list.files('../data/tmp', pattern='samples', full.names = TRUE), collapse = ' '), '> ../data/tmp/all.samples'))

sites <- read.csv('../data/tmp/all.sites', header = FALSE)
sites <- sites[-which(sites$V1=='seqnames'),]
names(sites) <- c('seqnames','strand','position','sampleName','refGenome')
sites$sampleName <- gsub('m\\d+_?\\d+\\.', '', sites$sampleName)
sites <- sites[!duplicated(sites),]
write.csv(sites, file='../data/tmp/all.sites', quote = FALSE, row.names = FALSE)

sites <- sites[sample(1:nrow(sites), 20000),]
write.csv(sites, file='../data/tmp/all.20000.sites', quote = FALSE, row.names = FALSE)

samples <- read.csv('../data/tmp/all.samples', header = FALSE)
samples <- samples[-which(samples$V1=='sampleName'),]
names(samples) <- c('sampleName','GTSP','patient')

samples$sampleName <- gsub('m\\d+_?\\d+\\.', '', samples$sampleName)
samples$GTSP <- gsub('m\\d+_?\\d+\\.', '', samples$GTSP)
samples$patient <- 'xxx'
samples <- samples[!duplicated(samples),]
write.csv(samples, file='../data/tmp/all.samples', quote = FALSE, row.names = FALSE)




comm <- paste(Rscript_path, '../software/genomicHeatmapMaker-from_input/genomic_heatmap_from_file.R', 
              '../data/tmp/all.samples',
              '-c ../software/genomicHeatmapMaker-from_input/INSPIIRED.yml',
              '-o', file.path('../data/tmp/GEN_all'),
              '-f ../data/tmp/all.20000.sites',
              '-r hg38')
write('echo all subjects', file='../data/tmp/commandScript.sh', append = TRUE)
write(comm, file='../data/tmp/commandScript.sh', append = TRUE)



comm <- paste('/home/opt/R-3.4.0/bin/Rscript ../software/EpigeneticHeatmapMaker-from_input/epi_heatmap_from_file.R', 
              '../data/tmp/all.samples',
              '-c  ../software/EpigeneticHeatmapMaker-from_input/INSPIIRED.yml',
              '-o ../data/tmp/EPI_all',
              '-t ../data/epiCellTypes',
              '-f ../data/tmp/all.20000.sites',
              '-r hg38')
write('echo all subjects', file='../data/tmp/epiCommandScript.sh', append = TRUE)
write(comm, file='../data/tmp/epiCommandScript.sh', append = TRUE)

dir.create('../data/tmp/EPI_all')
dir.create('../data/tmp/GEN_all')


#-----------------------------------------------v--------------------------------------------------#
#                 Manualy run the command scripts to create the genomic heat maps.                 #
#-----------------------------------------------^--------------------------------------------------#


