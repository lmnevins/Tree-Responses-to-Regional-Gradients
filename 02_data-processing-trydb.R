######################################################################
#
#  Tree Responses to Regional Gradients -- TRY-DB data processing
#
#    Rob Smith, phytomosaic@gmail.com, 29 Mar 2020
#
##      GNU General Public License, Version 3.0    ###################


### preamble
rm(list=ls())
# devtools::install_github('phytomosaic/ecole')
require(ecole)
require(data.table)
options(stringsAsFactors=F)
setwd('./data_raw/try_data/')

### function to read file
`f` <- function(x) {
        cat('File size (MB):',round(file.info(x)$size/1024^2,1),'--> ')
        cat('Time elapsed:',
            system.time({
                    out <- data.table::fread(x, header = T, sep = '\t',
                                             dec = '.', quote = '',
                                             fill = T, data.table = T)
            })[[3]], '\n')
        return(out)
}

# ### first, break up one ridiculously HUGE monster file (1.27 GB)
# cat('File size (MB):', round(file.info('./8540.txt')$size/1024^2),'\n')
# d <- f('./8540.txt') # ! ! ! TIMEWARN ! ! ! ~90 seconds
# nr  <- NROW(d) # n rows in this big file
# len <- 250000  # desired length of each smaller file
# nr / len       # how many smaller files to be produced
# i <- rep(1:ceiling(nr/len), ea=len)[1:nr] # create index
# # d[, indx := i]   # DONT add index to data.table!
# d[, fwrite(.SD, paste0('8540_',i,'.txt'),sep='\t'), by=i] # save small
# rm(list=ls())

### read ALL files in the directory ! ! ! TIMEWARN ! ! !
fnm <- list.files('.', pattern='.txt')   # all files in dir
fnm <- fnm[fnm != '_readme.txt']         # avoid the readme file
fnm <- fnm[fnm != '8540.txt']            # avoid the monster file
d   <- rbindlist(lapply(fnm, f)) ; rm(f) # TIMEWARN !
setnames(d, names(d), tolower(names(d))) # cleanup column names
rm(fnm, pth)                             # cleanup environment
d[,c('firstname','lastname','dataset',   # drop unneeded columns
     'accspeciesid','obsdataid','origlname','valuekindname',
     'origuncertaintystr','uncertaintyname','replicates',
     'reluncertaintypercent','reference','comment','v28') := NULL]
d <- d[order(accspeciesname, traitid)]   # sort for convenience

### clean species names
d[, speciesname    := ecole::clean_text(speciesname,lower=T)]
d[, accspeciesname := ecole::clean_text(accspeciesname,lower=T)]

### remove one blank row
d <- d[!c(traitname == '' & dataname == ''),]

### remove duplicates based on 'origobsdataid'
i <- which(duplicated(d,by='origobsdataid') & !is.na(d$origobsdataid))
length(i) / NROW(d) # 2% = propn of obsvns that are dupes
d <- d[-i,]  ;  rm(i)

### keep only species of interest
load('../../data/pnw_tree.rda', verbose=T)
(x <- sort(unique(pnw_tree$spp)))  ;  rm(pnw_tree)    # FIA species
y <- sort(unique(c(d$speciesname, d$accspeciesname))) # TRY species

### strict EXACT match of spp -- beware synonyms/orthovariants...
all(x %in% y)          # all 57 FIA spp are present in TRY
(keep <- y[y %in% x])  # TRY in FIA -- all 57 FIA species!
d <- d[accspeciesname %in% keep] # strict EXACT match

### summary
cat('Object size (MB):', round(object.size(d)/1024^2,2))
names(d)
hist(table(d$observationid), breaks=33, col='grey',
     xlab='Measurements per entity', main='')
table(d$dataid == 59)         # how many latitude observations?
sort(table(d$speciesname))    # old/mixed species names
sort(table(d$accspeciesname)) # accepted species names
hist(table(d$accspeciesname), breaks=33, col='grey',
     xlab='Measurements per species', main='')
table(d$datasetid)
table(d$unitname)
table(d$traitname)
sort(table(d$traitname))
sort(table(d$dataname))

### save
pnw_tra_trydb <- d
save(pnw_tra_trydb, file='../../data/pnw_tra_trydb.rda')

####    END    ####