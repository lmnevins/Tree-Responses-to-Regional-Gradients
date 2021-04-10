######################################################################
#
#  PNW-TREE vital rates -- query USDA PLANTS traits
#
#    Rob Smith, phytomosaic@gmail.com, 10 Apr 2020
#
##      GNU General Public License, Version 3.0    ###################

### preamble
rm(list=ls())
require(rvest)
require(taxize)
require(ecole)

### load tree data
load('./data/pnw_tree.rda', verbose=T)
tree <- pnw_tree  ;  rm(pnw_tree)
names(tree)
(x <- sort(unique(tree$spp)))  # species names
x[x == 'populus_trichocarpa'] <- 'populus_balsamifera_trichocarpa'

### resolve names in USDA PLANTS
gr <- gnr_resolve(names = x, preferred_data_sources = c(150),
                  canonical=T, http='post', fields='all')
gr$traiturl <- paste0('https://plants.sc.egov.usda.gov/',
                      'java/charProfile?symbol=',
                      sub('.*?symbol=', '', gr$url))
(notfound <- x[which(x %notin% gr$user_supplied_name)])  # expect 0

### get search URLs from lookup table
x     <- data.frame(orignames=x, url=NA, stringsAsFactors=F)
x$url <- gr$traiturl[match(x$orignames, gr$user_supplied_name)]
rm(gr)
(notfound <- x$orignames[which(is.na(x$url))])  # expect 0

### manual changes (chamaecyparis_nootkatensis)
x$url <- gsub('symbol=CUNO', 'symbol=CANO9', x$url)
x$url <- gsub('symbol=CHCH7', 'symbol=CHCHC4', x$url)

### TIME WARN ! ! ! loop for all species in the dataset ! ! !
nspp  <- length(x$orignames)
lst   <- vector('list', nspp)
names(lst) <- x$orignames
for (j in seq_along(x$orignames)) {
        # j <- 1
        cat('iter', j, 'of', nspp, '--', x$orignames[j], '\n')
        u <- x$url[j]
        ### try read html
        tryCatch({
                a <- rvest::html_nodes(x=xml2::read_html(x=u),
                                       css='table')
        },
        error=function(cond) {
                message(paste0(u, ' does not exist'))
        })
        ### grab the html table
        tryCatch({
                a <- html_nodes(a, css='table')[7]
                a <- rvest::html_table(a, fill=F)[[1]]
        },
        error=function(cond) {
                message(paste0(u, ' does not exist'))
                a <- data.frame(matrix(NA, nr=88, nc=2))
        })
        ### assign
        lst[[j]] <- a
}
### post-processing
y <- lapply(lst, function(i) if(all(dim(i) == c(88,2))) i[,2] else NA)
y <- t(do.call(cbind, y))
dimnames(y)[[2]] <- ecole::clean_text(lst[[1]][,1], lower=T) ; rm(lst)
y <- y[,-NCOL(y)] # last column is blank
y[y==''] <- NA    # fill NA
y <- as.data.frame(y, stringsAsFactors=F)
j <- sapply(y, function(i) !length(unique(na.omit(i))) %in% c(0,1))
y <- y[,j] # rm invariant traits
dim(y)     # 57 species x 72 traits
j <- sapply(y, function(j) !all(is.na(as.numeric(j))))
y[,j] <- sapply(y[,j], as.numeric)

### FIA stem wood traits   ########################################
r <- read.csv(file='./data_raw/fia_data/REF_SPECIES.csv',
              stringsAsFactors=F)
names(r) <- ecole::clean_text(names(r), lower=T)
r$speciesname <- ecole::clean_text(paste0(r$genus, '_',
                                          r$species, '_',
                                          r$variety, '_',
                                          r$subspecies), lower=T)
r$speciesname[r$speciesname == 'abies_shastensis'] <-
        'abies_magnifica_var_shastensis'
r$speciesname[r$speciesname ==
                      'chrysolepis_chrysophylla_chrysophylla'] <-
        'chrysolepis_chrysophylla'
ptrn <- 'spgrpcd|exists|sitetree|sftwd|jenkins|raile|modifi|created|cit'
r <- r[, -grep(ptrn, colnames(r))]
r <- r[r$speciesname %in% dimnames(y)[[1]],-c(1:15)]
dimnames(r)[[1]] <- r$speciesname
j <- c('wood_spgr_greenvol_drywt', 'bark_spgr_greenvol_drywt',
       'mc_pct_green_bark', 'mc_pct_green_wood',
       'wood_spgr_mc12vol_drywt', 'bark_vol_pct')
y <- cbind(y, r[match(dimnames(y)[[1]], r$speciesname),j])
rm(r,j)

###################################################################
####    save    ##############################################
###################################################################
pnw_tra_usda <- y
save(pnw_tra_usda, file='./data/pnw_tra_usda.rda')
###################################################################
##############################################################
###################################################################



###   plot USDA traits data as heatmap   ##########################
require(ecole)
rm(list=ls())
load('./data/pnw_tra_usda.rda', verbose=T)
d <- pnw_tra_usda  ;  rm(pnw_tra_usda)
names(d)
dd <- data.frame(sapply(d, function(x) {
        if (is.character(x)) {
                as.numeric(as.factor(x))
        } else {
                x
        }
}))
dd <- data.frame(sapply(dd, ecole::standardize))
rownames(dd) <- rownames(d)
`plot_heatmap` <- function (x, labcex=0.7, xexp=1, yexp=1, ...) {
        col <- viridis::inferno(99, begin=0.15, end=0.95)
        op <- par(mar = c(0, 2.5 * xexp, 4 * yexp, 0) + 0.3)
        image(1:ncol(x), 1:nrow(x), t(x), col = col, axes = F, xlab = "",
              ylab = "", ...)
        axis(3, at = 1:ncol(x), labels = colnames(x), las = 3, tick = F,
             cex.axis = labcex)
        axis(2, at = 1:nrow(x), labels = rownames(x), las = 1, cex.axis = labcex)
        par(op)
}
png('./fig/fig_04_usda_traits_heatmap.png',6.5,6.5,'in',
    bg='transparent', res=1000)
ecole::set_par_mercury(1, oma = c(0, 0, 0, 0), font=1, tcl=-0.2,
                       mgp=c(1.0,0.4,0), pty='m')
plot_heatmap(dd, xexp=3.4, yexp=1.9, labcex=0.5)
dev.off()

####    end USDA PLANTS    ####