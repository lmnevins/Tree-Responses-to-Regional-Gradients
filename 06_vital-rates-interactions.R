######################################################################
#
#  PNW-TREE vital rates -- linear model INTERACTIONS
#
#     Rob Smith, phytomosaic@gmail.com, 05 Oct 2020
#
##      GNU General Public License, Version 3.0    ###################

### preamble
rm(list=ls())
# devtools::install_github('phytomosaic/ecole')
require(ecole)
require(data.table)
require(ggplot2)

### plotting specs
# x1 <- expression(log[10] ~ CMD ~ (mm ~ y^{-1}))
# x2 <- expression(Fire ~ probability^{1/3})
# x3 <- 'Composite stress'
# # y1 <- 'Abundance (QMD)'
# y2 <- 'Growth (RGR)'
# y3 <- 'Survival'
# y4 <- 'Crowding (NCI)'
# y5 <- 'Func neigh (PC1)' # NCIS
# y6 <- 'Func neigh (PC2)' # NCIS

### load data
setwd("~/Documents/WSU/R Scripts")
load('pnw_tree-3.rda', verbose=T)
d <- pnw_tree  ;  rm(pnw_tree)
### read functional neighborhood SES per individual, from Jenny
#     NCIS = func neighborhood values =
#        pairwise trait deviations, multiplied by DBH^2, divided by dist^2
#     SES = deviation from random expectation, scaled by SD of random draws
#       problem: should scale and center traits first
#       problem: SES involving congeners always = 0 since traits identical
d[, c('pc1','pc2') := NULL] # rm climate PCs (avoid name conflict w trait PCs)
setwd("~/Dropbox/Trait_data/Results/Functional_neighborhood_SES")
fnm <- c('ses_PC1_alltraits.csv','ses_PC2_alltraits.csv')
ses <- lapply(fnm, fread, colClasses = c('character', 'numeric'))
ses <- Reduce(merge,ses)
names(ses) <- gsub('stand_size_', '', tolower(names(ses)))
names(ses)[names(ses)=='focal'] <- 'tre_cn'
d   <- merge(d, ses, by='tre_cn')
rm(ses)

### exclude juveniles < 5.0 inches (12.7 cm) per Fenn et al. (2020)
i <- !((!is.na(d$dia_y1) & d$dia_y1 < 5) |
           (!is.na(d$dia_y2) & d$dia_y2 < 5))
d <- d[i,]

### transformations
d$pai     <- d$dia_y2 - d$dia_y1  # annual increment (untransformed)
d$dia_y1  <- log10(1 + d$dia_y1)  # log10-transform
d$dia_y2  <- log10(1 + d$dia_y2)  # log10-transform
d$annincr <- d$dia_y2 - d$dia_y1  # annual increment (transformed)
d$rgr     <- d$annincr / d$dia_y1 # relative growth rate (relative to y1)

# change name format to remove underscore and capitalize genus 
d$spp <- gsub("_", " ", d$spp)
d$spp <- gsub("^(\\w)(\\w+)", "\\U\\1\\L\\2",
              d$spp, perl = TRUE)

# ### quadratic mean function
# `quad_mean` <- function(x, na.rm=TRUE, zero.rm=FALSE, ...){
#   if(any(x < 0, na.rm=TRUE)) stop('quadratic mean not defined for neg values')
#   if(zero.rm) x[x==0] <- NA
#   return(sqrt(mean((x^2), na.rm=na.rm)))
# }
### summarize means of SPP in PLOT (SPP/PLOT = sample unit)
p <- d[, list(
    # startdiam    = mean(dia_y1, na.rm=T),      # mean size
    # qmd          = quad_mean(dia_y1, na.rm=T), # quad mean size
    surv         = mean(surv, na.rm=T),
    rgr          = mean(rgr, na.rm=T),
    plt_nci      = median(i_nci, na.rm=T),
    plt_ses_pc1  = mean(pc1, na.rm=T),
    plt_ses_pc2  = mean(pc2, na.rm=T)
), by=list(plt_cn, spp)] ### <-- (SPP/PLOT = sample unit)
# remove all the INDIVIDUAL-level measures
j <- !names(d) %in% c('tre_cn','dist','subp','subpid','x','y',
                      'is_seed','is_edge','dia_y1','dia_y2','annincr',
                      'plt_cn','spp','qmd','surv','rgr','plt_nci',
                      'i_nci','i_ncilog','pc1','pc2','pai')
p <- cbind(p, d[match(p$plt_cn, d$plt_cn), ..j]) ; rm(j) # plot match

### trim extreme NCI values > 99.9 percentile
p$plt_nci[p$plt_nci > 2] <- NA
p$rgr[p$rgr > 0.10]      <- NA
is.na(p) <- is.na(p)  # force all NaN to NA
# rm(d)

### frequencies
frq <- rev(sort(table(d$spp)))
nm  <- names(frq)

### calculate 'composite stress' gradient (sum of fire/drought percentiles)
`ptile` <- function(x) ecdf(x)(x)
a  <- cbind(fire=ptile(p$firelog), cmd=ptile(p$cmdlog)) # percentiles
p$stress <- rowMeans(a, na.rm=T) ; rm(a, ptile) # max value = 1 is intuitive

### END data prep ###########################################################



###   test interactions    #################################################

# does the effect of functional neighborhoods depend on position
#     along gradients (CMD or FIRE)?
p$INT_CMD  <- p$cmdlog  * p$plt_ses_pc1  # interaction term (CMD)
p$INT_FIRE <- p$firelog * p$plt_ses_pc1  # interaction term (fire)
set.seed(88)
n_boot <- 99 # ! ! ! TIMEWARN ! ! !
n_top  <- 20 # take only top-20 most frequent spp
j      <- c('rgr','surv')        # responses
j_pred <- c('INT_CMD','INT_FIRE') # predictors
# 300 rows / 20 species / 2 response vars / 6 predictors
`sens` <- function(X, y) {
    s <- sensitivity::src(X, y, nboot=n_boot)$SRC
    c(sen = s[['original']] - s[['bias']],
      lwr = s[['min. c.i.']], upr = s[['max. c.i.']])
}
s <- lapply(j, function(yvar) { # apply for each combn of species, predictor and response
    sapply(nm[1:n_top], function(i) {
        pp <- p[p$spp==i,]        # species subset
        jj <- c(j_pred,yvar)      # predictor/response subset
        pp <- na.omit(pp[, ..jj]) # NA omit
        sens(data.frame(pp[,..j_pred]), c(pp[,..yvar])) # sensitivity
    })
})
s <- do.call(rbind.data.frame, s)
s$response  <- rep(j, ea=length(j_pred)*3)    # 3 metrics per predictor
s$predictor <- rep(j_pred, times=length(j)*3) # 3 metrics per response
s$metric    <- rep(rep(c('sen','lwr','upr'), ea=length(j_pred)), times=length(j))
s <- reshape(s,
             idvar=c('response','predictor','metric'),
             varying=list(nm[1:n_top]),
             v.names='val',
             timevar='spp',
             times=nm[1:n_top],
             direction='long',
             new.row.names = 1:9999)
s <- reshape(s,
             idvar=c('response','predictor','spp'),
             timevar='metric',
             v.names='val',
             direction='wide')
names(s)    <- gsub('val.','',names(s))
s$response  <- plyr::mapvalues(s$response, j, c('Growth','Survival'))
s$predictor <- plyr::mapvalues(s$predictor, j_pred,
                               c('CMD x Func neigh','Fire x Func neigh'))
s$predictor <- factor(s$predictor)
s$posneg    <- factor(ifelse(s$lwr > 0.0, 1, ifelse(s$upr < (-0.0), -1, 0)))
s$bootrng   <- s$upr - s$lwr  # <<--- sensitivity
### sort species by STRESS
x <- subset(p, spp %in% nm[1:20]) # subset top-20 species
x <- x[, lapply(.SD, mean, na.rm=T), by=spp, .SDcols='stress']
setorder(x, stress, spp) # order by STRESS
s$spp    <- factor(s$spp, levels=x$spp)
s$relfrq <- unlist(c(frq[match(s$spp, names(frq))])) / NROW(d) * 100
s$logfrq <- log10(s$relfrq) # log relative frequency

### plot interaction coefficients (w bootstraps)
set_par_mercury(1)

ecole::set_par(1)
png(paste0('./fig/tmp/fig_05_interaction_coef.png'),
    wid=6.5, hei=4.5, uni='in', bg='transparent', res=700)

ggplot(s, aes(x = sen, y = spp, colour=posneg)) + geom_point(size=1) +
    scale_color_manual(values=c('#DD0000','#99999980','#0000DD')) +
    guides(color=F) + geom_errorbarh(aes(xmin=lwr, xmax=upr, height=0), size=0.3) +
    xlab('Effect size (bootstrapped coefficient value)') + ylab('') +
    facet_grid(response ~ predictor) + geom_vline(xintercept=0) + theme_classic() +
    theme(text=element_text(colour='black'),
          axis.text=element_text(family='Routed Gothic', colour='black', size=7),
          axis.text.y=element_text(colour=colvec(1:20,begin=0.1,end=0.8,alpha=1)),
          strip.background = element_blank(),
          strip.text = element_text(size=8), axis.title = element_text(size=8))
dev.off()


# edited to get rid of the gradient for species colors 
ecole::set_par(1)
png(paste0('./fig/tmp/fig_05_interaction_coef.png'),
    wid=6.5, hei=4.5, uni='in', bg='transparent', res=700)

ggplot(s, aes(x = sen, y = spp, colour=posneg)) + geom_point(size=1) +
  scale_color_manual(values=c('#DD0000','#99999980','#0000DD')) +
  guides(color=F) + geom_errorbarh(aes(xmin=lwr, xmax=upr, height=0), size=0.3) +
  xlab('Effect size (bootstrapped coefficient value)') + ylab('') +
  facet_grid(response ~ predictor) + geom_vline(xintercept=0) + theme_classic() +
  theme(text=element_text(colour='black'),
        axis.text=element_text(family='Arial', colour='black', size=7),
        axis.text.y=element_text(face="italic", size=8),
        strip.background = element_blank(),
        strip.text = element_text(size=8), axis.title = element_text(size=8))
dev.off()

### plot interaction *surfaces*
sortstress <- rev(sort(tapply(p$stress, p$spp, FUN=mean)))
nm_top <- nm[1:n_top]
nm_top <- names(rev(sort(sortstress[names(sortstress) %in% nm_top])))
`plot_surface` <- function(y, x1, x2, ...) {
    x1_rng <- range(x1, na.rm=T)
    x2_rng <- range(x2, na.rm=T)
    m      <- lm(y ~ x1 * x2, data=p)
    x1_seq <- seq(x1_rng[1], x1_rng[2], len=77)
    x2_seq <- seq(x2_rng[1], x2_rng[2], len=77)
    z      <- outer(x1_seq, x2_seq, FUN = function(a, b)
        predict(m, data.frame(x1 = a, x2 = b)))
    image(x1_seq, x2_seq, z, main = '', col = viridis::inferno(99),
          cex.lab=0.5, cex.axis=0.3, ...)
    contour(x1_seq, x2_seq, z, add=T, drawlabels=F)
}
setDF(p)

setwd("~/Documents/WSU/R Scripts")
png(paste0('./fig/tmp/fig_06_intxn_surfaces.png'),
    wid=3.5, hei=15.5, uni='in', bg='transparent', res=400)
set_par_mercury(1)
ecole::set_par(1)
par(mfcol=c(1, 4), mar=c(2,2,0,0), oma=c(0,2.1,2,0), mgp=c(1,.1,0), tcl=0)
lapply(nm_8, function(i) {
    p_i <- p[p$spp==i,]
    plot_surface(y  = p_i$rgr,
                 x1 = p_i$cmdlog,
                 x2 = p_i$plt_ses_pc1,
                 xlab = 'CMD',
                 ylab = 'Func neigh')
    mtext(text = i, side = 2, line=2, cex=0.3, las=3, font=2)
    if(i == 'quercus_kelloggii')
        mtext(text = 'Growth', side = 3, line=1, cex=0.5) })
lapply(nm_8, function(i) {
    p_i <- p[p$spp==i,]
    plot_surface(y  = p_i$rgr,
                 x1 = p_i$firelog,
                 x2 = p_i$plt_ses_pc1,
                 xlab = 'Fire',
                 ylab = 'Func neigh')
    if(i == 'quercus_kelloggii')
        mtext(text = 'Growth', side = 3, line=1, cex=0.5) })
lapply(nm_8, function(i) {
    p_i <- p[p$spp==i,]
    plot_surface(y  = p_i$surv,
                 x1 = p_i$cmdlog,
                 x2 = p_i$plt_ses_pc1,
                 xlab = 'CMD',
                 ylab = 'Func neigh')
    if(i == 'quercus_kelloggii')
        mtext(text = 'Survival', side = 3, line=1, cex=0.5) })
lapply(nm_8, function(i) {
    p_i <- p[p$spp==i,]
    plot_surface(y  = p_i$surv,
                 x1 = p_i$firelog,
                 x2 = p_i$plt_ses_pc1,
                 xlab = 'Fire',
                 ylab = 'Func neigh')
    if(i == 'quercus_kelloggii')
        mtext(text = 'Survival', side = 3, line=1, cex=0.5) })
dev.off()



### plot interaction *surfaces* for top-8 species *individually* in loop
dev.off()


nm_8 <- c(
  'chrysolepis_chrysophylla',
  'pinus_lambertiana',
  'pinus_monticola',
  'thuja_plicata'
)
x <- p[p$spp %in% nm_8,] # subset
# plot and save
lapply(1:4, function(i) {
    # subset species
    i <- 1
    p_i <- x[x$spp==nm_8[i],]
    # plotting parameters
    fnm <- paste0('fig_', sprintf('%02d',i),
                  '_surface_', nm_8[i], '.png')
    setwd("~/Documents/WSU/R Scripts")
    png(fnm, wid=6.5, hei=2.25, uni='in', bg='transparent', res=700)
    ecole::set_par(1)
    par(mfcol=c(1,4), mar=c(1.7,1.7,0,1), oma=c(0.1,0.3,1,0.1),
        mgp=c(1,0.1,0), tcl=0)
    # survival ~ cmd*ncis
    plot_surface(y  = p_i$rgr,
                 x1 = p_i$cmdlog,
                 x2 = p_i$plt_ses_pc1,
                 xlab = 'CMD',
                 ylab = 'Func neigh')
    mtext(text = 'Growth', side = 3, line=1, cex=0.5)
    mtext(text = nm_8[i], side = 3, line=2, cex=0.75)
    # survival ~ cmd*ncis
    plot_surface(y  = p_i$rgr,
                 x1 = p_i$firelog,
                 x2 = p_i$plt_ses_pc1,
                 xlab = 'Fire',
                 ylab = 'Func neigh')
    mtext(text = 'Growth', side = 3, line=1, cex=0.5)
    # survival ~ cmd*ncis
    plot_surface(y  = p_i$surv,
                 x1 = p_i$cmdlog,
                 x2 = p_i$plt_ses_pc1,
                 xlab = 'CMD',
                 ylab = 'Func neigh')
    mtext(text = 'Survival', side = 3, line=1, cex=0.5)
    # survival ~ cmd*ncis
    plot_surface(y  = p_i$surv,
                 x1 = p_i$firelog,
                 x2 = p_i$plt_ses_pc1,
                 xlab = 'Fire',
                 ylab = 'Func neigh')
    mtext(text = 'Survival', side = 3, line=1, cex=0.5)
})


####    END    ####