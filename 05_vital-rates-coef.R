######################################################################
#
#  PNW-TREE vital rates -- linear model coefficients and sensitivities
#
#     Rob Smith, phytomosaic@gmail.com, 28 Aug 2020
#
##      GNU General Public License, Version 3.0    ###################

### preamble
rm(list=ls())
remotes::install_github('phytomosaic/ecole')
require(ecole)
require(data.table)
require(sensitivity)
require(ggplot2)

### plotting specs
x1 <- expression(log[10] ~ CMD ~ (mm ~ y^{-1}))
x2 <- expression(Fire ~ probability^{1/3})
x3 <- 'Composite stress'
y2 <- 'Growth (RGR)'
y3 <- 'Survival'
y4 <- 'Crowding (NCI)'
y5 <- 'Func neigh (PC1)' # NCIS
y6 <- 'Func neigh (PC2)' # NCIS

### load data
load('./data/pnw_tree.rda', verbose=T)
d <- pnw_tree  ;  rm(pnw_tree)
### load functional neighborhood SES per individual
#     NCIS = func neighborhood values =
#        pairwise trait deviations, multiplied by DBH^2, divided by dist^2
#     SES = deviation from random expectation, scaled by SD of random draws
#       problem: should scale and center traits first
#       problem: SES involving congeners always = 0 since traits identical
d[, c('pc1','pc2') := NULL] # rm climate PCs (avoid name conflict w trait PCs)
fnm <- c('./data_raw/func_neigh_ses/ses_PC1_alltraits.csv',
         './data_raw/func_neigh_ses/ses_PC2_alltraits.csv')
ses <- lapply(fnm, fread, colClasses = c('character', 'numeric'))
ses <- Reduce(merge,ses)
names(ses) <- gsub('stand_size_', '', tolower(names(ses)))
names(ses)[names(ses)=='focal'] <- 'tre_cn'
d   <- merge(d, ses, by='tre_cn')
rm(ses)

### exclude juveniles < 5.0 inches (12.7 cm) per GRM and Fenn et al. (2020)
i <- !((!is.na(d$dia_y1) & d$dia_y1 < 5) |
         (!is.na(d$dia_y2) & d$dia_y2 < 5))
d <- d[i,]

### transformations
d$pai     <- d$dia_y2 - d$dia_y1  # periodic annual increment (untransformed)
d$dia_y1  <- log10(1 + d$dia_y1)  # log10-transform
d$dia_y2  <- log10(1 + d$dia_y2)  # log10-transform
d$annincr <- d$dia_y2 - d$dia_y1  # annual increment (transformed)
d$rgr     <- d$annincr / d$dia_y1 # relative growth rate (relative to y1)

### summarize means of SPP in PLOT (SPP/PLOT = sample unit)
p <- d[, list(
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
is.na(p) <- is.na(p)

### frequencies
frq <- rev(sort(table(d$spp)))
nm  <- names(frq)

### calculate 'composite stress' gradient (sum of fire/drought percentiles)
`ptile` <- function(x) ecdf(x)(x)
a  <- cbind(fire=ptile(p$firelog), cmd=ptile(p$cmdlog)) # percentiles
p$stress <- rowMeans(a, na.rm=T) ; rm(a, ptile) # max value = 1 is intuitive

### END data prep ##############################################################


### BEGIN mapping ##############################################################

### map enviro gradients
png(filename='./fig/tmp/fig_00_stressgradient.png', wid=6.5,
    hei=6.5, uni='in', bg='transparent', res=700)
set_par_mercury(4)
cc <- ecole::colvec(p$stress, begin=0.1, end=0.95, alpha=0.80)
plot(p$lon, p$lat, pch=16, cex=0.4, xlab='', ylab='', xlim=c(-125,-115),
     col=ecole::colvec(p$cmdlog, begin=0.1, end=0.95, alpha=0.80), asp=1.1)
add_text(0.50,0.1,x1)
plot(p$lon, p$lat, pch=16, cex=0.4, xlab='', ylab='', xlim=c(-125,-115),
     col=ecole::colvec(p$firelog, begin=0.1, end=0.95, alpha=0.80), asp=1.1)
add_text(0.55,0.1,x2)
plot(p$cmdlog, p$firelog, pch=16, cex=0.4, col=cc, ylab=x2, xlab=x1)
add_text(0.55,0.1,x3)
plot(p$lon, p$lat, pch=16, cex=0.4, xlab='', ylab='', xlim=c(-125,-115),
     col=cc, asp=1.1)
add_text(0.55,0.1,x3)
dev.off()

### map tree responses
png(filename='./fig/tmp/fig_00_map_crossspecies_responses.png', wid=9.5,
    hei=6.5, uni='in', bg='transparent', res=700)
set_par_mercury(5)
cc <- ecole::colvec(p$stress, begin=0.1, end=0.95, alpha=0.80)
plot(p$lon, p$lat, pch=16, cex=0.4, xlab='', ylab='', xlim=c(-125,-115),
     col=ecole::colvec(p$rgr, begin=0.1, end=0.95, alpha=0.80), asp=1.1)
add_text(0.50,0.1,y2)
plot(p$lon, p$lat, pch=16, cex=0.4, xlab='', ylab='', xlim=c(-125,-115),
     col=ecole::colvec(p$surv, begin=0.1, end=0.95, alpha=0.80), asp=1.1)
add_text(0.55,0.1,y3)
plot(p$lon, p$lat, pch=16, cex=0.4, xlab='', ylab='', xlim=c(-125,-115),
     col=ecole::colvec(p$plt_nci, begin=0.1, end=0.95, alpha=0.80), asp=1.1)
add_text(0.55,0.1,y4)
plot(p$lon, p$lat, pch=16, cex=0.4, xlab='', ylab='', xlim=c(-125,-115),
     col=ecole::colvec(p$plt_ses_pc1, begin=0.1, end=0.95, alpha=0.80), asp=1.1)
add_text(0.55,0.1,y5)
plot(p$lon, p$lat, pch=16, cex=0.4, xlab='', ylab='', xlim=c(-125,-115),
     col=ecole::colvec(p$plt_ses_pc2, begin=0.1, end=0.95, alpha=0.80), asp=1.1)
add_text(0.55,0.1,y6)
dev.off()
### END mapping ################################################################


### BEGIN models ###############################################################

###  CROSS-SPECIES models --- univariate for CMD or FIRE
`plot_lm` <- function (x, y, col='#00000040', lcol='#FF0000BF', cex=0.8, pch=16,
                       xlab=NULL, ylab=NULL, args.func=list(), r2=T, ...) {
  if (is.null(xlab)) xlab <- deparse(substitute(x))
  if (is.null(ylab)) ylab <- deparse(substitute(y))
  plot(x=x, y=y, col=col, pch=pch, cex=cex, xlab=xlab, ylab=ylab, ...)
  f <- do.call(stats::lm, c(list(formula = y ~ x), args.func))
  abline(f, col = lcol, lwd = 2)
  if (r2) ecole::add_text(0.05, 0.90,
                          paste0('Adj-R2 = ',round(summary(f)$adj.r.squared,2)))
}
png(filename='./fig/tmp/fig_00_lm_across_spp.png',
    wid=10.5, hei=4.25, uni='in', bg='transparent', res=500)
set_par_mercury(10)  ;  par(mfcol=c(2,5))
cc <- ecole::colvec(p$stress, begin=0.1, end=0.95, alpha=0.80)
plot_lm(p$cmdlog,  p$rgr, cex=0.2, xlab=x1, ylab=y2, col=cc)
plot_lm(p$firelog, p$rgr, cex=0.2, xlab=x2, ylab=y2, col=cc)
plot_lm(p$cmdlog,  p$surv, cex=0.2, xlab=x1, ylab=y3, col=cc)
plot_lm(p$firelog, p$surv, cex=0.2, xlab=x2, ylab=y3, col=cc)
plot_lm(p$cmdlog,  p$plt_nci, cex=0.2, xlab=x1, ylab=y4, col=cc)
plot_lm(p$firelog, p$plt_nci, cex=0.2, xlab=x2, ylab=y4, col=cc)
plot_lm(p$cmdlog,  p$plt_ses_pc1, cex=0.2, xlab=x1, ylab=y5, col=cc)
plot_lm(p$firelog, p$plt_ses_pc1, cex=0.2, xlab=x2, ylab=y5, col=cc)
plot_lm(p$cmdlog,  p$plt_ses_pc2, cex=0.2, xlab=x1, ylab=y6, col=cc)
plot_lm(p$firelog, p$plt_ses_pc2, cex=0.2, xlab=x2, ylab=y6, col=cc)
dev.off()

###  WITHIN-SPECIES models --- all predictors additively
set.seed(88)
n_boot <- 99 # ! ! ! TIMEWARN ! ! !
n_top  <- 20 # take only top-20 most frequent spp
j      <- c('rgr','surv') # responses
j_pred <- c('cmdlog','firelog','plt_nci','plt_ses_pc1','plt_ses_pc2') # predictors
# 300 rows / 20 species / 2 response vars / 6 predictors
`sens` <- function(X, y) {
  s <- sensitivity::src(X, y, nboot=n_boot)$SRC
  c(sen = s[['original']] - s[['bias']],
    lwr = s[['min. c.i.']], upr = s[['max. c.i.']])
}
s <- lapply(j, function(yvar) { # apply for each combn of species/pred/resp
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
s$metric <- rep(rep(c('sen','lwr','upr'), ea=length(j_pred)), times=length(j))
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
                               c('CMD','FIRE','CROWDING',
                                 'Func neigh PC1','Func neigh PC2'))
s$predictor <- factor(s$predictor,
                      levels=c('CMD','FIRE','CROWDING',
                               'Func neigh PC1','Func neigh PC2'))
s$posneg    <- factor(ifelse(s$lwr > 0.0, 1, ifelse(s$upr < (-0.0), -1, 0)))
s$bootrng   <- s$upr - s$lwr  # <<--- sensitivity
### sort species by STRESS
x <- subset(p, spp %in% nm[1:20]) # subset top-20 species
x <- x[, lapply(.SD, mean, na.rm=T), by=spp, .SDcols='stress']
setorder(x, stress, spp) # order by STRESS
s$spp    <- factor(s$spp, levels=x$spp)
s$relfrq <- unlist(c(frq[match(s$spp, names(frq))])) / NROW(d) * 100
s$logfrq <- log10(s$relfrq) # log relative frequency
### EFFECT SIZES from standardized regression coefficients (w bootstraps)
set_par_mercury(1)
png(paste0('./fig/tmp/fig_04_effectsize.png'),
    wid=6.5, hei=4.5, uni='in', bg='transparent', res=700)
ggplot(s, aes(x = sen, y = spp, colour=posneg)) + geom_point(size=1) +
  scale_color_manual(values=c('#DD0000','#99999980','#0000DD')) +
  guides(color=F) + geom_errorbarh(aes(xmin=lwr, xmax=upr, height=0), size=0.3) +
  xlab('Effect size (bootstrapped coefficient value)') + ylab('') +
  facet_grid(response ~ predictor) + geom_vline(xintercept=0) + theme_classic() +
  theme(text=element_text(family='Routed Gothic', colour='black'),
        axis.text=element_text(family='Routed Gothic', colour='black', size=7),
        axis.text.y=element_text(colour=colvec(1:20,begin=0.1,end=0.8,alpha=1)),
        strip.background = element_blank(),
        strip.text = element_text(size=8), axis.title = element_text(size=8))
dev.off()
### SENSITVITY (same species order as above)
png(paste0('./fig/tmp/fig_05_sensitivity.png'),
    wid=9.5, hei=6.5, uni='in', bg='transparent', res=700)
ggplot(s, aes(x = bootrng, y = spp)) +
  geom_segment(aes(x=0,xend=bootrng,y=spp,yend=spp),color='black',size=0.01) +
  geom_point() + scale_color_manual(values=c('#DD0000','#99999980','#0000DD')) +
  guides(color=F) + xlab('Sensitivity (bootstrap range)') + ylab('') +
  facet_grid(response ~ predictor) + geom_vline(xintercept=0) + theme_classic()+
  theme(text=element_text(family='Routed Gothic', colour='black'),
        axis.text=element_text(family='Routed Gothic', colour='black'),
        axis.text.y=element_text(colour=colvec(1:20,begin=0.1,end=0.8,alpha=1)),
        strip.background = element_blank(),
        strip.text = element_text(size=12), axis.title = element_text(size=12))
dev.off()

####    END    ####