######################################################################
#
#  PNW-TREE vital rates -- get phylogeny for PNW trees
#
#    Rob Smith, phytomosaic@gmail.com, 02 April 2020
#
##      GNU General Public License, Version 3.0    ###################

### preamble
rm(list=ls())
require(ecole)
require(viridis)
require(brranching)

devtools::install_github("jinyizju/V.PhyloMaker")
require(V.PhyloMaker)

setwd("~/Documents/WSU/R Scripts")
load('pnw_tree-3.rda', verbose=T)
d <- pnw_tree  ;  rm(pnw_tree)

### get family / genus / species
(tax <- sort(unique(d$spp)))
tax  <- brranching::phylomatic_names(tax, db='apg')
tax[gsub('\\/.*','',tax) == 'NA']  # expect 0
tax  <- data.frame(do.call(rbind, strsplit(tax,'/')))
tax  <- data.frame(species = tax[,2],
                   genus = sub('\\_.*', '', tax[,2]),
                   family = tax[,1])

# change name format to remove underscore and capitalize genus 
tax$species <- gsub("_", " ", tax$species)
tax$species <- gsub("^(\\w)(\\w+)", "\\U\\1\\L\\2",
              tax$species, perl = TRUE)

### query V.PhyloMaker ! ! ! TIMEWARN ! ! !  57 taxa ~30 sec
system.time(p <-V.PhyloMaker::phylo.maker(sp.list=tax,scenarios='S3'))
phy <- p$scenario.3
phy$tip.label <- tolower(phy$tip.label)
table(p$species.list$status)
if (all(tax$species %in% phy$tip.label) &&
    all(phy$tip.label %in% tax$species) ) {
     cat('all taxa appear in `phy`\n')
     rm(tax, p)
}
pnw_phy <- phy
save(pnw_phy, file='./pnw_phy.rda')
rm(pnw_phy)

pnw_phy$tip.label <- gsub("_", " ", pnw_phy$tip.label)
pnw_phy$tip.label <- gsub("^(\\w)(\\w+)", "\\U\\1\\L\\2",
                    pnw_phy$tip.label, perl = TRUE)

d$spp <- gsub("_", " ", d$spp)
d$spp <- gsub("^(\\w)(\\w+)", "\\U\\1\\L\\2",
                          d$spp, perl = TRUE)

### plot phylogeny, colored by frequency
frq  <- aggregate(d$spp, list(d$spp), length)
(frq <- frq[match(pnw_phy$tip.label, frq[,1]),]) # order to match phy
cc <- ecole::colvec(log10(frq$x), pal=viridis::viridis(99), alpha=1)
png('./fig_00_phylogeny.png', wid=6.5, hei=4.5, unit='in',
    bg='transparent', res=900)
plot(ape::ladderize(pnw_phy), no.margin=TRUE, cex=0.7, tip.col=cc,
     label.offset=10, edge.width = 2)
tiplabels(pch=19, col=cc, cex=0.8)
dev.off()

####    END    ####