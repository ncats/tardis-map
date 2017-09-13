library(RMySQL)
library(dplyr)
library(pracma)
library(parallel)
library(ineq)
##library(multidplyr)
db = dbConnect(MySQL(),
               user='root',
               password='NC@TS-IFX@',
               dbname='tcrd460',
               host='127.0.0.1',
               port=3307)
disease <- dbGetQuery(db, "select d.*, t.name as tname from disease d, target t where t.id = d.target_id")
targets <- dbGetQuery(db, "select * from target")
proteins <- dbGetQuery(db, "select * from protein")
dop <- dbGetQuery(db, "select * from do_parent")
do <- dbGetQuery(db, "select * from do")

## Construct DO hierarchy as tabular form to support
## aggregation at different levels
get.parents <- function(idlist) {
  last.id <- idlist[[ length(idlist) ]]

  p <- subset(dop, doid %in% last.id$id)
  if (nrow(p) == 0) return(idlist)
  else {
    ## Not entirely correct since a given term can have
    ## multiple parents
    p <- p[1,]
    idlist[[ length(idlist)+1 ]] <- list(id=p$parent, level=last.id$level+1)
    return(get.parents(idlist))
  }
}
parent.paths <- lapply(unique(do$id),
                       function(x) get.parents(list(list(id=x, level=0))))

## Convert path list to path data.frames, where first row is starting doid
## and subsequent rows are parents. Last row should be root
parent.paths <- lapply(parent.paths, function(pl) {
  tmp <- do.call(rbind, lapply(pl, function(pe) data.frame(id=pe$id, level=pe$level)))
  merge(tmp, do[, c('id','name')], by='id') %>% arrange(level)
})

save(disease,targets,dop,proteins,do,parent.paths,
     file='work.Rda')

##################################################################################

dgrid <- disease %>%
    dplyr::filter(dtype == "JensenLab Text Mining" & !is.na(zscore)) %>%
    select(target_id, did, zscore) %>%
    mutate(z100 = scales::rescale(zscore, to=c(0,100)))
dgrid <- merge(dgrid, targets[,c('id','tdl','fam')], by.x='target_id', by.y='id')
dgrid <- merge(dgrid, dop, by.x='did', by.y='doid')
dgrid <- merge(dgrid, proteins[,c('id','uniprot','geneid','sym''chr')], by.x='target_id', by.y='id')
dgrid$fam[is.na(dgrid$fam)] <- 'Unspecified'


## Make subsets of the target disease grid at various
## levels of aggregation
agg1 <- dgrid %>%
    group_by(tdl, parent) %>%
    summarize(value = median(zscore))
agg2 <- dgrid %>%
    group_by(fam, parent) %>%
    summarize(value = median(zscore))

ggplot(agg2, aes(fam, parent))+
    geom_tile(aes(fill=value), colour='white')+
    scale_fill_gradient(low = "white",
                        high = "steelblue")


## Compute Gini Coefficient as a way to measure selectivity
## based on http://pubs.acs.org/doi/abs/10.1021/jm070562u
## designed for kinase inhibitors. Use the Gini method from
## the ineq package

## gini.coeff <- function(df, all.dids, plot=FALSE) {
##     #all.dids <- all.dids[ !(all.dids %in% df$did) ]
##     #tmp1 <- data.frame(target_id = -1, did=all.dids, zscore=NA, z100 = 0)
##                                         #tmp2 <- rbind(tmp1, df)
##     df$z100 <-   scales::rescale(df$zscore, to=c(0,100))              
##     tmp <- df %>%
##         arrange(z100) %>%
##         mutate(serial = 1:n()) %>%
##         mutate(sfrac = serial/max(serial),
##                zfrac = z100 / sum(z100)) %>%
##         mutate(cs = cumsum(zfrac))
##     print(tmp)
##     print(ggplot(tmp, aes(sfrac,cs))+ geom_line())
##     1-2*trapz(tmp$sfrac, tmp$zfrac)
## }

## gini.coeff(subset(dgrid, target_id == 19144), unique(dgrid$did), TRUE)
## gini.coeff(subset(dgrid, target_id == 9342), unique(dgrid$did), TRUE)

tsum <- mclapply(unique(dgrid$target_id), function(x)  {
    dx <- subset(dgrid, target_id == x)
    data.frame(target_id = x, GC=Gini(dx$z100, corr=TRUE), TI=Theil(dx$z100, parameter=NULL))
}, mc.cores=4)
tsum <- do.call(rbind, tsum)

tsum <- merge(tsum, proteins[,c('id','sym','uniprot')], by.x='target_id', by.y='id')

ggplot(tsum, aes(x=GC))+
    geom_histogram(colour='black', fill='beige')+
    xlab("Gini Coefficient")+ylab("Number of targets")+
    xlim(c(0,1))+ggtitle("All Targets")
ggsave("target-gini-distribution.pdf")

tmp <- merge(tsum, targets, by.x='target_id', by.y='id')
ggplot(tmp, aes(x=GC))+
    geom_histogram(colour='black', fill='beige')+
    xlab("Gini Coefficient")+ylab("Number of targets")+
    xlim(c(0,1))+facet_wrap(~tdl, scale='free')+
    theme(strip.text.x = element_text(size=14,face='bold'))
ggsave("target-gini-distribution-tdl.pdf")

tmp$fam[ is.na(tmp$fam) ] <- 'Unspecified'
ggplot(tmp, aes(x=GC))+
    geom_histogram(colour='black', fill='beige')+
    xlab("Gini Coefficient")+ylab("Number of targets")+
    xlim(c(0,1))+facet_wrap(~fam, scale='free')+
    theme(strip.text.x = element_text(size=12,face='bold'))
ggsave("target-gini-distribution-fam.pdf")


###  WHAMM, WISP2, TP53
tids <- subset(tsum, sym %in% c('WHAMM','WISP2','TP53'))$target_id
tmp <- subset(dgrid, target_id == 9342) %>%
    arrange(z100) %>%
    mutate(serial = 1:n()) %>%
    mutate(sfrac = serial/max(serial),
           zfrac = z100 / sum(z100)) %>%
    mutate(cs = cumsum(zfrac))
trapz(tmp$sfrac, tmp$zfrac)
ggplot(tmp, aes(sfrac,cs, colour=target_id))+ geom_line()
