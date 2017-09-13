shinyServer(function(input, output) {
       
    output$plotlyHM <- renderPlotly({

        tlevel <- input$targetAggLevel
        dlevel <- input$diseaseAggLevel
        syms <- input$syms
        
        tcol <- switch(tlevel,
                       TDL = 'tdl',
                       Family = 'fam')
        dcol = switch(dlevel,
                      Level1 = 'level1',
                      Level2 = 'level2',
                      Level3 = 'level3',
                      Raw = 'did')

        ## dgrid <- disease %>%
        ##     dplyr::filter(dtype == "JensenLab Text Mining" & !is.na(zscore)) %>%
        ##     select(target_id, did, zscore) %>%
        ##     mutate(z100 = scales::rescale(zscore, to=c(0,100)))
        ## dgrid <- merge(dgrid, targets[,c('id','tdl','fam')], by.x='target_id', by.y='id')
        ## dgrid <- merge(dgrid, dop, by.x='did', by.y='doid')
        ## dgrid <- merge(dgrid, proteins[,c('id','uniprot','geneid','sym','chr')], by.x='target_id', by.y='id')
        ## dgrid$fam[is.na(dgrid$fam)] <- 'Unspecified'

        ## tmp <- do.call(rbind, mclapply(parent.paths, function(pd) {
        ##     if (nrow(pd) == 1) return(NULL)
        ##     l1 <- pd$id[nrow(pd)-1]
        ##     l2 <- NA
        ##     if (nrow(pd) > 2) l2 <- pd$id[nrow(pd)-2]
        ##     l3 <- NA
        ##     if (nrow(pd) > 3) l3 <- pd$id[nrow(pd)-3]
        ##     data.frame(did=pd$id[1], level1=l1, level2=l2, level3=l3)
        ## }, mc.cores=6))
        ## dgrid <- merge(dgrid, tmp, by='did')
        
        syms <- sapply(str_split(syms,",")[[1]], str_trim)
        if (syms == "" || length(syms) == 0) {
            m <- dgrid %>%
                group_by_(tcol,dcol) %>%
                summarize(value = median(zscore))
            names(m) <- c('tcol', 'dcol', 'value')
        } else {
            m <- dgrid %>%
                dplyr::filter(sym %in% syms) %>%
                group_by_('sym',dcol) %>%
                summarize(value = median(zscore))
            names(m) <- c('tcol', 'dcol', 'value')
        }

        
        ## Convert to matrix
        tmp <- dcast(m, dcol ~ tcol) %>% dplyr::filter(!is.na(dcol))
        rns <- tmp$dcol
        cns <- names(tmp)[-1]
        tmp <- as.matrix(tmp[,-1])
        
        row.names(tmp) <- do$name[ do$id %in% rns ]# subset(dns, did == rns)$name #rns
        colnames(tmp) <- cns

        #myClust <- hclust(dist(tmp))
        #tmp <- tmp[myClust$order, rev(myClust$order)]

        txtmat <- matrix("", nrow=nrow(tmp), ncol=ncol(tmp))
        ## for (i in 1:nrow(tmp)){
        ##     for (j in 1:ncol(tmp)) {
        ##         txtmat[i,j] <- sprintf("%s: %s<br>%s<br>Z: %3.2f",
        ##                                row.names(tmp)[i],
        ##                                subset(disease, id == row.names(tmp)[i])$name[1],
        ##                                colnames(tmp)[j],
        ##                                tmp[i,j])
        ##     }
        ## }
        p <- plot_ly(z = tmp,
                x = colnames(tmp),
                y = row.names(tmp),
                type = "heatmap",
                showscale = TRUE
                #hoverinfo = 'text',
                #text = txtmat
                )
        
        layout(p, 
               title = "",
               autosize = FALSE,
               width = 850,
               height = 1600,
               margin = list(l = 200, r = 50, b = 250, t = 50, pad = 4),
               font = list(size = 12),
               xaxis = list(title = ""),
               yaxis = list(title = "")
        )
    })
      
})
