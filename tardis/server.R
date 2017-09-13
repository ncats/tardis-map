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
                      Raw = 'did')

        dgrid <- disease %>%
            dplyr::filter(dtype == "JensenLab Text Mining" & !is.na(zscore)) %>%
            select(target_id, did, zscore) %>%
            mutate(z100 = scales::rescale(zscore, to=c(0,100)))
        dgrid <- merge(dgrid, targets[,c('id','tdl','fam')], by.x='target_id', by.y='id')
        dgrid <- merge(dgrid, dop, by.x='did', by.y='doid')
        dgrid <- merge(dgrid, proteins[,c('id','uniprot','geneid','sym','chr')], by.x='target_id', by.y='id')
        dgrid$fam[is.na(dgrid$fam)] <- 'Unspecified'

        print(head(dgrid))
        
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
        tmp <- dcast(m, dcol ~ tcol)
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
