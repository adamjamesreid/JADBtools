library(d3heatmap)
library(shiny)
library(DT)
require()

ui <- shinyUI(navbarPage(
    'CorTracks', id='nvpage',
    tabPanel(
        'Input',
        actionButton('selFilt', 'Select filtered'),
        actionButton('selPage', 'Select visible on page'),
        actionButton('selNone', 'Select none'),
        textOutput('nselected', inline = TRUE),
        actionButton('run', 'Correlate', class='btn, btn-success', icon = icon('paly')),
        tags$hr(),
        dataTableOutput('ex1')
    ),
    tabPanel(
        'Output',
        selectInput("palette", "Palette", c("Spectral", "YlOrRd", "RdYlBu", "Greens", "Blues")),
        checkboxInput("cluster", "Apply clustering", value = TRUE),
        checkboxInput("zlim", "-1 to 1 color key limits", value = FALSE),
        d3heatmapOutput("heatmap"),
        plotOutput('legend', height = 150)
))
    

)

server <- function(input, output, session) {
    output$ex1 <- renderDataTable(
        datatable(
            all  %>% 
                dplyr::select(Factor, Antibody, ExtractID, Crosslinker, Strain, Stage, ContactExpID, dateCreated, dateUpdated) %>% 
                mutate(dateCreated=as.Date(dateCreated),dateUpdated=as.Date(dateUpdated)),
            filter = 'top',
            extensions = c('ColReorder', 'ColVis'),
            plugins = 'natural',
            options = list(
                pageLength = 10, 
                autoWidth = TRUE,
                searchHighlight = TRUE,
                dom = 'CRlfrtip', 
                colReorder = list(realtime = TRUE)
            )
        ) %>% formatDate(c('dateCreated', 'dateUpdated'), 'toLocaleDateString')
    )
    
    reactcor <- eventReactive(input$run, {
        IDs <- unlist(all[input$ex1_rows_selected,]$ContactExpID)
        names(IDs) <- IDs
        
        IDs %>% sapply(getFilePath, format = 'bw', processing=processing) -> paths
        require(rtracklayer)
        bwfl <- BigWigFileList(unlist(paths))
        withProgress(
            message = 'Calculation in progress',
                     detail = 'This may take a while...', value = 0, {
                         lst <- lapply(bwfl, function(bwf) {
                             incProgress(1/length(bwfl), detail = basename(path(bwf)))
                             extarct_vector(bwf, size = seqlengths(bwf)/res)
                        })
                     })
        
        M <- do.call(cbind, lst)
        C <- cor(M)
        
        return(C)
    })
    
    output$heatmap <- renderD3heatmap({
        pal <- colorRampPalette(rev( 
                RColorBrewer::brewer.pal(
                    brewer.pal.info[input$palette,]$maxcolor,
                    input$palette
                ) 
        ))
        
        colorFunc <- col_bin(pal(512), bins = seq(-1,1,length.out = 512))
        d3heatmap(
            reactcor(), 
            colors = if(input$zlim) colorFunc else pal(512),
            symm = TRUE, 
            scale = 'none', 
            #theme='dark', 
            revC = TRUE,
            dendrogram = if (input$cluster) "both" else "none"
        )
    })
    
    output$legend <- renderPlot({
        pal <- colorRampPalette(rev( 
            RColorBrewer::brewer.pal(
                RColorBrewer::brewer.pal.info[input$palette,]$maxcolor,
                input$palette
            ) 
        ))
        colorFunc <- col_bin(pal(512), bins = seq(-1,1,length.out = 512))
        
        plot.new()
        fields::image.plot(
            if (input$zlim) matrix(c(-1,1, 0, 0), 2,2) else reactcor(), 
            legend.only = TRUE, horizontal = TRUE, col = pal(512),
            smallplot = c(.1,.9,.5,.7), legend.cex = 3.0
        )
    })
    
    output$nselected <- renderText({
        paste(length(input$ex1_rows_selected), 'experiemnts selected')
    })
    
    observe({
        if(input$run==0) return()
        updateNavbarPage(session, inputId = 'nvpage', selected = 'Output')
    })
    
    observe({
        if(input$selFilt==0) return()
        isolate({
            proxy <- dataTableProxy('ex1')
            selectRows(proxy, NULL)
            selectRows(proxy, input$ex1_rows_all)
        })
    })
    
    observe({
        if(input$selPage==0) return()
        isolate({
            proxy <- dataTableProxy('ex1')
            selectRows(proxy, NULL)
            selectRows(proxy, input$ex1_rows_current)
        })
    })
    
    observe({
        if(input$selNone==0) return()
        isolate({
            proxy <- dataTableProxy('ex1')
            selectRows(proxy, NULL)
        })
    })
    
}
runCor <- function(){
    shinyApp(ui, server)
}