

#Load in Libraries and Functions Used

library(readr)
library(stringr)
library(webshot)
library(plyr)
library(shiny)
library(shinythemes)
library(searchable)
library(DBI)
library(shinyjs)
library(tidyverse)
library(plotly)
library(shinydashboard)
library(RColorBrewer)
library(rcytoscapejs)
source("functions.R")

server <- function(input, output, session) {
  # data <- reactive ({ #broken - could never get reactive code to work
  #   getdata(input)
  # })
  
  ploter <- function() {
    # function used to generate plot
    plot_checker <- plot_select(input)
    GENE <- toupper(input$GENE)
    GENE2 <- toupper(get_gene2(input))
    if (GENE2 == "")
    {
      ga_collect_event(event_category = "gene search2",
                       event_action = paste(GENE, " ", GENE2, sep = ""))
      GENE2 = GENE
    }
    data <- getdata(GENE, GENE2, input$libchoose, input) #Getting Data
    validate(need((nrow(data) > 0),#Used to make sure gene is found
                  paste(
                    "One of input genes not found in selected library or no tissue selected"
                  )
    )) 
    ga_collect_event(event_category = "gene search", event_action = GENE)
    data <- arrange(data, BF) #arranging data for plot
    names = str_split_fixed(data$CCLE_CELL, "_", 2) #Spliting CCLE CELL up, so that tissue is used for name
    data$Name = names[, 1]
    data <- transform(data, Name = reorder(Name,-BF)) #Reorder for BF
    data$Type = names[, 2] #Data Type, could use primary disease from database as well.
    data <- droplevels(data)
    Expression_tag = 'Log(CCLE RPKM): ' #different Exp Methods
    if (input$libchoose == 7)
    {
      Expression_tag = 'Log(TPM): '
    }
    if (input$libchoose == 6)
    {
      Expression_tag = 'Log(FPKM): '
    }
    ay <- list(
      #secondary y axis for Expression
      tickfont = list(color = "violet"),
      overlaying = "y",
      side = "right",
      title = "Expression",
      showgrid = FALSE,
      zeroline = FALSE,
      showline = TRUE
    )
    ay2 <- list(
      #secondary y axis for Mutations
      tickfont = list(color = "black"),
      overlaying = "y",
      side = "right",
      title = "",
      tickformat = ',d',
      showgrid = FALSE,
      zeroline = FALSE,
      showline = TRUE#, range = c(-0.15,1.15)
    )
    ay3 <- list(
      #secondary y axis for CN
      tickfont = list(color = "black"),
      overlaying = "y",
      side = "right",
      title = "Copy Number",
      tickformat = ',d',
      showgrid = FALSE,
      zeroline = FALSE,
      showline = TRUE#, range = c(-0.15,1.15)
    )
    ax <- list(
      #top x axis, just used for aesthetics
      overlaying = "x",
      side = "top",
      showgrid = FALSE,
      zeroline = FALSE,
      showline = TRUE,
      showticklabels = FALSE,
      showlabels = FALSE
    )
    
    #Initial Plotting of BFs
    p <-
      plot_ly() %>% add_markers(
        x = ~ data$Name,
        y = ~ data$BF,
        name = "BF",
        hoverinfo = 'text',
        marker = list(
          size = 15,
          opacity = 0.65,
          color = data$colors
        ),
        text = ~ paste(
          'BF: ',
          data$BF,
          '\n',
          'Cell Line:',
          data$Name,
          '\n',
          'Type:',
          data$Type
        )
      )
    #Plotting TSG Threshold if BF are present (plot is skewed otherwise)
    if (min(na.omit(data$BF)) < -49 && input$libchoose != 8) #Need custom range for 8
    {
      p <-
        add_lines(
          p,
          x = ~ data$Name,
          y = -50,
          name = "TSG Behavior Threshold",
          line = list(color = 'blue'),
          hoverinfo = 'text',
          text = ~ paste('Threshold for Tumor Suppressor Behavior, BF = ',-50)
        )
    }
    if (5 %in% plot_checker) {
      #Plot code for the Secondary BF
      data$BF2 = round(data$BF2, digits = 3)
      p <-
        add_markers(
          p,
          x = ~ data$Name,
          y = ~ data$BF2,
          name = "Compared Gene BFs",
          hoverinfo = 'text',
          marker = list(
            size = 8,
            opacity = 0.65,
            color = data$colors,
            symbol = 2
          ),
          text = ~ paste(
            '2nd Gene BF: ',
            data$BF2,
            '\n',
            'Cell Line:',
            data$Name,
            '\n',
            'Type:',
            data$Type,
            '\n'
          )
        )
      
      graphtitle = paste(GENE, " Essentiality with ", GENE2, " Essentiality", sep = "")
      #Adding Layout aesthetics
      p <- layout(
        p,
        title = graphtitle,
        xaxis2 = ax,
        yaxis2 = list(
          overlaying = "y",
          side = "right",
          showgrid = FALSE,
          zeroline = FALSE,
          showline = TRUE,
          showticklabels = FALSE,
          showlabels = FALSE
        ),
        yaxis = list(
          title = "Quantile Normalized BF",
          showgrid = TRUE,
          showline = TRUE,
          zeroline = FALSE
        ),
        xaxis = list(
          title = "Cell Lines",
          tickangle = 45,
          showgrid = TRUE,
          showline = TRUE,
          tickmode = 'array'
        ),
        margin = list(b = 100, r = 100)
      ) %>% config(displayModeBar = FALSE) %>%  layout(showlegend = TRUE)
    }
    
    if (6 %in% plot_checker) {
      #Plot code for Expression Metrics
      data$EXP = log(data$EXP + 0.001) #fixing EXP
      data$EXP[data$EXP < -2] = -2
      data$EXP = round(data$EXP, digits = 3)
      p <-
        add_markers(
          p,
          x = ~ data$Name,
          y = ~ data$EXP,
          name = "Expression",
          yaxis = "y2",
          marker = list(color = 'violet', opacity = .75),
          hoverinfo = 'text',
          text = ~ paste(
            Expression_tag,
            data$EXP,
            '\n',
            'Cell Line:',
            data$Name,
            '\n',
            'Type:',
            data$Type
          )
        )
      
      graphtitle = paste(GENE, " Essentiality with ", GENE2, " Expression", sep = "")
      
      #Adding Layout aesthetics
      p <- layout(
        p,
        title = graphtitle,
        xaxis2 = ax,
        yaxis2 = ay,
        yaxis = list(
          title = "Quantile Normalized BF",
          showgrid = TRUE,
          showline = TRUE,
          zeroline = FALSE
        ),
        xaxis = list(
          title = "Cell Lines",
          tickangle = 45,
          showgrid = TRUE,
          showline = TRUE,
          tickmode = 'array'
        ),
        margin = list(b = 100, r = 100)
      ) %>% config(displayModeBar = FALSE) %>%  layout(showlegend = TRUE)
    }
    
    if (7 %in% plot_checker) { #Plot code for mutations
      data$MUT[data$MUT != "Silent or Not Detected"] = "Mutant"
      data$MUT[data$MUT == "Silent or Not Detected"] = "Not-Detected \n or Silent"
      
      p <-
        add_markers(
          p,
          x = ~ data$Name,
          y = ~ data$MUT,
          name = "Mutations",
          yaxis = "y2",
          marker = list(
            size = 8,
            opacity = 0.33,
            color = 'black',
            symbol = 1
          ),
          hoverinfo = 'text',
          text = ~ paste(
            'Cell Line:',
            data$Name,
            '\n',
            'Type:',
            data$Type,
            '\n',
            'Mutation:',
            data$MUT
          )
        )
      graphtitle = paste(GENE, " Essentiality with ", GENE2, " Mutations", sep = "")
      p <- layout(
        p,
        title = graphtitle,
        xaxis2 = ax,
        yaxis2 = ay2,
        yaxis = list(
          title = "Quantile Normalized BF",
          showgrid = TRUE,
          showline = TRUE,
          zeroline = FALSE
        ),
        xaxis = list(
          title = "Cell Lines",
          tickangle = 45,
          showgrid = TRUE,
          showline = TRUE,
          tickmode = 'array'
        ),
        margin = list(b = 100, r = 100)
      ) %>% config(displayModeBar = FALSE) %>%  layout(showlegend = TRUE)
    }
    
    if (8 %in% plot_checker) { #CN Plot code option
      data$CN[data$CN < -6] = -6
      data$CN = round(data$CN, digits = 3)
      p <-
        add_markers(
          p,
          x = ~ data$Name,
          y = ~ data$CN,
          name = "Copy Number",
          yaxis = "y2",
          marker = list(color = 'orange', opacity = .75),
          hoverinfo = 'text',
          text = ~ paste(
            "Copy Number",
            data$CN,
            '\n',
            'Cell Line:',
            data$Name,
            '\n',
            'Type:',
            data$Type
          )
        )
      graphtitle = paste(GENE, " Essentiality with ", GENE2, " Copy Number", sep = "")
      p <-
        layout(
          p,
          title = graphtitle,
          xaxis2 = ax,
          yaxis2 = ay3,
          yaxis = list(
            title = "Quantile Normalized BF",
            showgrid = TRUE,
            showline = TRUE,
            zeroline = FALSE
          ),
          xaxis = list(
            title = "Cell Lines",
            tickangle = 45,
            showgrid = TRUE,
            showline = TRUE,
            tickmode = 'array'
          ),
          margin = list(b = 100, r = 100)
        ) %>% config(displayModeBar = FALSE) %>% layout(showlegend = TRUE)
    }
    
    p <- #Adding Essentiality Threshold Line, and other general primary y-axis aesthetics
      add_lines(
        p,
        x = ~ data$Name,
        y = 5,
        name = "Essentiality Threshold",
        line = list(color = 'red'),
        hoverinfo = 'text',
        text = ~ paste('Essentiality Threshold, BF = ', 5)
      )
    p <- layout(
      p,
      xaxis2 = ax,
      yaxis = list(
        title = "Quantile Normalized BF",
        showgrid = TRUE,
        showline = TRUE,
        zeroline = FALSE
      ),
      xaxis = list(
        title = "Cell Lines",
        tickangle = 45,
        showgrid = TRUE,
        showline = TRUE,
        tickmode = 'array'
      ),
      margin = list(b = 100, r = 100)
    ) %>%
      config(displayModeBar = FALSE) %>%
      layout(
        showlegend = TRUE,
        legend = list(
          orientation = 'h',
          xanchor = "center",
          x = 0.5,
          yanchor = "top",
          y = -0.3
        )
      )
    
    hide(id = "loading-content", #This is code for the loading page - ends when initial plot loads
         anim = TRUE,
         animType = "fade") #end resulting plotly image
    return(p)
  }
  
  output$plot <- renderPlotly({
    #rendering the output plot
    checkGroup <- check_groups(input)
    validate(need((length(checkGroup) > 0), "No tissue selected")) #Makes sure tissue selected
    
    ploter()
  })
  
  output$cytoplot <- renderRcytoscapejs({ #Code for the Cytoscape Plot on Avana Q3
    if (input$libchoose == 1) {
      # load in network if selected is network
      GENE <- toupper(input$GENE)
      clust <- find_cluster(GENE)
      validate(need((nrow(clust) > 0), "Gene not associated to other genes in Avana data"))
      geneset = unique(toupper(clust$GENE))
      network <- get_coess_network(GENE, geneset)
      if (nrow(network) > 40)
      {
        network <- network %>% filter(GENE1 == GENE | GENE2 == GENE)
        geneset = unique(c(network$GENE1,network$GENE2))
      }
      if (nrow(network) > 40)
      {
        network <- network[sample(nrow(network), 40),]
        geneset = unique(c(network$GENE1,network$GENE2))
      }
      id <- geneset
      name <- id
      nodeData <- data.frame(id, name, stringsAsFactors = FALSE)
      source <- network$GENE1
      target <- network$GENE2
      edgeData <- data.frame(source, target, stringsAsFactors = FALSE)
      cyNetwork <-
        createCytoscapeJsNetwork(nodeData, edgeData, edgeTargetShape = "none")
      rcytoscapejs(
        nodeEntries = cyNetwork$nodes,
        edgeEntries = cyNetwork$edges,
        highlightConnectedNodes = FALSE
      ) #resulting cytoscape
    }
    
  })
  
  output$cluster_sel <- #Code for the Cluster Text
    renderText({
      #code used to find cluster of interest at bottom
      GENE <- toupper(input$GENE)
      clust <- find_cluster(GENE)
      validate(need((nrow(clust) > 0), ""))
      clust = unique(toupper(clust$Cluster))
      string = paste("\nGene is found in Cluster ", clust, "\n")
      return(string)
    })
  
  onclick("cytoplot",
          updateTextInput(session, "GENE", value = input$clickedNode)) #when clicking a node the cytoscape reappears
  
  # output$compar_metric <- renderText({
  #   return("go")
  # })
  
  output$stat_test_text <- renderText({ # Calculates the statistic below the plot (code found in functions)
    stat_calc(input)
  })
  
  wilcoxon_calc <- eventReactive(input$run_wilcoxon_Tissue, { #code used to run wilcoxon rank sum test
    GENE <- toupper(input$GENE)
    GENE2 <- toupper(get_gene2(input))
    if (GENE2 == "")
    {
      GENE2 = GENE
    }
    data <- getdata(GENE, GENE2, input$libchoose, input)
    validate(need((nrow(data) > 3), paste("")))
    data <- droplevels(data)
    data <- data %>% select(Type, BF)
    vars <- unique(data$Type)
    results = c()
    for (var in vars) {
      X <- data %>% filter(Type == var) %>% pull(BF)
      Y <- data %>% filter(Type != var) %>% pull(BF)
      if (length(X) > 1) { # Specific code to find wilcox test p-val
        tissue_spec <-
          wilcox.test(X, Y, alternative = "two.sided", paired = FALSE)
        tissue_pval <- log(tissue_spec$p.value)
        # tissue_pval <- formatC(tissue_pval, format = "e", digits = 3)
        tissue_pval <- formatC(tissue_pval, digits = 3)
        results = c(results, tissue_pval)
      }
      else{
        #text = "Insufficient n"
        pval <- 1
        results <- c(results, pval)
      }
    }
    #This code can be easily converted to table output if needed
    #wilcoxontab = data.frame(Tiss.Type = vars,Log.P.Val = results)
    names(results) = vars
    minpval = min(results)
    var = names(results[results == minpval])
    wilcoxon_text = paste(
      "Min Log(P-Value) = ",
      as.character(minpval),
      "\n Null Hypothesis: ",
      var,
      " has the same BF distribution as the rest of the selected tissue types.",
      sep = ""
    )
    return(wilcoxon_text)
  })
  
  output$wilcoxon_tiss <- renderText({
    wilcoxon_calc()
  })
  
  # output$wilcoxon_tiss <- DT::renderDataTable(DT::datatable({
  #   data <- wilcoxon_calc()
  # }))
  
  output$downloadData <- downloadHandler(
    #downloads data text files
    filename = function() {
      paste(input$GENE, ".csv", sep = "")
    },
    content = function(file) {
      GENE <- toupper(input$GENE)
      GENE2 <- toupper(get_gene2(input))
      if (GENE2 == "")
      {
        GENE2 <- GENE
      }
      data <- getdata(GENE, GENE2, input$libchoose, input)
      if (GENE == GENE2) {
        data <- data %>% select(-GENE2, -BF2)
      }
      if (GENE != GENE2) {
        plot_checker <- plot_select(input)
        if (5 %in% plot_checker)
        {
          data <- data %>% select(CCLE_CELL, Type, GENE, BF, GENE2, BF2)
        }
        if (6 %in% plot_checker)
        {
          data$EXP = log(data$EXP + 0.001)
          data$EXP[data$EXP < -2] = -2
          data$EXP = round(data$EXP, digits = 3)
          data <- data %>% select(CCLE_CELL, Type, GENE, BF, GENE2, EXP)
        }
        if (7 %in% plot_checker)
        {
          data <- data %>% select(CCLE_CELL, Type, GENE, BF, GENE2, MUT)
        }
        if (8 %in% plot_checker)
        {
          data$CN[data$CN < -6] = -6
          data <- data %>% select(CCLE_CELL, Type, GENE, BF, GENE2, CN)
        }
        
      }
      write.csv(data, file, row.names = FALSE)
    }
  )
  
  
  # observe({ #for using different tests at one point...
  #   plot_selected = plot_select(input)
  #   plot_selected = as.numeric(plot_selected)
  #   plot_selected = max(plot_selected)
  #   if (plot_selected %in% c(5, 6, 8))
  #   {
  #     updateSelectInput(session, "stat_test", choices = c("Pearson", "Spearman"))
  #   }
  #   if (plot_selected %in% c(7))
  #   {
  #     updateSelectInput(session, "stat_test", choices = c("Wilcoxon Rank Sum Test"))
  #   }
  # })
  # 
  
  # No Longer reactive code - was used for updating selections of wilcox rank sum test options
  # observe({
  #   tiss_selected = check_groups(input)
  #   if(input$libchoose == 1){
  #     tissues = c("Bladder","Bone","Brain","Breast",
  #                 "Colon/Colorectal","Endometrial/Uterine","Esophageal",
  #                 "Fibroblast","Gastric","Head and Neck",
  #                 "Kidney","Leukemia","Liver","Lung",
  #                 "Lymphoma","Neuroblastoma","Ovary",
  #                 "Pancreas","Sarcoma","Skin","Thyroid")
  #     tissues = tissues[as.numeric(tiss_selected)]
  #   }
  #   if(input$libchoose == 2){
  #     tissues = c("Bone","Breast","Colon/Colorectal","Leukemia",
  #                 "Lung","Ovary","Pancreas","Prostate",
  #                 "Rhabdoid","Skin")
  #
  #     tissues = tissues[as.numeric(tiss_selected)]
  #   }
  #   if(input$libchoose == 3){
  #     tissues = c("Breast","Cervical","Colon/Colorectal",
  #                 "Immortalized","Leukemia","Pancreas")
  #     tissues = tissues[as.numeric(tiss_selected)]
  #   }
  #
  #   if(input$libchoose == 4){
  #     tissues = c("Leukemia","Lymphoma")
  #     tissues = tissues[as.numeric(tiss_selected)]
  #   }
  #   if(input$libchoose == 5){
  #     tissues = c("Leukemia")
  #     tissues = tissues[as.numeric(tiss_selected)]
  #   }
  #   if(input$libchoose == 6){
  #     tissues = c("Breast","Colon/Colorectal","Ovary","Pancreas")
  #     tissues = tissues[as.numeric(tiss_selected)]
  #   }
  #   if(input$libchoose == 7){
  #     tissues = c("Bile Duct", "Bladder", "Bone", "Brain", "Breast",
  #                 "Cervical", "Colon/Colorectal", "Endometrial/Uterine", "Esophageal",
  #                 "Fibroblast", "Gastric", "Head and Neck", "Kidney",
  #                 "Leukemia", "Liver", "Lung", "Lymphoma",
  #                 "Myeloma", "Neuroblastoma", "Ovary", "Pancreas",
  #                 "Rhabdoid", "Sarcoma", "Skin", "Thyroid")
  #     tissues = tissues[as.numeric(tiss_selected)]
  #   }
  #   if(input$libchoose == 8){
  #     tissues = c( "Bile Duct","Bone","Brain","Breast",
  #                  "Colon/Colorectal","Endometrium/Uterine","Esophagus",
  #                  "Eye","Gastric","Head and Neck","Kidney",
  #                  "Leukemia","Lung","Lymphoma", "Myeloma",
  #                  "Neuroblastoma","Non-Cancerous", "Ovary",
  #                  "Pancreas","Prostate","Sarcoma","Skin")
  #     tissues = tissues[as.numeric(tiss_selected)]
  #   }
  #   updateSelectInput(session, "wilcoxon_tissue",
  #                     choices = tissues
  #   )})
  
  observe({
    #reset button - code updates interface when clicked
    if (input$libchoose == 1) {
      if (input$reset %% 2 == 0)
      {
        updateCheckboxGroupInput(
          session,
          "checkGroup1",
          "Select Tissues",
          choices = list(
            "Bladder" = 1,
            "Bone" = 2,
            "Brain" = 3,
            "Breast" = 4,
            "Colon/Colorectal" = 5,
            "Endometrial/Uterine" = 6,
            "Esophagael" = 7,
            "Fibroblast" = 8,
            "Gastric" = 9,
            "Head and Neck" = 10,
            "Kidney" = 11,
            "Leukemia" = 12,
            "Liver" = 13,
            "Lung" = 14,
            "Lymphoma" = 15,
            "Neuroblastoma" = 16,
            "Ovary" = 17,
            "Pancreas" = 18,
            "Sarcoma" = 19,
            "Skin" = 20,
            "Thyroid" = 21
          ),
          selected = c(
            1,
            2,
            3,
            4,
            5,
            6,
            7,
            8,
            9,
            10,
            11,
            12,
            13,
            14,
            15,
            16,
            17,
            18,
            19,
            20,
            21
          )
        )
      }
      else
      {
        updateCheckboxGroupInput(
          session,
          "checkGroup1",
          "Select Tissues",
          choices = list(
            "Bladder" = 1,
            "Bone" = 2,
            "Brain" = 3,
            "Breast" = 4,
            "Colon/Colorectal" = 5,
            "Endometrial/Uterine" = 6,
            "Esophagael" = 7,
            "Fibroblast" = 8,
            "Gastric" = 9,
            "Head and Neck" = 10,
            "Kidney" = 11,
            "Leukemia" = 12,
            "Liver" = 13,
            "Lung" = 14,
            "Lymphoma" = 15,
            "Neuroblastoma" = 16,
            "Ovary" = 17,
            "Pancreas" = 18,
            "Sarcoma" = 19,
            "Skin" = 20,
            "Thyroid" = 21
          ),
          selected = c()
        )
      }
    }
    if (input$libchoose == 2) {
      if (input$reset2 %% 2 == 0)
      {
        updateCheckboxGroupInput(
          session,
          "checkGroup2",
          "Select Tissues",
          choices = list(
            "Bone" = 1,
            "Breast" = 2,
            "Colon/Colorectal" = 3,
            "Leukemia" = 4,
            "Lung" = 5,
            "Ovary" = 6,
            "Pancreas" = 7,
            "Prostate" = 8,
            "Rhabdoid" = 9,
            "Skin" = 10
          ),
          selected = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
        )
        
      }
      else
      {
        updateCheckboxGroupInput(
          session,
          "checkGroup2",
          "Select Tissues",
          choices = list(
            "Bone" = 1,
            "Breast" = 2,
            "Colon/Colorectal" = 3,
            "Leukemia" = 4,
            "Lung" = 5,
            "Ovary" = 6,
            "Pancreas" = 7,
            "Prostate" = 8,
            "Rhabdoid" = 9,
            "Skin" = 10
          ),
          selected = c()
        )
        
      }
    }
    if (input$libchoose == 3) {
      if (input$reset3 %% 2 == 0)
      {
        updateCheckboxGroupInput(
          session,
          "checkGroup3",
          "Select Tissues",
          choices = list(
            "Breast" = 1,
            "Cervical" = 2,
            "Colon/Colorectal" = 3,
            "Immortalized" = 4,
            "Leukemia" = 5,
            "Pancreas" = 6
          ),
          selected = c(1, 2, 3, 4, 5, 6)
        )
        
      }
      else
      {
        updateCheckboxGroupInput(
          session,
          "checkGroup3",
          "Select Tissues",
          choices = list(
            "Breast" = 1,
            "Cervical" = 2,
            "Colon/Colorectal" = 3,
            "Immortalized" = 4,
            "Leukemia" = 5,
            "Pancreas" = 6
          ),
          selected = c()
        )
      }
    }
    if (input$libchoose == 4) {
      if (input$reset4 %% 2 == 0)
      {
        updateCheckboxGroupInput(
          session,
          "checkGroup4",
          "Select Tissues",
          choices = list("Leukemia" = 1, "Lymphoma" = 2),
          selected = c(1, 2)
        )
      }
      else
      {
        updateCheckboxGroupInput(
          session,
          "checkGroup4",
          "Select Tissues",
          choices = list("Leukemia" = 1, "Lymphoma" = 2),
          selected = c()
        )
      }
    }
    if (input$libchoose == 5) {
      if (input$reset5 %% 2 == 0)
      {
        updateCheckboxGroupInput(
          session,
          "checkGroup5",
          "Select Tissues",
          choices = list("Leukemia" = 1),
          selected = c(1)
        )
      }
      else
      {
        updateCheckboxGroupInput(
          session,
          "checkGroup5",
          "Select Tissues",
          choices = list("Leukemia"),
          selected = c()
        )
      }
    }
    if (input$libchoose == 6) {
      if (input$reset6 %% 2 == 0)
      {
        updateCheckboxGroupInput(
          session,
          "checkGroup6",
          "Select Tissues",
          choices = list(
            "Breast" = 1,
            "Colon/Colorectal" = 2,
            "Ovary" = 3,
            "Pancreas" = 4
          ),
          selected = c(1, 2, 3, 4)
        )
      }
      else
      {
        updateCheckboxGroupInput(
          session,
          "checkGroup6",
          "Select Tissues",
          choices = list(
            "Breast" = 1,
            "Colon/Colorectal" = 2,
            "Ovary" = 3,
            "Pancreas" = 4
          ),
          selected = c()
        )
      }
    }
    if (input$libchoose == 7) {
      if (input$reset7 %% 2 == 0)
      {
        updateCheckboxGroupInput(
          session,
          "checkGroup7",
          "Select Tissues",
          choices = list(
            "Bile Duct" = 1,
            "Bladder" = 2,
            "Bone" = 3,
            "Brain" = 4,
            "Breast" = 5,
            "Cervical" = 6,
            "Colon/Colorectal" = 7,
            "Endometrial/Uterine" = 8,
            "Esophagael" = 9,
            "Fibroblast" = 10,
            "Gastric" = 11,
            "Head and Neck" = 12,
            "Kidney" = 13,
            "Leukemia" = 14,
            "Liver" = 15,
            "Lung" = 16,
            "Lymphoma" = 17,
            "Myeloma" = 18,
            "Neuroblastoma" = 19,
            "Ovary" = 20,
            "Pancreas" = 21,
            "Rhabdoid" = 22,
            "Sarcoma" = 23,
            "Skin" = 24,
            "Thyroid" = 25
          ),
          selected = c(
            1,
            2,
            3,
            4,
            5,
            6,
            7,
            8,
            9,
            10,
            11,
            12,
            13,
            14,
            15,
            16,
            17,
            18,
            19,
            20,
            21,
            22,
            23,
            24,
            25
          )
        )
        
      }
      else
      {
        updateCheckboxGroupInput(
          session,
          "checkGroup7",
          "Select Tissues",
          choices = list(
            "Bile Duct" = 1,
            "Bladder" = 2,
            "Bone" = 3,
            "Brain" = 4,
            "Breast" = 5,
            "Cervical" = 6,
            "Colon/Colorectal" = 7,
            "Endometrial/Uterine" = 8,
            "Esophagael" = 9,
            "Fibroblast" = 10,
            "Gastric" = 11,
            "Head and Neck" = 12,
            "Kidney" = 13,
            "Leukemia" = 14,
            "Liver" = 15,
            "Lung" = 16,
            "Lymphoma" = 17,
            "Myeloma" = 18,
            "Neuroblastoma" = 19,
            "Ovary" = 20,
            "Pancreas" = 21,
            "Rhabdoid" = 22,
            "Sarcoma" = 23,
            "Skin" = 24,
            "Thyroid" = 25
          ),
          selected = c()
        )
        
      }
    }
    if (input$libchoose == 8) {
      if (input$reset8 %% 2 == 0)
      {
        updateCheckboxGroupInput(
          session,
          "checkGroup8",
          "Select Tissues",
          choices = list(
            "Bile Duct" = 1,
            "Bone" = 2,
            "Brain" = 3,
            "Breast" = 4,
            "Colon/Colorectal" = 5,
            "Endometrium/Uterine" = 6,
            "Esophagus" = 7,
            "Eye" = 8,
            "Gastric" = 9,
            "Head and Neck" = 10,
            "Kidney" = 11,
            "Leukemia" = 12,
            "Lung" = 13,
            "Lymphoma" = 14,
            "Myeloma" = 15,
            "Neuroblastoma" = 16,
            "Non-Cancerous" = 17,
            "Ovary" = 18,
            "Pancreas" = 19,
            "Prostate" = 20,
            "Sarcoma" = 21,
            "Skin" = 22
          ),
          selected = c(
            1,
            2,
            3,
            4,
            5,
            6,
            7,
            8,
            9,
            10,
            11,
            12,
            13,
            14,
            15,
            16,
            17,
            18,
            19,
            20,
            21,
            22
          )
        )
        
      }
      else
      {
        updateCheckboxGroupInput(
          session,
          "checkGroup8",
          "Select Tissues",
          choices = list(
            "Bile Duct" = 1,
            "Bone" = 2,
            "Brain" = 3,
            "Breast" = 4,
            "Colon/Colorectal" = 5,
            "Endometrium/Uterine" = 6,
            "Esophagus" = 7,
            "Eye" = 8,
            "Gastric" = 9,
            "Head and Neck" = 10,
            "Kidney" = 11,
            "Leukemia" = 12,
            "Lung" = 13,
            "Lymphoma" = 14,
            "Myeloma" = 15,
            "Neuroblastoma" = 16,
            "Non-Cancerous" = 17,
            "Ovary" = 18,
            "Pancreas" = 19,
            "Prostate" = 20,
            "Sarcoma" = 21,
            "Skin" = 22
          ),
          selected = c()
        )
      }
    }
    #session$onSessionEnded(stopApp)
  })
  
  #Rest of the code is for downloading text files
  
  output$downloadNet_Ext <-
    downloadHandler(
      #Text download functions
      filename = function() {
        paste("CoessNetwork_CRISPR_v1_Extended.txt", sep = "")
      },
      content = function(file) {
        file.copy("data/CoessNetwork_CRISPR_v1_Extended.txt", file)
      }
    )
  output$downloadBehan_19 <-
    downloadHandler(
      #Text download functions
      filename = function() {
        paste("qbf_Behan_2019.txt", sep = "")
      },
      content = function(file) {
        file.copy("data/qbf_Behan_2019.txt", file)
      }
    )
  
  output$downloadavana_18 <-
    downloadHandler(
      #Text download functions
      filename = function() {
        paste("qbf_Avanadata_2018.txt", sep = "")
      },
      content = function(file) {
        file.copy("data/qbf_Avanadata_2018.txt", file)
      }
    )
  output$downloadavana_17 <-
    downloadHandler(
      #Text download functions
      filename = function() {
        paste("qbf_Avanadata_2017.txt", sep = "")
      },
      content = function(file) {
        file.copy("data/qbf_Avanadata_2017.txt", file)
      }
    )
  
  output$downloadgecko_data <-
    downloadHandler(
      #Text download functions
      filename = function() {
        paste("qbf_gecko.txt", sep = "")
      },
      content = function(file) {
        file.copy("data/qbf_gecko.txt", file)
      }
    )
  
  output$downloadwang_data <-
    downloadHandler(
      #Text download functions
      filename = function() {
        paste("qbf_wang.txt", sep = "")
      },
      content = function(file) {
        file.copy("data/qbf_wang.txt", file)
      }
    )
  
  output$downloadtzelepis_data <-
    downloadHandler(
      #Text download functions
      filename = function() {
        paste("qbf_tzelepis_yusa.txt", sep = "")
      },
      content = function(file) {
        file.copy("data/qbf_tzelepis_yusa.txt", file)
      }
    )
  
  output$downloadshRNA_data <-
    downloadHandler(
      #Text download functions
      filename = function() {
        paste("qbf_shrna.txt", sep = "")
      },
      content = function(file) {
        file.copy("data/qbf_shrna.txt", file)
      }
    )
  
  
  output$downloadtkov1_data <-
    downloadHandler(
      #Text download functions
      filename = function() {
        paste("qbf_tkov1.txt", sep = "")
      },
      content = function(file) {
        file.copy("data/qbf_tkov1.txt", file)
      }
    )
  
  output$downloadCitation_avana18 <-
    downloadHandler(
      #Text download functions
      filename = function() {
        paste("citations_Avana_q3_2017.txt", sep = "")
      },
      content = function(file) {
        file.copy("./citations/citations_Avana_q3_2017.txt", file)
      }
    )
  
  output$downloadCitation_Behan19 <-
    downloadHandler(
      #Text download functions
      filename = function() {
        paste("citations_Behan19.txt", sep = "")
      },
      content = function(file) {
        file.copy("./citations/citations_Behan19.txt", file)
      }
    )
  
  output$downloadCitation_avana <-
    downloadHandler(
      #Text download functions
      filename = function() {
        paste("citations_Avana_q3_2017.txt", sep = "")
      },
      content = function(file) {
        file.copy("./citations/citations_Avana_q3_2017.txt", file)
      }
    )
  
  output$downloadCitation_shRNA <-
    downloadHandler(
      #Text download functions
      filename = function() {
        paste("citations_shRNA.txt", sep = "")
      },
      content = function(file) {
        file.copy("citations/citations_shRNA.txt", file)
      }
    )
  
  output$downloadCitation_GeCKO <-
    downloadHandler(
      #Text download functions
      filename = function() {
        paste("citations_GeCKO.txt", sep = "")
      },
      content = function(file) {
        file.copy("citations/citations_GeCKO.txt", file)
      }
    )
  
  output$downloadCitation_TKOv1 <-
    downloadHandler(
      #Text download functions
      filename = function() {
        paste("citations_TKOv1.txt", sep = "")
      },
      content = function(file) {
        file.copy("citations/citations_TKOv1.txt", file)
      }
    )
  
  output$downloadCitation_Wang <-
    downloadHandler(
      #Text download functions
      filename = function() {
        paste("citations_Wang.txt", sep = "")
      },
      content = function(file) {
        file.copy("citations/citations_Wang.txt", file)
      }
    )
  
  output$downloadCitation_Tzelepis <-
    downloadHandler(
      #Text download functions
      filename = function() {
        paste("citations_Tzelepis.txt", sep = "")
      },
      content = function(file) {
        file.copy("citations/citations_Tzelepis.txt", file)
      }
    )
  
  output$downloadNet <- downloadHandler(
    filename = function() {
      paste("CoessNetwork_CRISPR_v1.txt", sep = "")
    },
    content = function(file) {
      file.copy("data/CoessNetwork_CRISPR_v1.txt", file)
    }
  )
  
  output$downloadClusters <- downloadHandler(
    filename = function() {
      paste("CoessNetwork_CRISPR_v1_Cluster.txt", sep = "")
    },
    content = function(file) {
      file.copy("data/CoessNetwork_CRISPR_v1_Cluster.txt", file)
    }
  )
  
  # output$downloadPlot <- downloadHandler( #BROKEN CODE
  #   filename = paste(input$GENE,"_BFplot.pdf",sep = ""),
  #   content = function(file) {
  #     EXPort(ploter(),file)
  #     browseURL(file)
  #   })
  
}
