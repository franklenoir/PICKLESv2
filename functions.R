#Functions for Pickles R Shiny App

'%!in%' <-
  function(x, y)
    ! ('%in%'(x, y)) # reverse of 'in' function - not in

#getdata outs data together from the sql db. GENE = input gene, GENE2 = input GENE2 libchoose = library choosen, input = dynamic user input.
getdata <- function(GENE, GENE2, libchoose, input) {
  con <-
    dbConnect(RSQLite::SQLite(), dbname = "pickles-database.db") #Connect to the database
  
  if (libchoose == 8) {
    # Each option is pretty consistent. Library selects the data set being used.
    library <- "Sanger_Q1_2019_dat"
    tissue_query <-
      "SELECT CCLE_Name,Primary_Disease FROM Master_Tissue_Key WHERE Sanger_Q1_19 == 1 " #Finds which tissue types are present in the database
    query <-
      paste("SELECT CCLE_CELL,GENE,BF FROM ",
            library,
            " WHERE GENE ='",
            GENE,
            "'",
            sep = "") # Selecting GENE1
    query2 <-
      paste(
        "SELECT CCLE_CELL,GENE as GENE2,BF AS BF2 FROM ",
        library,
        " WHERE GENE ='",
        GENE2,
        "'",
        sep = ""
      ) # Selects all metrics for GENE2
    tissue_check <-
      input$checkGroup8 #Checking to see which tissues were selected.
  }
  if (libchoose == 7) {
    library <- "Avana_Q4_2018_dat"
    tissue_query <-
      "SELECT CCLE_Name,Primary_Disease FROM Master_Tissue_Key WHERE Avana_Q4_18 == 1"
    query <-
      paste("SELECT CCLE_CELL,GENE,BF FROM ",
            library,
            " WHERE GENE ='",
            GENE,
            "'",
            sep = "")
    query2 <-
      paste(
        "SELECT CCLE_CELL,GENE as GENE2,BF AS BF2,MUT,EXP,CN FROM ",
        library,
        " WHERE GENE ='",
        GENE2,
        "'",
        sep = ""
      )
    tissue_check <- input$checkGroup7
  }
  if (libchoose == 6) {
    library <- "HART_2014_shRNA_dat"
    tissue_query <-
      "SELECT CCLE_Name,Primary_Disease FROM Master_Tissue_Key WHERE Hart_SHRNA == 1"
    query <-
      paste("SELECT CCLE_CELL,GENE,BF FROM ",
            library,
            " WHERE GENE ='",
            GENE,
            "'",
            sep = "")
    query2 <-
      paste(
        "SELECT CCLE_CELL,GENE as GENE2,BF AS BF2,EXP FROM ",
        library,
        " WHERE GENE ='",
        GENE2,
        "'",
        sep = ""
      )
    tissue_check <- input$checkGroup6
  }
  if (libchoose == 5) {
    library <- "Yusa_Tzelepis_AML_2016_dat"
    tissue_query <-
      "SELECT CCLE_Name,Primary_Disease FROM Master_Tissue_Key WHERE Tzelepis_AML == 1 "
    query <-
      paste("SELECT CCLE_CELL,GENE,BF FROM ",
            library,
            " WHERE GENE ='",
            GENE,
            "'",
            sep = "")
    query2 <-
      paste(
        "SELECT CCLE_CELL,GENE as GENE2,BF AS BF2 FROM ",
        library,
        " WHERE GENE ='",
        GENE2,
        "'",
        sep = ""
      )
    tissue_check <- input$checkGroup5
  }
  if (libchoose == 4) {
    library <- "Wang_AML_2016_dat"
    tissue_query <-
      "SELECT CCLE_Name,Primary_Disease FROM Master_Tissue_Key WHERE Wang_AML == 1"
    query <-
      paste("SELECT CCLE_CELL,GENE,BF FROM ",
            library,
            " WHERE GENE ='",
            GENE,
            "'",
            sep = "")
    query2 <-
      paste(
        "SELECT CCLE_CELL,GENE as GENE2,BF AS BF2 FROM ",
        library,
        " WHERE GENE ='",
        GENE2,
        "'",
        sep = ""
      )
    tissue_check <- input$checkGroup4
  }
  if (libchoose == 3) {
    library <- "Hart_TKOv1_dat"
    tissue_query <-
      "SELECT CCLE_Name,Primary_Disease FROM Master_Tissue_Key WHERE Hart_TKOv1 == 1"
    query <-
      paste("SELECT CCLE_CELL,GENE,BF FROM ",
            library,
            " WHERE GENE ='",
            GENE,
            "'",
            sep = "")
    query2 <-
      paste(
        "SELECT CCLE_CELL,GENE as GENE2,BF AS BF2 FROM ",
        library,
        " WHERE GENE ='",
        GENE2,
        "'",
        sep = ""
      )
    tissue_check <- input$checkGroup3
  }
  if (libchoose == 2) {
    library <- "Achilles_GeCKO_2016_dat"
    tissue_query <-
      "SELECT CCLE_Name,Primary_Disease FROM Master_Tissue_Key WHERE Achilles_GECKO == 1"
    query <-
      paste("SELECT CCLE_CELL,GENE,BF FROM ",
            library,
            " WHERE GENE ='",
            GENE,
            "'",
            sep = "")
    query2 <-
      paste(
        "SELECT CCLE_CELL,GENE as GENE2,BF AS BF2,EXP FROM ",
        library,
        " WHERE GENE ='",
        GENE2,
        "'",
        sep = ""
      )
    tissue_check <- input$checkGroup2
  }
  if (libchoose == 1) {
    library <- "Avana_Q3_2017_dat" #avana library
    tissue_query <-
      "SELECT CCLE_Name,Primary_Disease FROM Master_Tissue_Key WHERE Avana_Q3_17 == 1"
    query <-
      paste("SELECT CCLE_CELL,GENE,BF FROM ",
            library,
            " WHERE GENE ='",
            GENE,
            "'",
            sep = "")
    query2 <-
      paste(
        "SELECT CCLE_CELL,GENE as GENE2,BF AS BF2,MUT,EXP,CN FROM ",
        library,
        " WHERE GENE ='",
        GENE2,
        "'",
        sep = ""
      )
    tissue_check <- input$checkGroup1
  }
  
  data <-
    dbGetQuery(con, statement = query) #Calls GENE1 data from the database.
  data_gene2 <-
    dbGetQuery(con, statement = query2) #Calls GENE2 data from the database.
  data <- merge(data, data_gene2) #Merging data together
  tiss_data <- dbGetQuery(con, statement = tissue_query)
  rm(data_gene2, library, tissue_query, query, query2) #Removing variables no longer used
  dbDisconnect(con)
  
  p <-
    brewer.pal(12, "Paired") #tissue colors, picking the palette in order to standardize colors across
  p2 <- brewer.pal(8, "Dark2")
  p3 <- brewer.pal(9, "Set1")
  p4 <- brewer.pal(8, "Set2")
  p5 <- brewer.pal(12, "Set3")
  
  pal <-  c(p, p2, p3, p4, p5)
  rm(p, p2, p3, p4, p5)
  #All possible tissues (not all are used)
  tissues <-
    c(
      "BLADDER CANCER",
      "BONE CANCER",
      "BRAIN CANCER",
      "BREAST CANCER",
      "COLON/COLORECTAL CANCER",
      "ENDOMETRIAL/UTERINE CANCER",
      "ESOPHAGEAL CANCER",
      "FIBROBLAST",
      "GASTRIC CANCER",
      "HEAD AND NECK CANCER",
      "KIDNEY CANCER",
      "LEUKEMIA",
      "LIVER CANCER",
      "LUNG CANCER",
      "LYMPHOMA",
      "NEUROBLASTOMA",
      "OVARIAN CANCER",
      "PANCREATIC CANCER",
      "SARCOMA",
      "SKIN CANCER",
      "THYROID CANCER",
      "PROSTATE CANCER",
      "RHABDOID",
      "CERVICAL CANCER",
      "IMMORT",
      "IMMORTALIZED",
      "BILE DUCT CANCER",
      "MYELOMA",
      "EYE CANCER",
      "LIPOSARCOMA",
      "NON-CANCEROUS",
      "ADRENAL CANCER",
      "EMBRYONAL CANCER",
      "GALLBLADDER CANCER",
      "PRIMARY CELLS"
    )
  
  
  colnames(tiss_data) = c("CCLE_CELL", "Type")
  colors <-
    pal[1:length(tissues)] #tissue colors, code is necessary to keep colors consistent
  names(colors) <- tissues
  tiss_data$colors <-
    colors[tiss_data$Type] #assigning the colors to the tissue type
  
  tiss_data <- tissue_filter(tiss_data, tissue_check)
  data <- merge(tiss_data, data)
  rm(tiss_data, pal, colors, tissues)
  
  return(data)
}

get_coess_network <-
  function(GENE, geneset) {
    #Get network data from the same SQL database
    con <-
      dbConnect(RSQLite::SQLite(), dbname = "pickles-database.db") #Connect to the database
    geneset <-
      paste(geneset, collapse = "\',\'") #Parsing geneset vector to be read by SQL
    geneset <- paste("\'", geneset, "\'", sep = "")
    query <-
      paste(
        "SELECT * FROM Coess_Network WHERE GENE1 in (",
        geneset,
        ") AND GENE2 in (",
        geneset,
        ")",
        sep = ""
      )
    data <- dbGetQuery(con, statement = query)
    dbDisconnect(con)
    return(data)
  }

find_cluster <-
  function(GENE, input) {
    #Function used for finding specific clusters of genes
    con <-
      dbConnect(RSQLite::SQLite(), dbname = "pickles-database.db") #Connect to the database
    query <-
      paste("SELECT * FROM Coess_Network_Clusters WHERE GENE ='",
            GENE,
            "'",
            sep = "")
    data <- dbGetQuery(con, statement = query)
    clust <- data$Cluster[1]
    query <-
      paste("SELECT * FROM Coess_Network_Clusters WHERE Cluster ='",
            clust,
            "'",
            sep = "")
    data <- dbGetQuery(con, statement = query)
    dbDisconnect(con)
    return(data)
  }

tissue_filter <-
  function (tissue, tissue_check) {
    #tissue filter function (selects tissues that are desired from buttons)
    types <- levels(factor(tissue$Type))
    types <- types[as.numeric(tissue_check)]
    tissue <- tissue %>% filter(Type %in% types)
    return(tissue)
  }

plot_select <-
  function(input) {
    #function used for bf, expression, CN, Mut buttons
    library <- input$libchoose
    if (library == 1) {
      plot_select <- input$Plot_selection1
      checkGroup <- input$checkGroup1
    }
    if (library == 2) {
      plot_select <- input$Plot_selection2
      checkGroup <- input$checkGroup2
    }
    if (library == 3) {
      plot_select <- input$Plot_selection3
      checkGroup <- input$checkGroup3
    }
    if (library == 4) {
      plot_select <- input$Plot_selection4
      checkGroup <- input$checkGroup4
    }
    if (library == 5) {
      plot_select <- input$Plot_selection5
      checkGroup <- input$checkGroup5
    }
    if (library == 6) {
      plot_select <- input$Plot_selection6
      checkGroup <- input$checkGroup6
    }
    if (library == 7) {
      plot_select <- input$Plot_selection7
      checkGroup <- input$checkGroup7
    }
    if (library == 8) {
      plot_select <- input$Plot_selection8
      checkGroup <- input$checkGroup8
    }
    return(plot_select)
  }

check_groups <-
  function(input) {
    #old code that was used for tissue selection, now just used to identify if a tissue is selected
    library <- input$libchoose
    if (library == 1) {
      checkGroup <- input$checkGroup1
    }
    if (library == 2) {
      checkGroup <- input$checkGroup2
    }
    if (library == 3) {
      checkGroup <- input$checkGroup3
    }
    if (library == 4) {
      checkGroup <- input$checkGroup4
    }
    if (library == 5) {
      checkGroup <- input$checkGroup5
    }
    if (library == 6) {
      checkGroup <- input$checkGroup6
    }
    if (library == 7) {
      checkGroup <- input$checkGroup7
    }
    if (library == 8) {
      checkGroup <- input$checkGroup8
    }
    return(checkGroup)
  }

get_gene2 <- function(input) {
  #identifying second gene input
  library <- input$libchoose
  if (library == 1) {
    GENE2 <- input$GENE2_1
  }
  if (library == 2) {
    GENE2 <- input$GENE2_2
  }
  if (library == 3) {
    GENE2 <- input$GENE2_3
  }
  if (library == 4) {
    GENE2 <- input$GENE2_4
  }
  if (library == 5) {
    GENE2 <- input$GENE2_5
  }
  if (library == 6) {
    GENE2 <- input$GENE2_6
  }
  if (library == 7) {
    GENE2 <- input$GENE2_7
  }
  if (library == 8) {
    GENE2 <- input$GENE2_8
  }
  return(GENE2)
}

stat_calc <-
  function(input) {
    #Code used to identify and generate the statistic below the plot code
    text <- ""
    GENE <- toupper(input$GENE)
    GENE2 <- get_gene2(input)
    GENE2 <- toupper(GENE2)
    if (GENE2 == "")
    {
      GENE2 <- GENE
    }
    plot_selected  <-
      as.numeric(plot_select(input)) #checking which plot is being used
    stat_test <-
      as.character(toupper(input$stat_test)) # at one point there was going to be multiple options. Not the case right now
    data <-
      getdata(GENE, GENE2, input$libchoose, input) #getting code. Probably can reduce code by only calling this once
    validate(need((nrow(data) > 1), paste(""))) #make sure sufficient data is present
    data <-
      droplevels(data) #drop extra levels of data that isn't selected tissue.
    if (plot_selected == 5) {
      #BF
      data <- data %>% select(BF, BF2)
      data <- data %>% na.omit
      validate(need((nrow(data) > 2), paste("Insuffcient Data")))
      cor_res <-
        cor.test(data$BF, data$BF2, method = "pearson") #Correlation Test of BFs
      text <-
        paste(
          "Pearson Correlation of ",
          GENE,
          " BF vs ",
          GENE2,
          " BF = ",
          round(cor_res$estimate, 3),
          sep = ""
        )
    }
    if (plot_selected == 6) {
      #EXP
      data <- data %>% select(BF, EXP)
      data <- data %>% na.omit
      validate(need((nrow(data) > 2), paste("Insuffcient Exp Data")))
      cor_res <-
        cor.test(data$BF, data$EXP, method = "pearson") #Correlation Test of BF vs EXP
      text <-
        paste(
          "Pearson Correlation of ",
          GENE,
          " BF vs ",
          GENE2,
          " Expression = ",
          round(cor_res$estimate, 3),
          sep = ""
        )
    }
    if (plot_selected == 7) {
      #Mutations
      data <- data %>% select(BF, MUT)
      X <- data %>% filter(MUT == "Mut")
      Y <- data %>% filter(MUT != "Mut")
      validate(need((nrow(X) > 2), paste("Insuffcient Mut Data")))
      wilcox_res <-
        wilcox.test(Y$BF, X$BF, alternative = "two.sided", paired = FALSE) #Rank Sum Tests comparing mut vs nonmut stats
      pval <- wilcox_res$p.value
      pval <- formatC(pval, format = "e", digits = 3)
      text <-
        paste(
          "Wilcoxon Rank Sum Test P-Value = ",
          as.character(pval),
          "\n Null Hypothesis: ",
          GENE,
          " has the same distribtion of BFs when ",
          GENE2,
          " is Mutated",
          sep = ""
        )
      
    }
    if (plot_selected == 8) {
      #CN
      data <- data %>% select(BF, CN)
      data <- data %>% na.omit
      validate(need((nrow(data) > 2), paste("Insuffcient CN Data")))
      cor_res <-
        cor.test(data$BF, data$CN, method = "pearson") #Correlation Test of BF vs CN
      text <-
        paste(
          "Pearson Correlation of ",
          GENE,
          " BF vs ",
          GENE2,
          " Copy Number = ",
          round(cor_res$estimate, 3),
          sep = ""
        )
    }
    return(text) #text to return for test
  }
