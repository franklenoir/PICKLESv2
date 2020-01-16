library(shiny)
library(shinyjs)
library(rcytoscapejs)
library(shinyWidgets)
library(shinythemes)
library(shinyBS)
library(plotly)
library(searchable)
library(shinydashboard)
library(GAlogger)

#Google Analytics 
ga_set_tracking_id("UA-124182410-1")
ga_set_approval(consent = TRUE)

#Load in page
appCSS <- "
#loading-content {
position: absolute;
background: #000000;
opacity: 0.9;
z-index: 100;
left: 0;
right: 0;
height: 100%;
text-align: center;
color: #FFFFFF;
}
"
#Main UI code
fluidPage(
  useShinyjs(),
  inlineCSS(appCSS),
  
  # Loading message
  div(id = "loading-content",
      h2("Loading...")),
  
  theme = shinytheme("flatly"),
  tags$head(includeScript("google-analytics.js")),
  #tags$a(img(src='logo.png', align = "right",height='150px',width='150px'),href="http://hart-lab.org/"),
  tags$div(
    class = "h2",
    checked = NA,
    "Welcome to PICKLES -",
    a("Citation", href = "https://academic.oup.com/nar/article/46/D1/D776/4564803")
  ),
  tags$div(
    class = "h5",
    checked = NA,
    "Developed by members of the",
    a("Hart Lab.", href = "http://hart-lab.org/"),
    "\tPlease report bugs to 'pickles at hart-lab.org'"
  ),
  # headerPanel("\t Welcome to PICKLES, the database of Pooled In vitro CRISPR Knockout Library Essentiality Screens -",
  #             a("Citation.", href = "https://academic.oup.com/nar/article/46/D1/D776/4564803"), "Developed by members of the",
  #             a("Hart Lab.", href = "http://hart-lab.org/"),"\tPlease report bugs to 'pickles at hart-lab.org'"),
  
  #Side bar - for all the gene search and library selection
  sidebarPanel(
    textInput("GENE", "Search a Gene", value = "ERBB2"),
    selectInput(
      "libchoose",
      label = h5("Select Screen Library"),
      choices = list(
        "Behan 2019" = 8,
        "Avana 2018 q4" = 7,
        "Avana 2017 q3 - Coessential Network" = 1,
        "GeCKO" = 2,
        "TKOv1" = 3,
        "Wang AML" = 4,
        "Tzelepis" = 5,
        "shRNA" = 6
      ),
      selected = 7
    ),
    conditionalPanel(
      condition = "input.libchoose == 1",
      textInput("GENE2_1", "Compare Selected Gene to Another Gene:", value = ""),
      selectInput(
        "Plot_selection1",
        label = "Select 2nd Gene Feature to be Plotted on Figure",
        choices = list(
          "Gene Normalized BF" = 5,
          "Gene Expression Log(RPKM)" = 6,
          "Gene Mutations" = 7
        ),
        selected = c(6)
      ),
      actionButton("reset", "Select All Tissues/Reset"),
      checkboxGroupInput(
        "checkGroup1",
        h5("Select Tissues"),
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
        selected = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21)
      ),
      downloadButton("downloadavana_17", "Download All BF Data", style =
                       'padding:4px; font-size:80%'),
      br(),
      downloadButton("downloadCitation_avana", "Citation Information", style =
                       'padding:4px; font-size:80%'),
      br(),
      downloadButton("downloadNet", "Download Network Core" , style = 
                       'padding:4px; font-size:80%'),
      br(),
      downloadButton("downloadNet_Ext", "Download Network Extended", style =
                       'padding:4px; font-size:80%'),
      br(),
      downloadButton("downloadClusters", "Download Gene Clusters", style =
                       'padding:4px; font-size:80%')
    ),
    conditionalPanel(
      textInput("GENE2_2", "Compare Searched Gene to:", value = ""),
      selectInput(
        "Plot_selection2",
        label = "Select 2nd Gene Feature to be Plotted on Figure",
        choices = list(
          "Gene Normalized BF" = 5,
          "Gene Expression Log(RPKM)" = 6
        ),
        selected = c(6)
      ),
      actionButton("reset2", "Select All/Reset"), #Plot Reset Code
      condition = "input.libchoose == 2",
      checkboxGroupInput(
        "checkGroup2",
        h5("Select Tissues"),
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
      ),
      downloadButton("downloadgecko_data", "Download All BF Data", style =
                       'padding:4px; font-size:80%'),
      br(),
      downloadButton("downloadCitation_GeCKO", "Citation Information", style =
                       'padding:4px; font-size:80%')
    ),
    conditionalPanel(
      textInput("GENE2_3", "Compare Searched Gene to:", value = ""),
      selectInput(
        "Plot_selection3",
        label = "Select 2nd Gene Feature to be Plotted on Figure",
        choices = list("Gene Normalized BF" = 5),
        selected = c(5)
      ),
      actionButton("reset3", "Select All/Reset"),
      condition = "input.libchoose == 3",
      checkboxGroupInput(
        "checkGroup3",
        h5("Select Tissues"),
        choices = list(
          "Breast" = 1,
          "Cervical" = 2,
          "Colon/Colorectal" = 3,
          "Immortalized" = 4,
          "Leukemia" = 5,
          "Pancreas" = 6
        ),
        selected = c(1, 2, 3, 4, 5, 6)
      ),
      downloadButton("downloadtkov1_data", "Download All BF Data", style =
                       'padding:4px; font-size:80%'),
      br(),
      downloadButton("downloadCitation_TKOv1", "Citation Information", style =
                       'padding:4px; font-size:80%')
    ),
    conditionalPanel(
      textInput("GENE2_4", "Compare Searched Gene to:", value = ""),
      selectInput(
        "Plot_selection4",
        label = "Select 2nd Gene Feature to be Plotted on Figure",
        choices = list("Gene Normalized BF" = 5),
        selected = c(5)
      ),
      actionButton("reset4", "Select All/Reset"),
      condition = "input.libchoose == 4",
      checkboxGroupInput(
        "checkGroup4",
        h5("Select Tissues"),
        choices = list("Leukemia" = 1, "Lymphoma" = 2),
        selected = c(1, 2)
      ),
      downloadButton("downloadwang_data", "Download All BF Data", style =
                       'padding:4px; font-size:80%'),
      br(),
      downloadButton("downloadCitation_Wang", "Citation Information", style =
                       'padding:4px; font-size:80%')
    ),
    conditionalPanel(
      textInput("GENE2_5", "Compare Searched Gene to:", value = ""),
      selectInput(
        "Plot_selection5",
        label = "Select 2nd Gene Feature to be Plotted on Figure",
        choices = list("Gene Normalized BF" = 5),
        selected = c(5)
      ),
      actionButton("reset5", "Select All/Reset"),
      condition = "input.libchoose == 5",
      checkboxGroupInput(
        "checkGroup5",
        h5("Select Tissues"),
        choices = list("Leukemia" = 1),
        selected = c(1)
      ),
      downloadButton("downloadtzelepis_data", "Download All BF Data", style =
                       'padding:4px; font-size:80%'),
      br(),
      downloadButton("downloadCitation_Tzelepis", "Citation Information" , style =
                       'padding:4px; font-size:80%')
    ),
    conditionalPanel(
      textInput("GENE2_6", "Compare Searched Gene to:", value = ""),
      selectInput(
        "Plot_selection6",
        label = "Select 2nd Gene Feature to be Plotted on Figure",
        choices = list("Gene Normalized BF" = 5, "Log (RPKM)" = 6),
        selected = c(6)
      ),
      actionButton("reset6", "Select All/Reset"),
      condition = "input.libchoose == 6",
      checkboxGroupInput(
        "checkGroup6",
        h5("Select Tissues"),
        choices = list(
          "Breast" = 1,
          "Colon/Colorectal" = 2,
          "Ovary" = 3,
          "Pancreas" = 4
        ),
        selected = c(1, 2, 3, 4)
      ),
      downloadButton("downloadshRNA_data", "Download All BF Data", style =
                       'padding:4px; font-size:80%'),
      br(),
      downloadButton("downloadCitation_shRNA", "Citation Information", style =
                       'padding:4px; font-size:80%')
    ),
    conditionalPanel(
      condition = "input.libchoose == 7",
      textInput("GENE2_7", "Compare Selected Gene to Another Gene:", value = ""),
      selectInput(
        "Plot_selection7",
        label = "Select 2nd Gene Feature to be Plotted on Figure",
        choices = list(
          "Gene Normalized BF" = 5,
          "Gene Expression Log (TPM)" = 6,
          "Gene Mutations" = 7,
          "Gene Copy Number" = 8
        ),
        selected = c(6)
      ),
      actionButton("reset7", "Select All Tissues/Reset"),
      checkboxGroupInput(
        "checkGroup7",
        h5("Select Tissues"),
        choices = list(
          "Bile Duct" = 1,
          "Bladder" = 2,
          "Bone" = 3,
          "Brain" = 4,
          "Breast" = 5,
          "Cervical" = 6,
          "Colon/Colorectal" = 7,
          "Endometrial/Uterine" = 8,
          "Esophageal" = 9,
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
      ),
      downloadButton("downloadavana_18", "Download All BF Data", style =
                       'padding:4px; font-size:80%'),
      br(),
      downloadButton("downloadCitation_avana18", "Citation Information", style =
                       'padding:4px; font-size:80%'),
      br()
    ),
    conditionalPanel(
      condition = "input.libchoose == 8",
      textInput("GENE2_8", "Compare Selected Gene to Another Gene:", value = ""),
      selectInput(
        "Plot_selection8",
        label = "Select 2nd Gene Feature to be Plotted on Figure",
        choices = list("Gene Normalized BF" = 5),
        selected = c(5)
      ),
      actionButton("reset8", "Select All Tissues/Reset"),
      checkboxGroupInput(
        "checkGroup8",
        h5("Select Tissues"),
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
        selected = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22)
      ),
      downloadButton("downloadBehan_19", "Download All BF Data", style =
                       'padding:4px; font-size:80%'),
      br(),
      downloadButton("downloadCitation_Behan19", "Citation Information", style =
                       'padding:4px; font-size:80%'),
      br()
    )
  ),
  shiny::hr(),
  #Main Panel - Plots, Correlations, and network
  mainPanel(
    plotlyOutput("plot", width = '100%'),
    textOutput("stat_test_text", container = span),
    br(),
    br(),
    downloadButton("downloadData", "Download Displayed Data", style =
                     'padding:4px; font-size:80%'),
    br(),
    br(),
    splitLayout(
      cellWidths = 350,
      cellArgs = list(style = "padding: 10px"),
      #selectInput('wilcoxon_tissue', "BF Tissue Specificity Wilcoxon Rank Sum Test", choices = ""),
      #Ommitted for now - code currently slows down interface
      tags$b(
        "BF Tissue Specificity Wilcoxon Rank Sum Test",
        tags$br(),
        "Output = most signficant tissue"
      ),
      
      bsButton(
        "run_wilcoxon_Tissue",
        "\n Calculate",
        icon = icon("angle-double-right"),
        style = "info",
        size = "medium"
      ),
      tags$head(tags$style(
        HTML(".shiny-split-layout > div {overflow: visible;}")
      ))
    ),
    
    #DT::dataTableOutput("wilcoxon_tiss"), #for table results instead of single
    textOutput("wilcoxon_tiss"),
    br(),
    br(),
    verbatimTextOutput("event")
  ),
  shiny::hr(),
  useShinyjs(),
  conditionalPanel(
    condition = "input.libchoose == 1",
    br(),
    mainPanel(rcytoscapejsOutput("cytoplot", height = "600px")),
    br(),
    br(),
    h3(textOutput("cluster_sel", container = span), align = "center")
  )
)
