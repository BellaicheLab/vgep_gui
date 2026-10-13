# --- Libraries ---
library(shiny)
library(shinyjs)
library(shinyWidgets)
library(shinythemes)
library(ggplot2)
library(ggiraph)
library(DT)

# --- Data Loading ---

# check.names = FALSE keeps gene symbols such as Su(fu) or 14-3-3epsilon as written
AtlasData <- read.csv("data/Atlas.csv", check.names = FALSE)
ValidationData <- read.csv("data/Validation.csv", check.names = FALSE)
vGEPData <- read.csv("data/vGEP.csv", check.names = FALSE)
macrochaetaePosition <- read.csv("data/macrochaete.csv", sep = ";")
mdDomains <- read.csv("data/md_domains.csv")

# --- Data Preprocessing ---

"%notin%" <- Negate("%in%")

# Extract spatial coordinates, then drop them from data frames
vGEPcoord <- data.frame(coordX = vGEPData$coordX, coordY = vGEPData$coordY)
Atlascoord <- data.frame(coordX = AtlasData$coordX, coordY = AtlasData$coordY)
ValidationCoord <- data.frame(coordX = ValidationData$coordX, coordY = ValidationData$coordY)
coord_cols <- c("coordX", "coordY")
vGEPData <- vGEPData[, !(colnames(vGEPData) %in% coord_cols)]
AtlasData <- AtlasData[, !(colnames(AtlasData) %in% coord_cols)]
ValidationData <- ValidationData[, !(colnames(ValidationData) %in% coord_cols)]

# Annotate validation genes with measurement type
protein_val <- c("C15", "caup", "Dl", "Wnt6", "svp", "vvl", "18w", "Tollo", "yrt", "stan")
fish_val <- c("slow", "sli", "pdm3", "toe", "crp", "nw", "hbs", "hth", "robo1", "Qsox2", "fj", "dpp", "Pvr")
val_names <- colnames(ValidationData)
val_names[val_names %in% protein_val] <- paste0(val_names[val_names %in% protein_val], " (protein)")
val_names[val_names %in% fish_val] <- paste0(val_names[val_names %in% fish_val], " (RNA)")
colnames(ValidationData) <- val_names

# Gene name lists for UI selectors
names_atlas <- sort(colnames(AtlasData))
names_validation <- sort(colnames(ValidationData))
names_vGEP <- sort(colnames(vGEPData))

# Build base plot for the Draw tab ROI selection
draw.data <- data.frame(expression = 0, vGEPcoord)
draw.data$ID <- paste0("X", draw.data$coordX, "Y", draw.data$coordY)
draw.selection <- ggplot(draw.data, aes(x = coordY, y = coordX, color = expression)) +
  geom_point_interactive(aes(data_id = ID), size = 3) +
  scale_colour_gradient(low = "white", high = "azure4") +
  coord_fixed(ratio = 1) +
  scale_y_reverse() +
  theme_void() +
  theme(legend.position = "none")

# Morphodynamic domains (5-MD1..5), coloured as in the paper
mdColours <- c("#892A81", "#84CA69", "#FF5720", "#C1BBB1", "#004A20")
mdTextColours <- c("white", "black", "black", "black", "white")

# Domain border segments (plot space: x = coordY, y = coordX) wherever a 4-neighbour is another domain or none
md_key <- paste(mdDomains$coordX, mdDomains$coordY)
mdBorders <- unique(do.call(rbind, lapply(list(c(0, 1), c(0, -1), c(1, 0), c(-1, 0)), function(d) {
  nb <- mdDomains$domain[match(paste(mdDomains$coordX + d[1], mdDomains$coordY + d[2]), md_key)]
  b <- mdDomains[is.na(nb) | nb != mdDomains$domain, ]
  data.frame(
    x = b$coordY + ifelse(d[2] == 0, -0.5, d[2] / 2),
    xend = b$coordY + ifelse(d[2] == 0, 0.5, d[2] / 2),
    y = b$coordX + ifelse(d[1] == 0, -0.5, d[1] / 2),
    yend = b$coordX + ifelse(d[1] == 0, 0.5, d[1] / 2)
  )
})))

# One label per domain, on its bin closest to the domain centroid
mdLabels <- do.call(rbind, lapply(split(mdDomains, mdDomains$domain), function(d) {
  i <- which.min((d$coordX - mean(d$coordX))^2 + (d$coordY - mean(d$coordY))^2)
  data.frame(domain = d$domain[1], coordX = d$coordX[i], coordY = d$coordY[i])
}))

# --- Helper Functions ---

# Add the morphodynamic-domain and macrochaetae layers on top of a notum map.
addLandmarks <- function(g, showMacrochaetae, showMD, borderColour = "black") {
  if (showMD) {
    g <- g +
      geom_segment(
        data = mdBorders, aes(x = x, y = y, xend = xend, yend = yend),
        inherit.aes = FALSE, colour = borderColour, linewidth = 0.5
      ) +
      geom_label(
        data = mdLabels, aes(coordY, coordX, label = paste0("5-MD", domain)),
        inherit.aes = FALSE, fill = mdColours[mdLabels$domain],
        colour = mdTextColours[mdLabels$domain], size = 2.5, fontface = "bold"
      )
  }
  if (showMacrochaetae) {
    g <- g + geom_point(
      data = macrochaetaePosition, aes(coordY, coordX), inherit.aes = FALSE,
      shape = 21, color = "black", fill = "yellow", size = 4.5
    )
  }
  g
}

# Plot a gene expression heatmap (white-to-red) on the notum spatial grid.
plotMap <- function(inputGene, dataset, coordinates, showMacrochaetae = TRUE, showMD = FALSE) {
  df <- data.frame(gene_data = dataset[, inputGene], coordinates)
  df <- df[!is.na(df[, 1]), ]
  g <- ggplot(data = df, aes(x = df[, 3], y = df[, 2])) +
    geom_tile(aes(fill = df[, 1])) +
    coord_fixed(ratio = 1) +
    scale_y_reverse() +
    scale_fill_gradient2(low = "white", high = "#FF0000",
      name = "Normalized expression (AU)") +
    theme_void() +
    theme(legend.position = "right",
      legend.title = element_text(angle = 90, hjust = 0.5, size = 14),
      legend.text = element_text(size = 11),
      legend.key.height = unit(1.5, "cm"),
      legend.key.width = unit(0.4, "cm")) +
    guides(fill = guide_colorbar(title.position = "right", title.hjust = 0.5))
  addLandmarks(g, showMacrochaetae, showMD)
}

# Plot a two-color (red/green RGB) overlay of two gene expression patterns.
plotMapOverlay <- function(inputGene1, inputGene2, dataset, coordinates, showMacrochaetae = TRUE, showMD = FALSE) {
  gene_data <- dataset[, c(inputGene1, inputGene2)]
  df <- as.data.frame(cbind(gene_data, coordinates))
  keep <- !is.na(df[, 1]) & !is.na(df[, 2])
  df <- as.data.frame(df[keep, ])
  df$rgb <- rgb(df[, 1] / max(df[, 1]), df[, 2] / max(df[, 2]), 0)
  g <- ggplot(df, aes(x = coordY, y = coordX, fill = rgb)) +
    geom_raster() +
    coord_fixed(ratio = 1) +
    scale_y_reverse() +
    scale_fill_identity() +
    theme_void()
  addLandmarks(g, showMacrochaetae, showMD, borderColour = "white")
}

# Compute Pearson correlation between a user-drawn ROI and all genes in the dataset.
corrWithROI <- function(draw.data, ROI, dataset) {
  out <- draw.data
  selection <- out$ID %in% ROI
  out$expression[selection] <- 1
  if (sum(out$expression) < 10) return(NULL)
  correlations <- apply(
    data.frame(correlation = out$expression),
    MARGIN = 2, FUN = cor, dataset
  )
  df <- data.frame(gene = colnames(dataset), correlation = round(correlations, 2))
  df <- df[order(df$correlation, decreasing = T), ]
  row.names(df) <- NULL
  df
}

# Compute Pearson correlation between a gene of interest and all other genes.
corrWithGOI <- function(inputGene, dataset) {
  correlations <- apply(
    data.frame(correlation = dataset[, inputGene]),
    MARGIN = 2, FUN = cor, dataset
  )
  df <- data.frame(gene = colnames(dataset), correlation = round(correlations, 2), row.names = NULL)
  df
}

# --- UI ---

ui <- shinyUI(tagList(
  useShinyjs(),
  fluidPage(
    theme = shinytheme("cosmo"),
    navbarPage("Notum vGEP", selected = "About",

      # Atlas tab
      tabPanel("Atlas",
        tags$head(tags$style(
          type = "text/css",
          "#imageAtlas img {max-width: 100%; width: 100%; height: auto}"
        )),
        sidebarLayout(
          sidebarPanel(
            selectInput("atlasGene", "Atlas Gene", choices = names_atlas, selected = "wg"),
            switchInput(inputId = "toggleAtlas", value = TRUE, label = "Macrochaetae"),
            switchInput(inputId = "toggleAtlasMD", value = FALSE, label = "MD domains")
          ),
          mainPanel(imageOutput("imageAtlas"))
        )
      ),

      # Validation tab
      tabPanel("Validation",
        tags$head(tags$style(
          type = "text/css",
          "#imageValidation img {max-width: 100%; width: 100%; height: auto}"
        )),
        sidebarLayout(
          sidebarPanel(
            selectInput("ValidationGene", "Validation Gene",
              choices = names_validation, selected = "fj (RNA)"),
            switchInput(inputId = "toggleValidation", value = TRUE, label = "Macrochaetae"),
            switchInput(inputId = "toggleValidationMD", value = FALSE, label = "MD domains")
          ),
          mainPanel(imageOutput("imageValidation"))
        )
      ),

      # vGEP tab
      tabPanel("vGEP",
        tags$head(tags$style(
          type = "text/css",
          "#imagevGEP img {max-width: 100%; width: 100%; height: auto}"
        )),
        sidebarLayout(
          sidebarPanel(
            selectizeInput("vGEP", "vGEP expression pattern",
              choices = names_vGEP, selected = "dpp",
              options = list(maxOptions = 17000)),
            switchInput(inputId = "togglevGEP", value = TRUE, label = "Macrochaetae"),
            switchInput(inputId = "togglevGEPMD", value = FALSE, label = "MD domains")
          ),
          mainPanel(imageOutput("imagevGEP"))
        )
      ),

      # Correlate tab
      tabPanel("Correlate",
        tags$head(tags$style(
          type = "text/css",
          "#imageCorrGOI img {max-width: 100%; width: 100%; height: auto}",
          "#imageCorrSelectedFromTable img {max-width: 100%; width: 100%; height: auto}",
          "#comp2vGEP img {max-width: 100%; width: 100%; height: auto}",
          "#correlateTextGOI {font-size: 20px}",
          "#correlateTextSelectedFromTable {font-size: 20px}"
        )),
        fluidRow(
          column(width = 4,
            wellPanel(
              selectizeInput("Correlate", "vGEP",
                choices = names_vGEP, selected = "dpp",
                options = list(maxOptions = 17000))
            ),
            DTOutput("corrWithvGEP"),
            downloadButton("corrWithvGEPDownload", "Download Table")
          ),
          column(width = 4,
            htmlOutput("correlateTextGOI", align = "center"),
            imageOutput("imageCorrGOI")
          ),
          column(width = 4,
            htmlOutput("correlateTextSelectedFromTable", align = "center"),
            imageOutput("imageCorrSelectedFromTable")
          )
        )
      ),

      # Draw tab
      tabPanel("Draw",
        tags$head(tags$style(
          type = "text/css",
          "#DrawPlot img {max-width: 100%; width: 100%; height: 100%}"
        )),
        titlePanel(
          h3("Draw a pattern to discover gene expression in ROI", align = "center")),
        h4(em("Select region by clicking on dots or by using the lasso"), align = "center"),
        fluidRow(
          column(width = 7,
            actionButton("reset", label = "Reset selection"),
            hidden(downloadButton("DrawDTdownload", "Download Table")),
            div(style = "display: inline-block; vertical-align: top; margin-left: 20px",
              switchInput(inputId = "toggleDraw", value = TRUE, label = "Macrochaetae")),
            div(style = "display: inline-block; vertical-align: top",
              switchInput(inputId = "toggleDrawMD", value = FALSE, label = "MD domains")),
            ggiraph::girafeOutput("DrawPlot"),
            h5(em(
              "Select at least 10 pixels. Selected pixels will turn red.",
              br(),
              "Use blue lasso tool to select pixels and red lasso to unselect.",
              br(),
              "Lasso tools will become visible when the mouse is moved over the plot."
            ))
          ),
          column(width = 5,
            h5("Correlation (high-low)", align = "center"),
            DTOutput("corrWithROI")
          )
        )
      ),

      # Compare tab
      tabPanel("Compare",
        tags$head(tags$style(
          type = "text/css",
          "#comp1img img {max-width: 100%; width: 100%; height: auto}",
          "#comp2img img {max-width: 100%; width: 100%; height: auto}",
          "#comp1name {font-size: 20px}",
          "#comp2name {font-size: 20px}"
        )),
        fluidRow(
          column(width = 4,
            wellPanel(
              h4(strong("Which genes would you like to compare?")),
              br(),
              radioButtons("toCompare1", "Compare", choices = c("Atlas", "Validation", "vGEP")),
              radioButtons("toCompare2", "To", choices = c("Atlas", "Validation", "vGEP")),
              hidden(selectInput("atlas1", "Atlas", choices = names_atlas, selected = "sr")),
              hidden(selectInput("atlas2", "Atlas", choices = names_atlas, selected = "hh")),
              hidden(selectInput("validation1", "Validation",
                choices = names_validation, selected = "dpp (RNA)")),
              hidden(selectInput("validation2", "Validation",
                choices = names_validation, selected = "nw (RNA)")),
              hidden(selectizeInput("vGEP1", "vGEP", choices = names_vGEP, selected = "hh", options = list(maxOptions = 17000))),
              hidden(selectizeInput("vGEP2", "vGEP", choices = names_vGEP, selected = "dpp", options = list(maxOptions = 17000)))
            )
          ),
          column(width = 4,
            htmlOutput("comp1name", align = "center"),
            imageOutput("comp1img")
          ),
          column(width = 4,
            htmlOutput("comp2name", align = "center"),
            imageOutput("comp2img")
          )
        )
      ),

      # Overlay tab
      tabPanel("Overlay",
        tags$head(tags$style(
          type = "text/css",
          "#overlay_vGEP img {max-width: 100%; width: 100%; height: auto}",
          "#overlay_protein img {max-width: 100%; width: 100%; height: auto}",
          "#overlay_validation img {max-width: 100%; width: 100%; height: auto}"
        )),
        sidebarLayout(
          sidebarPanel(
            h4(strong("What would you like to overlay?")),
            br(),
            radioButtons("overlay", "Overlay", choices = c("Atlas", "Validation", "vGEP")),
            hidden(
              selectInput("overlay_atlas1", "Gene 1 (Red)", choices = names_atlas, selected = "hh"),
              selectInput("overlay_atlas2", "Gene 2 (Green)", choices = names_atlas, selected = "wg")
            ),
            hidden(
              selectInput("overlay_validation1", "Gene 1 (Red)",
                choices = names_validation, selected = "dpp (RNA)"),
              selectInput("overlay_validation2", "Gene 2 (Green)",
                choices = names_validation, selected = "Dl (protein)")
            ),
            hidden(
              selectizeInput("overlay_vGEP1", "Gene 1 (Red)", choices = names_vGEP, selected = "dpp", options = list(maxOptions = 17000)),
              selectizeInput("overlay_vGEP2", "Gene 2 (Green)", choices = names_vGEP, selected = "hh", options = list(maxOptions = 17000))
            ),
            switchInput(inputId = "toggleOverlay", value = TRUE, label = "Macrochaetae"),
            switchInput(inputId = "toggleOverlayMD", value = FALSE, label = "MD domains")
          ),
          mainPanel(
            conditionalPanel(condition = "input.overlay == 'Atlas'",
              imageOutput("overlay_protein")),
            conditionalPanel(condition = "input.overlay == 'Validation'",
              imageOutput("overlay_validation")),
            conditionalPanel(condition = "input.overlay == 'vGEP'",
              imageOutput("overlay_vGEP"))
          )
        )
      ),

      # About tab
      tabPanel("About",
        h4(strong("About this app")),
        p("The Notum vGEP application was developed to visualize and explore spatial gene expression data in the",
          em("Drosophila melanogaster"), "pupal notum epithelium (see Figure). All images reflect the right hemi-notum where yellow dots reflect landmark macrochaetae. The \"MD domains\" switch outlines the five morphodynamic domains (5-MD1–5). To generate the spatial transcriptome of the notum tissue, we first constructed a reference atlas based on the expression patterns of 24 proteins. Using vGEP, we then mapped scRNA-seq data onto this atlas, resulting in a spatial transcriptome comprising 4,312 gene expression profiles."),
        p("The accuracy of the spatial reconstruction was validated using expression profiles from an additional set of 7 proteins and 13 RNA target genes (measured by smiFISH), plus 3 additional protein patterns (Tollo, yrt, stan)."),
        br(),
        p(em("Note:"), "Atlas and Validation expression profiles represent average profiles derived from multiple animals."),
        p("The scRNA-seq data are available under GEO accession number GSE328433."),
        br(),
        h4(strong("This app contains")),
        tags$ul(
          tags$li(strong("Atlas"), " \u2014 The user can enter a gene name or browse the list of alphabetically ordered atlas gene expression profiles (24 protein patterns)."),
          tags$li(strong("Validation"), " \u2014 The user can enter a gene name or browse the list of alphabetically ordered validation gene expression profiles (23 patterns)."),
          tags$li(strong("vGEP"), " \u2014 The user can enter a gene name or browse the list of alphabetically ordered vGEP-predicted gene expression profiles (4,312 RNA patterns)."),
          tags$li(strong("Correlate"), " \u2014 The user can enter a gene name or browse an alphabetically ordered list of vGEP-predicted expression profiles, then click on a listed gene to identify correlated gene expression profiles."),
          tags$li(strong("Draw"), " \u2014 The user can define a ROI within the tissue to list the vGEP-predicted gene expression profiles ranked by their correlation with this ROI."),
          tags$li(strong("Compare"), " \u2014 The user can enter two gene names or browse an alphabetically ordered list of atlas, validation, and vGEP-predicted expression profiles to compare side by side their expression profiles."),
          tags$li(strong("Overlay"), " \u2014 The user can select two genes from atlas, validation, or vGEP-predicted expression profiles to visualize their co-expression as a two-color (red/green) overlay.")
        ),
        br(),
        img(src = "cover_image.png", width = "80%", align = "center"),
        br(), br(),
        h4(strong("Funding")),
        p("This work was supported by Institut Curie, CNRS, INSERM, ERC Advanced TIMORPH (340784) and Scaling-Sensitivity (101020243), CANCERO-INCA (PLBIO2020/BELLAICHE), ANR Labex DEEP (11-LBX-0044, ANR-10-IDEX-0001-02), Marie Sklodowska-Curie Innovative Training Network 'PolarNet' (675407) and FRM (FDT201904008163)."),
        br(),
        h4(strong("Contact")),
        p("For further information contact Yohanns Bellaïche: yohanns.bellaiche@curie.fr."),
        br(),
        h4(strong("Reference")),
        p("Please acknowledge the original publication when using this website or its data:"),
        p("Leroy, Van Leen et al., Cell Reports (2026). Code and data: ",
          a("10.5281/zenodo.19633285", href = "https://doi.org/10.5281/zenodo.19633285", target = "_blank"))
      )
    )
  )
))

# --- Server ---

server <- function(input, output, session) {

  # Atlas tab
  output$imageAtlas <- renderPlot({
    validate(need(input$atlasGene, "choose a gene"))
    plotMap(input$atlasGene, AtlasData, Atlascoord, input$toggleAtlas, input$toggleAtlasMD)
  }, height = 1600, width = 2400, res = 300)

  # Validation tab
  output$imageValidation <- renderPlot({
    validate(need(input$ValidationGene, "choose a gene"))
    plotMap(input$ValidationGene, ValidationData, ValidationCoord, input$toggleValidation, input$toggleValidationMD)
  }, height = 1600, width = 2400, res = 300)

  # vGEP tab
  output$imagevGEP <- renderPlot({
    validate(need(input$vGEP, "choose a gene"))
    plotMap(input$vGEP, vGEPData, vGEPcoord, input$togglevGEP, input$togglevGEPMD)
  }, height = 1600, width = 2400, res = 300)

  # Correlate tab
  corrOutput <- reactive({
    validate(need(input$Correlate, "choose a gene"))
    corrWithGOI(input$Correlate, vGEPData)
  })

  output$corrWithvGEP <- DT::renderDT({
    DT::datatable(corrOutput(),
      selection = "single", rownames = F,
      options = list(
        server = F,
        columnDefs = list(list(className = "dt-center", targets = 0:1)),
        pageLength = 10, responsive = TRUE, autoWidth = TRUE,
        order = list(list(1, "desc"))
      )
    )
  })

  selected_correlated_gene <- reactive({
    input$corrWithvGEP_rows_selected
  })

  output$imageCorrGOI <- renderPlot({
    validate(need(input$Correlate, ""))
    validate(need(selected_correlated_gene(), "Select gene in the list to compare"))
    plotMap(input$Correlate, vGEPData, vGEPcoord)
  }, height = 1600, width = 2400, res = 300)

  output$imageCorrSelectedFromTable <- renderPlot({
    validate(need(input$Correlate, ""))
    validate(need(selected_correlated_gene(), ""))
    plotMap(selected_correlated_gene(), vGEPData, vGEPcoord)
  }, height = 1600, width = 2400, res = 300)

  output$correlateTextGOI <- renderText({
    validate(need(input$Correlate, ""))
    validate(need(selected_correlated_gene(), ""))
    paste0("<b>", input$Correlate, "</b>")
  })

  output$correlateTextSelectedFromTable <- renderText({
    validate(need(selected_correlated_gene(), ""))
    paste0("<b>", colnames(vGEPData)[selected_correlated_gene()], "</b>")
  })

  output$corrWithvGEPDownload <- downloadHandler(
    filename = function() {
      paste("data-", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      write.csv(corrOutput(), file)
    }
  )

  # Draw tab
  selected_ROI <- reactive({
    input$DrawPlot_selected
  })

  output$DrawPlot <- renderGirafe({
    # Re-rendering on a toggle keeps the current selection
    girafe(
      code = print(addLandmarks(draw.selection, input$toggleDraw, input$toggleDrawMD)),
      width_svg = 8, height_svg = 6,
      options = list(opts_selection(
        selected = isolate(input$DrawPlot_selected),
        type = "multiple", css = "fill:#FF3333;stroke:black;"
      ))
    )
  })

  observeEvent(input$reset, {
    session$sendCustomMessage(type = "DrawPlot_set", message = character(0))
  })

  outputCorrWithROI <- reactive({
    corrWithROI(draw.data, selected_ROI(), vGEPData)
  })

  output$corrWithROI <- DT::renderDT({
    DT::datatable(outputCorrWithROI(),
      selection = "none",
      options = list(
        columnDefs = list(list(className = "dt-center", targets = 0:2)),
        pageLength = 10, responsive = TRUE, autoWidth = TRUE,
        order = list(list(2, "desc"))
      )
    )
  })

  output$DrawDTdownload <- downloadHandler(
    filename = function() {
      paste("data-", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      write.csv(outputCorrWithROI(), file)
    }
  )

  # Compare tab
  output$comp1name <- renderText({
    if (input$toCompare1 == "Atlas") {
      paste0("Atlas Gene ", "<b>", input$atlas1, "</b>")
    } else if (input$toCompare1 == "Validation") {
      paste0("Validation Gene ", "<b>", input$validation1, "</b>")
    } else {
      paste0("Predicted Spatial Expression for ", "<b>", input$vGEP1, "</b>")
    }
  })

  output$comp2name <- renderText({
    if (input$toCompare2 == "Atlas") {
      paste0("Atlas Gene ", "<b>", input$atlas2, "</b>")
    } else if (input$toCompare2 == "Validation") {
      paste0("Validation Gene ", "<b>", input$validation2, "</b>")
    } else {
      paste0("Predicted Spatial Expression for ", "<b>", input$vGEP2, "</b>")
    }
  })

  output$comp1img <- renderPlot({
    if (input$toCompare1 == "Atlas") {
      plotMap(input$atlas1, AtlasData, Atlascoord)
    } else if (input$toCompare1 == "Validation") {
      plotMap(input$validation1, ValidationData, ValidationCoord)
    } else {
      plotMap(input$vGEP1, vGEPData, vGEPcoord)
    }
  }, height = 1600, width = 2400, res = 300)

  output$comp2img <- renderPlot({
    if (input$toCompare2 == "Atlas") {
      plotMap(input$atlas2, AtlasData, Atlascoord)
    } else if (input$toCompare2 == "Validation") {
      plotMap(input$validation2, ValidationData, ValidationCoord)
    } else {
      plotMap(input$vGEP2, vGEPData, vGEPcoord)
    }
  }, height = 1600, width = 2400, res = 300)

  # Overlay tab
  output$overlay_vGEP <- renderPlot({
    plotMapOverlay(input$overlay_vGEP1, input$overlay_vGEP2,
      vGEPData, vGEPcoord, input$toggleOverlay, input$toggleOverlayMD)
  }, height = 1600, width = 2400, res = 300)

  output$overlay_protein <- renderPlot({
    plotMapOverlay(input$overlay_atlas1, input$overlay_atlas2,
      AtlasData, Atlascoord, input$toggleOverlay, input$toggleOverlayMD)
  }, height = 1600, width = 2400, res = 300)

  output$overlay_validation <- renderPlot({
    plotMapOverlay(input$overlay_validation1, input$overlay_validation2,
      ValidationData, ValidationCoord, input$toggleOverlay, input$toggleOverlayMD)
  }, height = 1600, width = 2400, res = 300)

  # Conditional UI elements (shinyjs toggles)
  observe({
    toggleElement(id = "atlas1", condition = input$toCompare1 == "Atlas")
    toggleElement(id = "atlas2", condition = input$toCompare2 == "Atlas")
    toggleElement(id = "validation1", condition = input$toCompare1 == "Validation")
    toggleElement(id = "validation2", condition = input$toCompare2 == "Validation")
    toggleElement(id = "vGEP1", condition = input$toCompare1 == "vGEP")
    toggleElement(id = "vGEP2", condition = input$toCompare2 == "vGEP")
    toggleElement(id = "overlay_atlas1", condition = input$overlay == "Atlas")
    toggleElement(id = "overlay_atlas2", condition = input$overlay == "Atlas")
    toggleElement(id = "overlay_vGEP1", condition = input$overlay == "vGEP")
    toggleElement(id = "overlay_vGEP2", condition = input$overlay == "vGEP")
    toggleElement(id = "overlay_validation1", condition = input$overlay == "Validation")
    toggleElement(id = "overlay_validation2", condition = input$overlay == "Validation")
    toggleElement("DrawDTdownload", condition = !is.null(outputCorrWithROI()))
  })
}

shinyApp(ui = ui, server = server)
