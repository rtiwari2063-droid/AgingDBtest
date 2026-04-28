# ==========================================
# 1. LOAD REQUIRED LIBRARIES
# ==========================================
library(shiny)
library(bslib)
library(DT)
library(ggplot2)
library(dplyr)
library(tidyr)
library(visNetwork) 

# ==========================================
# 2. PRE-LOAD HEAVY DATA (FETCH FROM GITHUB)
# ==========================================
# This helper safely downloads large binary files into WebAssembly memory
fetch_rds_from_github <- function(file_name) {
  # Points directly to your raw GitHub repository files
  base_url <- "https://raw.githubusercontent.com/rtiwari2063-droid/AgingDBtest/main/"
  file_url <- paste0(base_url, file_name)
  
  if(!file.exists(file_name)) {
    print(paste("Downloading", file_name, "from GitHub..."))
    # mode = "wb" ensures the binary .rds file is not corrupted during download
    download.file(file_url, destfile = file_name, mode = "wb", quiet = TRUE)
  }
  return(readRDS(file_name))
}

# ---------------------------------------------------------
# A. LOAD OPTIMIZED UNIVERSE & DICTIONARIES
# ---------------------------------------------------------
print("Setting up Background Universe...")
universe_df <- fetch_rds_from_github("optimized_universe.rds")

universe_genes <- unique(universe_df$ENTREZID)
N_total_db_genes <- length(universe_genes)
entrez_to_symbol <- setNames(universe_df$SYMBOL, universe_df$ENTREZID)
symbol_to_entrez <- setNames(universe_df$ENTREZID, universe_df$SYMBOL)

# ---------------------------------------------------------
# B. LOAD THE SINGLE MASTER DATABASE
# ---------------------------------------------------------
print("Loading Master Database...")
master_db <- fetch_rds_from_github("Master_Hallmark_Wide.rds")

unique_hallmarks <- sort(unique(as.character(master_db$HallMarks)))
combined_exo_genes <- master_db %>% dplyr::filter(In_Exosome == "Yes") %>% dplyr::pull(Gene_ID) %>% unique()
reserved_cols <- c("Gene_ID", "HallMarks", "In_Exosome", "Consensus", "Extended", "evidence", "no_evidence", "Whole Body")
cell_class_cols <- setdiff(colnames(master_db), reserved_cols)
cell_class_choices <- setNames(cell_class_cols, tools::toTitleCase(gsub("_", " ", cell_class_cols)))

# ---------------------------------------------------------
# C. LOAD OPTIMIZED NETWORKS
# ---------------------------------------------------------
print("Loading Optimized Network Data...")
db_tf <- fetch_rds_from_github("optimized_tf_network.rds")
db_kinase <- fetch_rds_from_github("optimized_kinase_network.rds")

# ---------------------------------------------------------
# D. TISSUE MAPPING (Hardcoded fallback for speed)
# ---------------------------------------------------------
tissue_mapping_df <- data.frame(Tissue = character(), Cell_Type_Class = character(), stringsAsFactors = FALSE)
available_tissues <- unique(na.omit(tissue_mapping_df$Tissue))
available_tissues <- available_tissues[available_tissues != ""]

# ==========================================
# 3. CUSTOM THEME & CSS 
# ==========================================
my_theme <- bs_theme(
  version = 5, primary = "#0EA5E9", info = "#3B82F6", base_font = font_google("Inter"), "navbar-bg" = "#1E3A8A" 
)

custom_css <- tags$head(tags$style(HTML("
  .navbar-brand { font-size: 1.6rem !important; font-weight: 700; color: #ffffff !important; }
  .nav-link { font-weight: 500; font-size: 1.05rem; color: #e2e8f0 !important; }
  .nav-link:hover { color: #ffffff !important; }
  .nav-link.active { border-bottom: 3px solid #38BDF8 !important; color: #ffffff !important; }
  .dataTables_wrapper { padding: 10px; }
  table.dataTable.no-footer { border-bottom: none !important; }
  .plot-container { background-color: #ffffff !important; border-radius: 8px; padding: 20px; box-shadow: 0 4px 6px -1px rgba(0,0,0,0.1); color: #000000 !important; }
  .filter-bar { background-color: var(--bs-tertiary-bg); border: 1px solid var(--bs-border-color); border-radius: 8px; padding: 20px 20px 5px 20px; margin-bottom: 25px; }
  .summary-box { background-color: var(--bs-secondary-bg); color: var(--bs-body-color); border-left: 4px solid #3B82F6; padding: 15px; border-radius: 6px; margin-bottom: 20px; }
  
  .hallmark-box { background: linear-gradient(145deg, #1e293b, #3b82f6); color: #f8fafc; border-radius: 16px; padding: 30px 15px; width: 100%; min-height: 190px; border: 1px solid rgba(255,255,255,0.1); box-shadow: 0 10px 25px -5px rgba(0, 0, 0, 0.1), 0 8px 10px -6px rgba(0, 0, 0, 0.1); transition: all 0.3s cubic-bezier(0.4, 0, 0.2, 1); display: flex; flex-direction: column; align-items: center; justify-content: center; text-align: center; margin-bottom: 20px; }
  .hallmark-box:hover { transform: translateY(-8px); box-shadow: 0 20px 25px -5px rgba(59, 130, 246, 0.4), 0 8px 10px -6px rgba(59, 130, 246, 0.1); background: linear-gradient(145deg, #1e40af, #60a5fa); color: white; }
  .hallmark-custom-icon { width: 64px; height: 64px; margin-bottom: 18px; object-fit: contain; transition: transform 0.3s ease; background-color: rgba(255,255,255,0.1); border-radius: 8px; }
  .hallmark-box:hover .hallmark-custom-icon { transform: scale(1.15); }
  .hallmark-text { font-size: 1.15rem; font-weight: 700; line-height: 1.3; letter-spacing: 0.3px; }
")))

# ==========================================
# 4. DEFINE THE USER INTERFACE (UI)
# ==========================================
ui <- page_navbar(
  id = "main_nav",
  title = span(icon("dna", style="color: #38BDF8; margin-right: 12px;"), "AgingHallmarkDB"),
  theme = my_theme, fillable = FALSE, header = custom_css, 
  
  nav_spacer(),
  nav_item(input_dark_mode(id = "dark_mode")), 
  
  nav_panel(title = "Browse Hallmarks", value = "tab_home",
            conditionalPanel(condition = "output.home_page === 'grid'", div(class="container text-center mt-5 mb-5", style = "max-width: 1400px;", h2("Browse Genes by Aging Hallmark", class="fw-bold mb-5"), uiOutput("hallmark_grid_ui"))),
            conditionalPanel(condition = "output.home_page === 'result'", div(class="container mt-4 mb-5", style = "max-width: 1000px;", actionButton("back_home_btn", "← Back to Grid", class="btn btn-outline-secondary fw-bold mb-4"), card(class = "shadow-sm border-0 mb-4", card_body(h4("Filter Selection", class="text-primary fw-bold mb-3"), layout_columns(col_widths = c(4, 4, 4), selectInput("home_db_choice", label = "Gene Set Database:", width = "100%", choices = list("Extended" = "Extended", "Consensus" = "Consensus")), radioButtons("home_ev_choice", "Evidence Filter:", inline = TRUE, choices = list("With Evidence" = "evidence", "No Evidence" = "no_evidence")), radioButtons("home_spec_level", "Specificity Level:", inline = TRUE, choices = list("Whole Body" = "all", "By Tissue" = "tissue", "By Cell Type Class" = "cell"))), conditionalPanel("input.home_spec_level == 'tissue'", selectizeInput("home_tissue_choice", "Select Tissue(s):", width = "100%", choices = available_tissues, multiple = TRUE), uiOutput("home_tissue_cell_ui")), conditionalPanel("input.home_spec_level == 'cell'", selectizeInput("home_cell_choice", "Select Cell Type(s):", width = "100%", choices = cell_class_choices, multiple = TRUE)))), card(class = "shadow-sm border-0", card_header(class = "bg-transparent", h3(textOutput("selected_hallmark_title"), class="text-primary fw-bold m-0")), card_body(uiOutput("home_error_msg"), DTOutput("hallmark_genes_table")))))
  ),
  
  nav_panel(title = "Hallmark Enrichment", value = "tab_enrich_inputs",
            conditionalPanel(condition = "output.enrich_page === 'input'", div(class = "container mt-5 mb-5", style = "max-width: 900px;", card(class = "shadow-sm border-0", card_body(h3("Hallmarks Enrichment", class = "text-primary fw-bold"), p("Identify significantly overrepresented aging hallmarks in your input gene set.", class="text-muted"), hr(), h4(tags$span(class = "badge bg-primary rounded-circle", "1"), " Enter gene list (Entrez IDs)", class = "fw-bold mt-2"), div(class = "mb-2", actionButton("sample_btn", "Load sample genes", class = "btn btn-sm btn-info text-white fw-bold"), actionButton("reset_btn", "Reset list", class = "btn btn-outline-secondary btn-sm ms-2")), textAreaInput("gene_input", label = NULL, height = "150px", placeholder = "e.g. 7007\n9370\n652"), h4(tags$span(class = "badge bg-primary rounded-circle", "2"), " Select parameters", class = "fw-bold mt-4"), layout_columns(col_widths = c(6, 6), selectInput("db_choice", label = "Gene Set Database:", width = "100%", choices = list("Extended" = "Extended", "Consensus" = "Consensus")), radioButtons("enrich_ev_choice", "Evidence Filter:", inline = TRUE, choices = list("With Evidence" = "evidence", "No Evidence" = "no_evidence"))), radioButtons("enrich_spec_level", "Specificity Level:", inline = TRUE, choices = list("Whole Body" = "all", "By Tissue" = "tissue", "By Cell Type Class" = "cell")), conditionalPanel("input.enrich_spec_level == 'tissue'", selectizeInput("enrich_tissue_choice", "Select Tissue(s):", width = "100%", choices = available_tissues, multiple = TRUE), uiOutput("enrich_tissue_cell_ui")), conditionalPanel("input.enrich_spec_level == 'cell'", selectizeInput("enrich_cell_choice", "Select Cell Type(s):", width = "100%", choices = cell_class_choices, multiple = TRUE)), h4(tags$span(class = "badge bg-primary rounded-circle", "3"), " Submit", class = "fw-bold mt-4"), uiOutput("enrich_error_msg"), actionButton("run_btn", "Run analysis →", class = "btn btn-lg btn-primary w-100 fw-bold mt-2"))))),
            conditionalPanel(condition = "output.enrich_page === 'result'", div(class = "container mt-4 mb-5", style = "max-width: 1200px;", layout_columns(col_widths = c(6, 6), actionButton("back_enrich_btn", "← Back to Inputs", class = "btn btn-outline-secondary fw-bold"), div(style = "text-align: right;", downloadButton("dl_enrich_csv", "📥 Download Table", class = "btn btn-success fw-bold"))), h3("Analysis Results", class = "fw-bold mt-4 mb-3 text-primary border-bottom border-primary pb-2"), uiOutput("mapping_summary_ui"), div(class = "filter-bar", checkboxInput("enable_filters", tags$strong("Enable Significance Filters (FDR & Minimum Overlap)"), value = FALSE), conditionalPanel(condition = "input.enable_filters == true", layout_columns(col_widths = c(6, 6), sliderInput("fdr_cutoff", "FDR Threshold:", min = 0.01, max = 0.25, value = 0.05, step = 0.01, width = "100%"), sliderInput("min_overlap", "Minimum Gene Overlap:", min = 1, max = 20, value = 1, step = 1, width = "100%")))), uiOutput("kpi_boxes_enrichment"), br(), card(class = "shadow-sm border-0 mb-4", card_body(DTOutput("enrichmentTable"))), card(class = "shadow-sm border-0", card_header(class = "d-flex justify-content-between align-items-center bg-transparent", radioButtons("plot_type", label = "Visualization:", inline = TRUE, choices = list("Bar plot" = "bar", "Bubble plot" = "bubble", "Circular plot" = "circular")), downloadButton("dl_enrich_plot", "🖼️ Download Image", class = "btn btn-outline-primary fw-bold")), div(class = "plot-container", plotOutput("enrichmentPlot", height = "650px")))))
  ),
  
  nav_menu(title = "Gene Mapping",
           nav_panel(title = "Map to Hallmarks", value = "tab_map_inputs",
                     conditionalPanel(condition = "output.map_page === 'input'", div(class = "container mt-5 mb-5", style = "max-width: 900px;", card(class = "shadow-sm border-0", card_body(h3("Map Gene Sets to Hallmarks", class = "text-primary fw-bold"), hr(), layout_columns(col_widths = c(6, 6), selectInput("mapping_db_choice", label = "1. Select Database:", width = "100%", choices = list("Extended" = "Extended", "Consensus" = "Consensus")), radioButtons("map_ev_choice", "Evidence Filter:", inline = TRUE, choices = list("With Evidence" = "evidence", "No Evidence" = "no_evidence"))), radioButtons("map_spec_level", "2. Specificity Level:", inline = TRUE, choices = list("Whole Body" = "all", "By Tissue" = "tissue", "By Cell Type Class" = "cell")), conditionalPanel("input.map_spec_level == 'tissue'", selectizeInput("map_tissue_choice", "Select Tissue(s):", width = "100%", choices = available_tissues, multiple = TRUE), uiOutput("map_tissue_cell_ui")), conditionalPanel("input.map_spec_level == 'cell'", selectizeInput("map_cell_choice", "Select Cell Type(s):", width = "100%", choices = cell_class_choices, multiple = TRUE)), div(class="d-flex justify-content-between align-items-end mt-3", tags$label("3. Paste Gene IDs:", class="fw-bold mb-0"), div(actionButton("sample_map_btn", "Load sample", class = "btn btn-sm btn-info text-white"), actionButton("reset_map_btn", "Reset", class = "btn btn-outline-secondary btn-sm ms-2"))), textAreaInput("mapping_input", label = NULL, height = "220px", placeholder = "Enter Entrez IDs..."), uiOutput("mapping_error_msg"), actionButton("run_mapping_btn", "Map to Hallmarks →", class = "btn btn-lg btn-primary w-100 fw-bold mt-3"))))),
                     conditionalPanel(condition = "output.map_page === 'result'", div(class = "container mt-4 mb-5", style = "max-width: 1000px;", layout_columns(col_widths = c(6, 6), actionButton("back_map_btn", "← Back to Inputs", class = "btn btn-outline-secondary fw-bold"), div(style = "text-align: right;", downloadButton("dl_map_csv", "📥 Download Mapping", class = "btn btn-success fw-bold"))), h3("Gene Mapping Results", class = "fw-bold mt-4 mb-3 text-primary border-bottom border-primary pb-2"), uiOutput("kpi_boxes_mapping"), br(), card(class = "shadow-sm border-0", card_body(DTOutput("mappingTable")))))
           ),
           nav_panel(title = "Exosome Presence Check", value = "tab_exo_inputs",
                     conditionalPanel(condition = "output.exo_page === 'input'", div(class = "container mt-5 mb-5", style = "max-width: 900px;", card(class = "shadow-sm border-0", card_body(h3("Exosome Presence Check", class = "text-primary fw-bold"), p("Check if your genes are present in exosome datasets.", class="text-muted"), hr(), layout_columns(col_widths = c(6, 6), selectInput("exo_db_choice", label = "1. Select Hallmark Database:", width = "100%", choices = list("Extended" = "Extended", "Consensus" = "Consensus")), radioButtons("exo_ev_choice", "Evidence Filter:", inline = TRUE, choices = list("With Evidence" = "evidence", "No Evidence" = "no_evidence"))), radioButtons("exo_spec_level", "2. Specificity Level:", inline = TRUE, choices = list("Whole Body" = "all", "By Tissue" = "tissue", "By Cell Type Class" = "cell")), conditionalPanel("input.exo_spec_level == 'tissue'", selectizeInput("exo_tissue_choice", "Select Tissue(s):", width = "100%", choices = available_tissues, multiple = TRUE), uiOutput("exo_tissue_cell_ui")), conditionalPanel("input.exo_spec_level == 'cell'", selectizeInput("exo_cell_choice", "Select Cell Type(s):", width = "100%", choices = cell_class_choices, multiple = TRUE)), div(class="d-flex justify-content-between align-items-end mt-4", tags$label("3. Paste Gene IDs (Entrez):", class="fw-bold mb-0"), div(actionButton("sample_exo_btn", "Load sample", class = "btn btn-sm btn-info text-white"), actionButton("reset_exo_btn", "Reset", class = "btn btn-outline-secondary btn-sm ms-2"))), textAreaInput("exo_input", label = NULL, height = "180px", placeholder = "Enter Entrez IDs..."), uiOutput("exo_error_msg"), actionButton("run_exo_btn", "Check Exosomes →", class = "btn btn-lg btn-primary w-100 fw-bold mt-3"))))),
                     conditionalPanel(condition = "output.exo_page === 'result'", div(class = "container mt-4 mb-5", style = "max-width: 1000px;", layout_columns(col_widths = c(6, 6), actionButton("back_exo_btn", "← Back to Inputs", class = "btn btn-outline-secondary fw-bold"), div(style = "text-align: right;", downloadButton("dl_exo_csv", "📥 Download Table", class = "btn btn-success fw-bold"))), h3("Exosome Presence Results", class = "fw-bold mt-4 mb-3 text-primary border-bottom border-primary pb-2"), uiOutput("exo_summary_ui"), br(), card(class = "shadow-sm border-0", card_body(DTOutput("exoTable")))))
           )
  ),
  
  nav_menu(title = "Network Analysis",
           nav_panel(title = "TF-Target Interactions", value = "tab_tf_net",
                     div(class = "container mt-5 mb-5", style = "max-width: 1200px;", card(class = "shadow-sm border-0", card_body(h3("Transcription Regulatory Network", class = "text-primary fw-bold"), p("Paste your gene list. The tool will map internal regulatory interactions.", class="text-muted"), hr(), layout_columns(col_widths = c(4, 4, 4), div(tags$label("1. Select Database:"), selectInput("tf_db_choice", label = NULL, width = "100%", choices = list("Extended" = "Extended", "Consensus" = "Consensus"))), div(tags$label("Evidence Filter:"), radioButtons("tf_ev_choice", label = NULL, inline = TRUE, choices = list("With Evidence" = "evidence", "No Evidence" = "no_evidence"))), div(tags$label("2. Specificity Level:"), radioButtons("tf_spec_level", label = NULL, inline = TRUE, choices = list("Whole Body" = "all", "By Tissue" = "tissue", "By Cell Type Class" = "cell")))), conditionalPanel("input.tf_spec_level == 'tissue'", selectizeInput("tf_tissue_choice", "Select Tissue(s):", width = "100%", choices = available_tissues, multiple = TRUE), uiOutput("tf_tissue_cell_ui")), conditionalPanel("input.tf_spec_level == 'cell'", selectizeInput("tf_cell_choice", "Select Cell Type(s):", width = "100%", choices = cell_class_choices, multiple = TRUE)), div(class="d-flex justify-content-between align-items-end mt-4", tags$label("3. Paste Gene IDs (Entrez IDs):", class="fw-bold mb-0 text-primary"), div(actionButton("sample_tf_btn", "Load sample", class = "btn btn-sm btn-info text-white fw-bold"), actionButton("reset_tf_btn", "Reset", class = "btn btn-outline-secondary btn-sm ms-2"))), textAreaInput("tf_input", label = NULL, height = "150px", placeholder = "e.g. 7124\n2\n1029\n5970\n2099\n4609"), uiOutput("tf_status_msg"), actionButton("run_tf_btn", "Map Internal TF Network →", class = "btn btn-lg btn-primary w-100 fw-bold mt-3"))), br(), uiOutput("tf_results_panel")
                     )
           ),
           nav_panel(title = "Kinase-Substrate Interactions", value = "tab_kinase_net",
                     div(class = "container mt-5 mb-5", style = "max-width: 1200px;", card(class = "shadow-sm border-0", card_body(h3("Substrate Kinase Network", class = "text-primary fw-bold"), p("Paste your gene list. The tool will map internal phosphorylation interactions.", class="text-muted"), hr(), layout_columns(col_widths = c(4, 4, 4), div(tags$label("1. Select Database:"), selectInput("kin_db_choice", label = NULL, width = "100%", choices = list("Extended" = "Extended", "Consensus" = "Consensus"))), div(tags$label("Evidence Filter:"), radioButtons("kin_ev_choice", label = NULL, inline = TRUE, choices = list("With Evidence" = "evidence", "No Evidence" = "no_evidence"))), div(tags$label("2. Specificity Level:"), radioButtons("kin_spec_level", label = NULL, inline = TRUE, choices = list("Whole Body" = "all", "By Tissue" = "tissue", "By Cell Type Class" = "cell")))), conditionalPanel("input.kin_spec_level == 'tissue'", selectizeInput("kin_tissue_choice", "Select Tissue(s):", width = "100%", choices = available_tissues, multiple = TRUE), uiOutput("kin_tissue_cell_ui")), conditionalPanel("input.kin_spec_level == 'cell'", selectizeInput("kin_cell_choice", "Select Cell Type(s):", width = "100%", choices = cell_class_choices, multiple = TRUE)), div(class="d-flex justify-content-between align-items-end mt-4", tags$label("3. Paste Gene IDs (Entrez IDs):", class="fw-bold mb-0 text-primary"), div(actionButton("sample_kin_btn", "Load sample", class = "btn btn-sm btn-info text-white fw-bold"), actionButton("reset_kin_btn", "Reset", class = "btn btn-outline-secondary btn-sm ms-2"))), textAreaInput("kin_input", label = NULL, height = "150px", placeholder = "e.g. 5594\n5595\n4193\n836\n207\n3551"), uiOutput("kin_status_msg"), actionButton("run_kin_btn", "Map Internal Kinase Network →", class = "btn btn-lg btn-primary w-100 fw-bold mt-3"))), br(), uiOutput("kin_results_panel")
                     )
           )
  ),
  nav_panel(title = "Methodology", div(class="container mt-5", h3("Methodology", class="text-primary fw-bold"), p("Rigorous Over-Representation Analysis (ORA) statistical framework...")))
)

# ==========================================
# 5. SERVER LOGIC
# ==========================================
server <- function(input, output, session) {
  
  rv <- reactiveValues(home_page = "grid", enrich_page = "input", map_page = "input", exo_page = "input")
  
  raw_enrichment_data <- reactiveVal(NULL) 
  enrich_mapping_stats <- reactiveVal(list()) 
  mapping_data <- reactiveVal(NULL)
  exo_data <- reactiveVal(NULL)
  exo_summary_stats <- reactiveVal(list()) 
  tf_results <- reactiveVal(NULL)
  kin_results <- reactiveVal(NULL)
  selected_hm <- reactiveVal(NULL)
  
  output$home_page <- reactive({ rv$home_page })
  outputOptions(output, "home_page", suspendWhenHidden = FALSE)
  output$enrich_page <- reactive({ rv$enrich_page })
  outputOptions(output, "enrich_page", suspendWhenHidden = FALSE)
  output$map_page <- reactive({ rv$map_page })
  outputOptions(output, "map_page", suspendWhenHidden = FALSE)
  output$exo_page <- reactive({ rv$exo_page })
  outputOptions(output, "exo_page", suspendWhenHidden = FALSE)
  
  get_sym <- function(entrez_vec) {
    if(length(entrez_vec) == 0) return(character(0))
    syms <- entrez_to_symbol[as.character(entrez_vec)]
    syms[is.na(syms)] <- "Unknown"
    return(unname(syms))
  }
  
  get_sym_entrez <- function(entrez_vec) {
    if(length(entrez_vec) == 0) return(character(0))
    syms <- get_sym(entrez_vec)
    return(paste0(syms, " (", entrez_vec, ")"))
  }
  
  get_custom_icon_path <- function(hm_name) {
    safe_name <- gsub(" ", "_", tolower(hm_name))
    return(paste0("icons/", safe_name, ".png"))
  }
  
  output$hallmark_grid_ui <- renderUI({
    cols <- lapply(unique_hallmarks, function(hm) {
      btn_id <- paste0("btn_hm_", make.names(hm))
      img_path <- get_custom_icon_path(hm)
      div(class = "col-lg-3 col-md-4 col-sm-6 mb-4", actionButton(btn_id, HTML(paste0("<img src='", img_path, "' class='hallmark-custom-icon' alt=''><div class='hallmark-text'>", hm, "</div>")), class = "hallmark-box"))
    })
    div(class = "row justify-content-center", cols)
  })
  
  lapply(unique_hallmarks, function(hm) { btn_id <- paste0("btn_hm_", make.names(hm)); observeEvent(input[[btn_id]], { selected_hm(hm); rv$home_page <- "result" }) })
  observeEvent(input$back_home_btn, { rv$home_page <- "grid" })
  output$selected_hallmark_title <- renderText({ req(selected_hm()); paste("Genes associated with:", selected_hm()) })
  
  output$hallmark_genes_table <- renderDT({
    req(selected_hm())
    target_classes <- if(input$home_spec_level == "tissue") input$home_tissue_cell_choice else input$home_cell_choice
    filtered_db <- get_aggregated_database(input$home_db_choice, input$home_ev_choice, input$home_spec_level, target_classes)
    
    if(is.null(filtered_db)) { output$home_error_msg <- renderUI({ HTML("<div class='alert alert-warning mb-4'>No data available for these specific filter settings.</div>") }); return(NULL) }
    
    hm_genes <- filtered_db %>% dplyr::filter(HallMarks == selected_hm()) %>% dplyr::select(Gene_ID) %>% dplyr::distinct()
    if(nrow(hm_genes) == 0) { output$home_error_msg <- renderUI({ HTML("<div class='alert alert-warning mb-4'>No genes found for this specific hallmark under the current filter settings.</div>") }); return(NULL) }
    output$home_error_msg <- renderUI({ NULL })
    
    hm_genes <- hm_genes %>% 
      dplyr::mutate(
        Symbol = get_sym(Gene_ID),
        Symbol = paste0("<a href='https://www.ncbi.nlm.nih.gov/gene/", Gene_ID, "' target='_blank' style='text-decoration:none;'>", Symbol, "</a>"),
        Entrez_ID = paste0("<a href='https://www.ncbi.nlm.nih.gov/gene/", Gene_ID, "' target='_blank' style='text-decoration:none;'>", Gene_ID, "</a>")
      ) %>% 
      dplyr::select(Symbol, Entrez_ID)
    
    datatable(hm_genes, class = 'table table-striped table-hover', escape = FALSE, options = list(pageLength = 15), rownames = FALSE) %>% formatStyle('Symbol', color = '#0EA5E9', fontWeight = 'bold')
  })
  
  build_tissue_subselector <- function(tissue_input, output_id) {
    req(tissue_input)
    matched_classes_raw <- tissue_mapping_df %>% dplyr::filter(Tissue %in% tissue_input) %>% dplyr::pull(Cell_Type_Class) %>% unique() %>% sort()
    matched_classes_formatted <- gsub(" ", "_", tolower(matched_classes_raw))
    valid_classes <- intersect(matched_classes_formatted, cell_class_cols)
    if(length(valid_classes) > 0) { selectizeInput(output_id, "Included Cell Types (Click backspace to remove any):", choices = valid_classes, selected = valid_classes, multiple = TRUE, width = "100%") }
  }
  
  output$home_tissue_cell_ui <- renderUI({ build_tissue_subselector(input$home_tissue_choice, "home_tissue_cell_choice") })
  output$enrich_tissue_cell_ui <- renderUI({ build_tissue_subselector(input$enrich_tissue_choice, "enrich_tissue_cell_choice") })
  output$map_tissue_cell_ui <- renderUI({ build_tissue_subselector(input$map_tissue_choice, "map_tissue_cell_choice") })
  output$exo_tissue_cell_ui <- renderUI({ build_tissue_subselector(input$exo_tissue_choice, "exo_tissue_cell_choice") })
  output$tf_tissue_cell_ui <- renderUI({ build_tissue_subselector(input$tf_tissue_choice, "tf_tissue_cell_choice") })
  output$kin_tissue_cell_ui <- renderUI({ build_tissue_subselector(input$kin_tissue_choice, "kin_tissue_cell_choice") })
  
  get_aggregated_database <- function(db_choice, ev_choice, spec_level, target_classes) {
    filtered_db <- master_db %>% dplyr::filter(.data[[db_choice]] == "Yes", .data[[ev_choice]] == "Yes")
    if (spec_level == "all") { filtered_db <- filtered_db %>% dplyr::filter(.data[["Whole Body"]] == "Yes")
    } else {
      if (length(target_classes) == 0) return(NULL)
      filtered_db <- filtered_db %>% dplyr::filter(if_any(all_of(target_classes), ~ . == "Yes"))
    }
    if (nrow(filtered_db) == 0) return(NULL)
    return(filtered_db %>% dplyr::select(Gene_ID, HallMarks) %>% dplyr::distinct())
  }
  
  observeEvent(input$sample_btn, { updateTextAreaInput(session, "gene_input", value = "7007\n9370\n652\n6347\n6366\n940\n1437\n2920\n3576\n355") })
  observeEvent(input$reset_btn, { updateTextAreaInput(session, "gene_input", value = "") })
  observeEvent(input$back_enrich_btn, { rv$enrich_page <- "input" })
  
  observeEvent(input$run_btn, {
    raw_input <- unlist(strsplit(input$gene_input, "[\n, ]+")); input_genes <- unique(na.omit(as.character(raw_input))); input_genes <- input_genes[input_genes != ""] 
    valid_input_genes <- intersect(input_genes, universe_genes); invalid_input_genes <- setdiff(input_genes, universe_genes)
    if(length(valid_input_genes) == 0) { output$enrich_error_msg <- renderUI({ HTML("<div class='alert alert-danger'>No valid Entrez IDs found in background universe.</div>") }); return() }
    
    enrich_mapping_stats(list(total = length(input_genes), valid = length(valid_input_genes), invalid_genes = invalid_input_genes))
    target_classes <- if(input$enrich_spec_level == "tissue") input$enrich_tissue_cell_choice else input$enrich_cell_choice
    hallmarks_db <- get_aggregated_database(input$db_choice, input$enrich_ev_choice, input$enrich_spec_level, target_classes)
    
    if(is.null(hallmarks_db)) { output$enrich_error_msg <- renderUI({ HTML("<div class='alert alert-warning'>No data available for these specific selections.</div>") }); return() }
    
    results <- hallmarks_db %>% dplyr::distinct(HallMarks, Gene_ID) %>% dplyr::group_by(HallMarks) %>%
      dplyr::summarise(Hallmark_Size = dplyr::n(), Overlap_Count = sum(Gene_ID %in% valid_input_genes), Expected_Overlap = (Hallmark_Size / N_total_db_genes) * length(valid_input_genes), Enrichment_Fold = Overlap_Count / Expected_Overlap, Overlap_Symbols = paste(get_sym(unique(Gene_ID[Gene_ID %in% valid_input_genes])), collapse = ", "), Overlap_Entrez = paste(unique(Gene_ID[Gene_ID %in% valid_input_genes]), collapse = ", ")) %>%
      dplyr::filter(Overlap_Count > 0) %>%
      dplyr::mutate(P_Value = phyper(q = Overlap_Count - 1, m = Hallmark_Size, n = N_total_db_genes - Hallmark_Size, k = length(valid_input_genes), lower.tail = FALSE), FDR = p.adjust(P_Value, method = "BH")) %>% dplyr::arrange(FDR, desc(Overlap_Count))
    
    if(nrow(results) == 0) { output$enrich_error_msg <- renderUI({ HTML("<div class='alert alert-warning'>Zero Overlap.</div>") }); return() }
    output$enrich_error_msg <- renderUI({ NULL }); raw_enrichment_data(results); rv$enrich_page <- "result" 
  })
  
  output$mapping_summary_ui <- renderUI({ req(enrich_mapping_stats()); stats <- enrich_mapping_stats(); invalid_text <- if(length(stats$invalid_genes) > 0) paste(get_sym_entrez(stats$invalid_genes), collapse = ", ") else "None"; div(class = "summary-box shadow-sm mb-4", h5(class = "text-primary fw-bold mb-3", icon("info-circle"), " Input Mapping Summary"), tags$ul(class = "mb-0", tags$li(tags$strong("Total Input Genes: "), stats$total), tags$li(tags$strong("Mapped to Universe (Valid): "), stats$valid), tags$li(tags$strong("Unmapped Entrez IDs: "), tags$span(class="text-danger", invalid_text)))) })
  filtered_enrichment_data <- reactive({ req(raw_enrichment_data()); df <- raw_enrichment_data(); if (isTRUE(input$enable_filters)) df %>% dplyr::filter(FDR <= input$fdr_cutoff, Overlap_Count >= input$min_overlap) else df })
  output$kpi_boxes_enrichment <- renderUI({ req(filtered_enrichment_data()); df <- filtered_enrichment_data(); if(nrow(df) == 0) return(div(class="alert alert-warning", "No hallmarks pass cutoffs.")); unique_count <- length(unique(unlist(strsplit(df$Overlap_Entrez, ", ")))); layout_columns(value_box(title = "Enriched Hallmarks", value = nrow(df), showcase = icon("chart-pie"), theme = "success"), value_box(title = "Most Significant", value = tags$span(style="font-size: 1.1rem;", df$HallMarks[1]), showcase = icon("trophy"), theme = "primary"), value_box(title = "Unique Genes Mapped", value = unique_count, showcase = icon("dna"), theme = "info")) })
  output$enrichmentTable <- renderDT({ req(filtered_enrichment_data()); df <- filtered_enrichment_data(); if(nrow(df)==0) return(NULL); names(df) <- gsub("_", " ", names(df)); datatable(df, class = 'table table-striped', options = list(pageLength = 10, scrollX = TRUE), rownames = FALSE) %>% formatSignif(columns = c("P Value", "FDR", "Enrichment Fold"), digits = 3) %>% formatStyle('FDR', color = styleInterval(0.05, c('#10B981', '#EF4444')), fontWeight = 'bold') %>% formatStyle('HallMarks', color = '#0EA5E9', fontWeight = 'bold') })
  current_plot <- reactive({ req(filtered_enrichment_data()); df <- filtered_enrichment_data(); if(nrow(df) == 0) return(NULL); if(input$plot_type == "bar") { ggplot(df, aes(x = reorder(HallMarks, -log10(FDR)), y = -log10(FDR))) + geom_bar(stat = "identity", fill = "#0EA5E9", width = 0.7) + coord_flip() + theme_minimal(base_size = 14) + labs(x = "", y = "-log10(FDR)") + theme(axis.text.y = element_text(face = "bold")) } else if (input$plot_type == "bubble") { ggplot(df, aes(x = Overlap_Count, y = reorder(HallMarks, -log10(FDR)))) + geom_point(aes(size = -log10(FDR), color = Enrichment_Fold), alpha = 0.85) + scale_size_continuous(range = c(4, 15), name = "Significance\n(-log10 FDR)") + scale_color_gradient(low = "#3B82F6", high = "#EF4444", name = "Enrichment Fold") + labs(x = "Overlap Count (Number of Genes)", y = "") + theme_minimal(base_size = 14) + theme(axis.text.y = element_text(face = "bold", size = 12), axis.text.x = element_text(size = 12), axis.title.x = element_text(face = "bold", margin = margin(t = 15)), panel.grid.major.y = element_line(color = "#E2E8F0", linetype = "dashed"), panel.grid.major.x = element_line(color = "#E2E8F0"), panel.grid.minor = element_blank(), legend.position = "right", legend.box.background = element_rect(color = "#E2E8F0", fill = "white", linewidth = 0.5), legend.margin = margin(10, 10, 10, 10), legend.title = element_text(face = "bold", size = 11)) } else { ggplot(df, aes(x = reorder(gsub(" ", "\n", HallMarks), -log10(FDR)), y = -log10(FDR), fill = -log10(FDR))) + geom_bar(stat = "identity", alpha = 0.9) + geom_hline(yintercept = -log10(input$fdr_cutoff), color = "red", linetype = "dashed") + scale_fill_viridis_c(option = "turbo") + coord_polar(start = 0) + theme_minimal() + labs(x = NULL, y = NULL) + theme(axis.text.y = element_blank()) } })
  output$enrichmentPlot <- renderPlot({ current_plot() })
  output$dl_enrich_csv <- downloadHandler(filename = function() { paste("Hallmark_Enrichment_", Sys.Date(), ".csv", sep="") }, content = function(file) { req(filtered_enrichment_data()); write.csv(filtered_enrichment_data(), file, row.names = FALSE) })
  output$dl_enrich_plot <- downloadHandler(filename = function() { paste("Enrichment_Plot_", Sys.Date(), ".png", sep="") }, content = function(file) { req(current_plot()); ggsave(file, plot = current_plot(), width = 10, height = 8, dpi = 300, bg = "white") })
  
  observeEvent(input$sample_map_btn, { updateTextAreaInput(session, "mapping_input", value = "7007\n9370\n652\n6347\n6366") })
  observeEvent(input$reset_map_btn, { updateTextAreaInput(session, "mapping_input", value = "") })
  observeEvent(input$back_map_btn, { rv$map_page <- "input" })
  
  observeEvent(input$run_mapping_btn, {
    genes <- unique(na.omit(as.character(unlist(strsplit(input$mapping_input, "[\n, ]+"))))); genes <- genes[genes != ""]
    if(length(genes) == 0) { output$mapping_error_msg <- renderUI({ HTML("<div class='alert alert-danger'>Enter Entrez ID.</div>") }); return() }
    target_classes <- if(input$map_spec_level == "tissue") input$map_tissue_cell_choice else input$map_cell_choice
    target_db <- get_aggregated_database(input$mapping_db_choice, input$map_ev_choice, input$map_spec_level, target_classes)
    if(is.null(target_db)) { output$mapping_error_msg <- renderUI({ HTML("<div class='alert alert-warning'>No data available.</div>") }); return() }
    matches <- target_db %>% dplyr::filter(Gene_ID %in% genes) %>% dplyr::select(Gene_ID, HallMarks) %>% dplyr::distinct() %>% dplyr::arrange(Gene_ID)
    if(nrow(matches) == 0) { output$mapping_error_msg <- renderUI({ HTML("<div class='alert alert-warning'>Zero Overlap.</div>") }); return() }
    matches$Symbol <- get_sym(matches$Gene_ID)
    matches <- matches %>% dplyr::select(Symbol, Entrez_ID = Gene_ID, HallMarks)
    output$mapping_error_msg <- renderUI({ NULL }); mapping_data(matches); rv$map_page <- "result" 
  })
  output$kpi_boxes_mapping <- renderUI({ req(mapping_data()); df <- mapping_data(); layout_columns(value_box(title = "Genes Mapped", value = length(unique(df$Entrez_ID)), showcase = icon("check-circle"), theme = "success"), value_box(title = "Unique Hallmarks", value = length(unique(df$HallMarks)), showcase = icon("layer-group"), theme = "primary")) })
  output$mappingTable <- renderDT({ req(mapping_data()); df <- mapping_data(); names(df) <- gsub("_", " ", names(df)); datatable(df, class = 'table table-striped', options = list(pageLength = 15), rownames = FALSE) %>% formatStyle('HallMarks', color = '#0EA5E9', fontWeight = 'bold') })
  output$dl_map_csv <- downloadHandler(filename = function() { paste("Gene_Hallmark_Mapping_", Sys.Date(), ".csv", sep="") }, content = function(file) { req(mapping_data()); write.csv(mapping_data(), file, row.names = FALSE) })
  
  observeEvent(input$sample_exo_btn, { updateTextAreaInput(session, "exo_input", value = "7007\n9370\n652\n6347\n6366\n9999999") })
  observeEvent(input$reset_exo_btn, { updateTextAreaInput(session, "exo_input", value = "") })
  observeEvent(input$back_exo_btn, { rv$exo_page <- "input" })
  
  observeEvent(input$run_exo_btn, {
    genes <- unique(na.omit(as.character(unlist(strsplit(input$exo_input, "[\n, ]+"))))); genes <- genes[genes != ""]
    if(length(genes) == 0) { output$exo_error_msg <- renderUI({ HTML("<div class='alert alert-danger'>Enter Entrez IDs.</div>") }); return() }
    target_classes <- if(input$exo_spec_level == "tissue") input$exo_tissue_cell_choice else input$exo_cell_choice
    target_hallmark_db <- get_aggregated_database(input$exo_db_choice, input$exo_ev_choice, input$exo_spec_level, target_classes)
    if(is.null(target_hallmark_db) || nrow(target_hallmark_db) == 0) { output$exo_error_msg <- renderUI({ HTML("<div class='alert alert-warning'>No data available for the selected hallmark context.</div>") }); return() }
    dynamic_hallmark_genes <- unique(na.omit(as.character(target_hallmark_db$Gene_ID)))
    hallmark_filtered_genes <- intersect(genes, dynamic_hallmark_genes)
    invalid_genes <- setdiff(genes, dynamic_hallmark_genes)
    if(length(hallmark_filtered_genes) == 0) { output$exo_error_msg <- renderUI({ HTML("<div class='alert alert-warning'>None of the inputted genes exist in your selected background aging hallmarks database context.</div>") }); return() }
    exo_summary_stats(list(total_input = length(genes), valid = length(hallmark_filtered_genes), invalid_genes = invalid_genes))
    results_df <- data.frame(Symbol = get_sym(hallmark_filtered_genes), Entrez_ID = hallmark_filtered_genes, stringsAsFactors = FALSE)
    results_df$In_Exosome <- ifelse(hallmark_filtered_genes %in% combined_exo_genes, "Yes", "No")
    output$exo_error_msg <- renderUI({ NULL }); exo_data(results_df); rv$exo_page <- "result" 
  })
  output$exo_summary_ui <- renderUI({ req(exo_summary_stats()); stats <- exo_summary_stats(); invalid_text <- if(length(stats$invalid_genes) > 0) paste(get_sym_entrez(stats$invalid_genes), collapse = ", ") else "None"; div(class = "summary-box shadow-sm mb-4", h5(class = "text-primary fw-bold mb-3", icon("info-circle"), " Input Processing Summary"), tags$ul(class = "mb-0", tags$li(tags$strong("Total Inputted Genes: "), stats$total_input), tags$li(tags$strong("Valid (Present in Selected DB): "), stats$valid), tags$li(tags$strong("Discarded (Not in Selected DB): "), tags$span(class="text-danger", invalid_text)))) })
  output$exoTable <- renderDT({ req(exo_data()); df <- exo_data(); names(df) <- gsub("_", " ", names(df)); datatable(df, class = 'table table-striped', options = list(pageLength = 15), rownames = FALSE) %>% formatStyle('Symbol', color = '#0EA5E9', fontWeight = 'bold') %>% formatStyle('In Exosome', color = styleEqual(c("Yes", "No"), c('#10B981', '#EF4444')), fontWeight = 'bold') })
  output$dl_exo_csv <- downloadHandler(filename = function() { paste("Exosome_Presence_Check_", Sys.Date(), ".csv", sep="") }, content = function(file) { req(exo_data()); write.csv(exo_data(), file, row.names = FALSE) })
  
  process_internal_network <- function(raw_input, db_network, target_db, network_type) {
    raw_entrez <- unique(na.omit(trimws(unlist(strsplit(raw_input, "[\n, ]+")))))
    raw_entrez <- raw_entrez[raw_entrez != ""]
    if(length(raw_entrez) == 0) return(list(error = "Please enter valid Entrez IDs."))
    syms <- entrez_to_symbol[raw_entrez]
    mapping_df <- data.frame(Entrez = raw_entrez[!is.na(syms)], Symbol = syms[!is.na(syms)], stringsAsFactors = FALSE)
    input_genes <- unique(mapping_df$Symbol)
    if(length(input_genes) == 0) return(list(error = "Could not convert Entrez IDs to Gene Symbols."))
    if(nrow(db_network) == 0) return(list(error = paste("The", network_type, "network file is missing or empty.")))
    
    context_entrez <- unique(as.character(na.omit(target_db$Gene_ID)))
    context_symbols <- unique(na.omit(entrez_to_symbol[context_entrez]))
    valid_input_genes <- intersect(input_genes, context_symbols)
    
    bio_edges <- db_network %>% dplyr::filter(source %in% valid_input_genes & target %in% valid_input_genes) %>% drop_na(source, target)
    invalid_entrez <- setdiff(raw_entrez, context_entrez)
    invalid_context_genes_ui <- get_sym_entrez(invalid_entrez)
    valid_no_edge_symbols <- setdiff(valid_input_genes, unique(c(bio_edges$source, bio_edges$target)))
    valid_no_edge_entrez <- mapping_df$Entrez[mapping_df$Symbol %in% valid_no_edge_symbols]
    valid_no_edge_genes_ui <- get_sym_entrez(valid_no_edge_entrez)
    if(nrow(bio_edges) == 0) return(list(error = NULL, edges = data.frame(), table = data.frame(), bio_edges = data.frame(), hm_mapping = data.frame(), net_entrez_map = c(), invalid_context_genes = invalid_context_genes_ui, valid_no_edge_genes = valid_no_edge_genes_ui))
    
    net_genes <- unique(c(bio_edges$source, bio_edges$target))
    net_entrez_map <- symbol_to_entrez[net_genes]
    net_entrez_chars <- as.character(net_entrez_map)
    hm_mapping <- target_db %>% dplyr::filter(Gene_ID %in% net_entrez_chars) %>% dplyr::select(Gene_ID, HallMarks) %>% dplyr::distinct()
    hm_edges <- data.frame(source = character(), target = character(), type = character(), stringsAsFactors=FALSE)
    if(nrow(hm_mapping) > 0) {
      hm_mapping$Symbol <- names(net_entrez_map)[match(hm_mapping$Gene_ID, net_entrez_map)]
      hm_edges <- data.frame(source = hm_mapping$Symbol, target = hm_mapping$HallMarks, type = "hallmark_association", stringsAsFactors=FALSE) %>% drop_na(source, target)
    }
    bio_edges$type <- "regulation"
    all_edges <- bind_rows(bio_edges, hm_edges)
    regulators <- unique(bio_edges$source)
    table_data <- data.frame(Regulator_Symbol = regulators, Regulator_Entrez = sapply(regulators, function(s) { net_entrez_map[s] }), Target_Symbol = "", Target_Entrez = "", Target_Count = 0, stringsAsFactors = FALSE)
    if(nrow(table_data) > 0) {
      for(i in 1:nrow(table_data)) {
        reg <- table_data$Regulator_Symbol[i]; targs <- bio_edges %>% dplyr::filter(source == reg) %>% dplyr::pull(target)
        table_data$Target_Symbol[i] <- paste(targs, collapse = ", "); table_data$Target_Entrez[i] <- paste(sapply(targs, function(t) net_entrez_map[t]), collapse = ", "); table_data$Target_Count[i] <- length(targs)
      }
      table_data <- table_data %>% dplyr::arrange(desc(Target_Count))
    }
    return(list(error = NULL, edges = all_edges, table = table_data, bio_edges = bio_edges, hm_mapping = hm_mapping, net_entrez_map = net_entrez_map, input_genes = input_genes, invalid_context_genes = invalid_context_genes_ui, valid_no_edge_genes = valid_no_edge_genes_ui))
  }
  
  observeEvent(input$sample_tf_btn, { updateTextAreaInput(session, "tf_input", value = "7124\n2\n1029\n5970\n2099\n4609\n1958\n3725") })
  observeEvent(input$reset_tf_btn, { updateTextAreaInput(session, "tf_input", value = "") })
  observeEvent(input$run_tf_btn, {
    output$tf_status_msg <- renderUI({ HTML("<div class='alert alert-info mt-3'>Loading database and mapping...</div>") })
    target_classes <- if(input$tf_spec_level == "tissue") input$tf_tissue_cell_choice else input$tf_cell_choice
    target_db <- get_aggregated_database(input$tf_db_choice, input$tf_ev_choice, input$tf_spec_level, target_classes)
    if(is.null(target_db)) { output$tf_status_msg <- renderUI({ HTML("<div class='alert alert-warning mt-3'>No data available.</div>") }); return() }
    res <- process_internal_network(input$tf_input, db_tf, target_db, "Transcription Factors")
    if(!is.null(res$error)) { output$tf_status_msg <- renderUI({ HTML(paste0("<div class='alert alert-danger mt-3'>", res$error, "</div>")) }); return() }
    output$tf_status_msg <- renderUI({ NULL }); tf_results(res)
  })
  output$dl_tf_csv <- downloadHandler(filename = function() { paste("TF_Network_", Sys.Date(), ".csv", sep="") }, content = function(file) { req(tf_results()); write.csv(tf_results()$table, file, row.names = FALSE) })
  
  output$tf_results_panel <- renderUI({
    req(tf_results()); res <- tf_results(); missing_ui <- tagList()
    num_missing_db <- length(res$invalid_context_genes); num_missing_edge <- length(res$valid_no_edge_genes); total_missing <- num_missing_db + num_missing_edge
    if (total_missing > 0) {
      missing_db_text <- if(num_missing_db > 0) paste0("<strong class='text-danger'>Not in Background DB (", num_missing_db, "):</strong> ", paste(res$invalid_context_genes, collapse=", "), "<br><br>") else ""
      missing_edge_text <- if(num_missing_edge > 0) paste0("<strong class='text-warning-emphasis'>No Internal Edges (", num_missing_edge, "):</strong> ", paste(res$valid_no_edge_genes, collapse=", ")) else ""
      missing_ui <- HTML(paste0("<details style='background-color: var(--bs-warning-bg-subtle); border: 1px solid var(--bs-warning-border-subtle); border-radius: 6px; padding: 10px 15px; margin-bottom: 20px;'><summary style='font-weight: bold; cursor: pointer; color: var(--bs-warning-text-emphasis); outline: none;'>⚠️ ", total_missing, " genes excluded from network (Click to view details)</summary><div style='margin-top: 12px; font-size: 0.9rem; max-height: 120px; overflow-y: auto; padding-right: 10px;'>", missing_db_text, missing_edge_text, "</div></details>"))
    }
    if (nrow(res$edges) == 0) { return(div(missing_ui, br(), div(class="alert alert-warning", "No internal edges exist between the supplied genes."))) }
    div(missing_ui, card(class = "shadow-sm border-0 mb-4", card_header(class = "bg-transparent text-primary fw-bold", "Strict Internal TF Regulatory Map"), card_body(HTML("<div style='background: var(--bs-tertiary-bg); border: 1px solid var(--bs-border-color); border-radius: 6px; padding: 12px 20px; margin-bottom: 15px; display: flex; gap: 20px; flex-wrap: wrap; font-size: 0.95rem; align-items: center;'><strong style='color: var(--bs-body-color);'>Legend:</strong><span><span style='color: #FBBF24; font-size: 1.3rem; line-height: 1;'>★</span> Aging Hallmark</span><span><span style='color: #EF4444; font-size: 1.2rem; line-height: 1;'>▲</span> Transcription Factor</span><span><span style='color: #9333EA; font-size: 1.4rem; line-height: 1;'>■</span> Target & TF</span><span><span style='color: #0EA5E9; font-size: 1.2rem; line-height: 1;'>●</span> Target Gene</span><span style='margin-left: auto;'><span style='color: #CBD5E1; font-weight: bold; letter-spacing: 2px;'>---</span> Hallmark Association</span><span><span style='color: #10B981; font-weight: bold;'>→</span> Stimulation</span><span><span style='color: #EF4444; font-weight: bold; font-size: 1.2rem;'>—</span> Inhibition</span><span><span style='color: #9CA3AF; font-weight: bold;'>→</span> Unknown</span></div>"), div(class = "plot-container", visNetworkOutput("tfNetwork", height = "750px")))), card(class = "shadow-sm border-0", card_header(class = "bg-transparent d-flex justify-content-between align-items-center", span(class="h4 text-primary fw-bold m-0", "TF-Target Breakdown"), downloadButton("dl_tf_csv", "📥 Download Table", class = "btn btn-sm btn-success fw-bold")), card_body(DTOutput("tfTable"))))
  })
  
  output$tfNetwork <- renderVisNetwork({
    req(tf_results()); res <- tf_results(); edges <- res$edges; hm_map <- res$hm_mapping
    node_ids <- unique(na.omit(c(edges$source, edges$target))); nodes <- data.frame(id = node_ids, stringsAsFactors = FALSE); nodes$label <- nodes$id
    is_hm <- nodes$id %in% hm_map$HallMarks; is_src <- nodes$id %in% res$bio_edges$source; is_tgt <- nodes$id %in% res$bio_edges$target
    nodes$group <- dplyr::case_when(is_hm ~ "Aging Hallmark", is_src & is_tgt ~ "TF & Target", is_src ~ "Transcription Factor", TRUE ~ "Target Gene")
    nodes$shape <- dplyr::case_when(nodes$group == "Aging Hallmark" ~ "star", nodes$group == "Transcription Factor" ~ "triangle", nodes$group == "TF & Target" ~ "box", TRUE ~ "dot")
    nodes$size <- ifelse(nodes$group == "Aging Hallmark", 45, ifelse(nodes$group == "Target Gene", 20, 35))
    nodes$title <- sapply(nodes$id, function(n) { if(n %in% hm_map$HallMarks) return(paste0("<b>Hallmark:</b> ", n)); ent <- res$net_entrez_map[n]; paste0("<b>Gene:</b> ", n, "<br><b>Entrez ID:</b> ", ent) })
    
    vis_edges_base <- data.frame(from = edges$source, to = edges$target, type = edges$type, is_stimulation = edges$is_stimulation, is_inhibition = edges$is_inhibition, stringsAsFactors = FALSE) %>% drop_na(from, to)
    vis_edges_base$stim <- grepl("true|1", as.character(vis_edges_base$is_stimulation), ignore.case=TRUE)
    vis_edges_base$inhib <- grepl("true|1", as.character(vis_edges_base$is_inhibition), ignore.case=TRUE)
    vis_edges_base$stim[is.na(vis_edges_base$stim)] <- FALSE
    vis_edges_base$inhib[is.na(vis_edges_base$inhib)] <- FALSE
    
    both_idx <- vis_edges_base$stim & vis_edges_base$inhib
    df_normal <- vis_edges_base[!both_idx, ]
    df_both_stim <- vis_edges_base[both_idx, ]; df_both_stim$inhib <- FALSE
    df_both_inhib <- vis_edges_base[both_idx, ]; df_both_inhib$stim <- FALSE
    
    final_edges <- bind_rows(df_normal, df_both_stim, df_both_inhib)
    final_edges$dashes <- (final_edges$type == "hallmark_association")
    final_edges$color <- dplyr::case_when(final_edges$type == "hallmark_association" ~ "#CBD5E1", final_edges$stim & !final_edges$inhib ~ "#10B981", !final_edges$stim & final_edges$inhib ~ "#EF4444", TRUE ~ "#9CA3AF")
    final_edges$arrows <- dplyr::case_when(final_edges$type == "hallmark_association" ~ "", !final_edges$stim & final_edges$inhib ~ "", TRUE ~ "to")
    final_edges$smooth <- TRUE
    
    visNetwork(nodes, final_edges) %>% visGroups(groupname = "Aging Hallmark", color = list(background = "#FBBF24", border = "#D97706")) %>% visGroups(groupname = "Transcription Factor", color = list(background = "#EF4444", border = "#B91C1C")) %>% visGroups(groupname = "TF & Target", color = list(background = "#9333EA", border = "#6D28D9")) %>% visGroups(groupname = "Target Gene", color = list(background = "#0EA5E9", border = "#0284C7")) %>% visExport(type = "png", name = "TF_Network", label = "Export PNG", style = "background-color: #0EA5E9; color: white; border: none; padding: 5px 10px; border-radius: 4px; cursor: pointer; float: right; margin-top: 10px;") %>% visOptions(highlightNearest = list(enabled = TRUE, degree = 1, hover = TRUE), nodesIdSelection = TRUE) %>% visEdges(smooth = list(enabled = TRUE, type = "curvedCW")) %>% visPhysics(solver = "forceAtlas2Based", forceAtlas2Based = list(gravitationalConstant = -100, springLength = 200), stabilization = list(enabled = TRUE, iterations = 150)) %>% visEvents(type = "once", stabilizationIterationsDone = "function() { this.setOptions({ physics: false }); }") %>% visInteraction(dragNodes = TRUE)
  })
  output$tfTable <- renderDT({ req(tf_results()); df <- tf_results()$table; names(df) <- gsub("_", " ", names(df)); datatable(df, class = 'table table-striped table-hover', options = list(pageLength = 10, scrollX = TRUE), rownames = FALSE) %>% formatStyle('Regulator Symbol', color = '#EF4444', fontWeight = 'bold') })
  
  observeEvent(input$sample_kin_btn, { updateTextAreaInput(session, "kin_input", value = "5594\n5595\n4193\n836\n207\n3551") })
  observeEvent(input$reset_kin_btn, { updateTextAreaInput(session, "kin_input", value = "") })
  observeEvent(input$run_kin_btn, {
    output$kin_status_msg <- renderUI({ HTML("<div class='alert alert-info mt-3'>Loading database and mapping...</div>") })
    target_classes <- if(input$kin_spec_level == "tissue") input$kin_tissue_cell_choice else input$kin_cell_choice
    target_db <- get_aggregated_database(input$kin_db_choice, input$kin_ev_choice, input$kin_spec_level, target_classes)
    if(is.null(target_db)) { output$kin_status_msg <- renderUI({ HTML("<div class='alert alert-warning mt-3'>No data available.</div>") }); return() }
    res <- process_internal_network(input$kin_input, db_kinase, target_db, "Kinases")
    if(!is.null(res$error)) { output$kin_status_msg <- renderUI({ HTML(paste0("<div class='alert alert-danger mt-3'>", res$error, "</div>")) }); return() }
    output$kin_status_msg <- renderUI({ NULL }); kin_results(res)
  })
  output$dl_kin_csv <- downloadHandler(filename = function() { paste("Kinase_Network_", Sys.Date(), ".csv", sep="") }, content = function(file) { req(kin_results()); write.csv(kin_results()$table, file, row.names = FALSE) })
  
  output$kin_results_panel <- renderUI({
    req(kin_results()); res <- kin_results(); missing_ui <- tagList()
    num_missing_db <- length(res$invalid_context_genes); num_missing_edge <- length(res$valid_no_edge_genes); total_missing <- num_missing_db + num_missing_edge
    if (total_missing > 0) {
      missing_db_text <- if(num_missing_db > 0) paste0("<strong class='text-danger'>Not in Background DB (", num_missing_db, "):</strong> ", paste(res$invalid_context_genes, collapse=", "), "<br><br>") else ""
      missing_edge_text <- if(num_missing_edge > 0) paste0("<strong class='text-warning-emphasis'>No Internal Edges (", num_missing_edge, "):</strong> ", paste(res$valid_no_edge_genes, collapse=", ")) else ""
      missing_ui <- HTML(paste0("<details style='background-color: var(--bs-warning-bg-subtle); border: 1px solid var(--bs-warning-border-subtle); border-radius: 6px; padding: 10px 15px; margin-bottom: 20px;'><summary style='font-weight: bold; cursor: pointer; color: var(--bs-warning-text-emphasis); outline: none;'>⚠️ ", total_missing, " genes excluded from network (Click to view details)</summary><div style='margin-top: 12px; font-size: 0.9rem; max-height: 120px; overflow-y: auto; padding-right: 10px;'>", missing_db_text, missing_edge_text, "</div></details>"))
    }
    if (nrow(res$edges) == 0) { return(div(missing_ui, br(), div(class="alert alert-warning", "No internal edges exist between the supplied genes."))) }
    div(missing_ui, card(class = "shadow-sm border-0 mb-4", card_header(class = "bg-transparent text-primary fw-bold", "Internal Kinase Phosphorylation Map"), card_body(HTML("<div style='background: var(--bs-tertiary-bg); border: 1px solid var(--bs-border-color); border-radius: 6px; padding: 12px 20px; margin-bottom: 15px; display: flex; gap: 20px; flex-wrap: wrap; font-size: 0.95rem; align-items: center;'><strong style='color: var(--bs-body-color);'>Legend:</strong><span><span style='color: #FBBF24; font-size: 1.3rem; line-height: 1;'>★</span> Aging Hallmark</span><span><span style='color: #8B5CF6; font-size: 1.4rem; line-height: 1;'>◆</span> Kinase</span><span><span style='color: #EC4899; font-size: 1.4rem; line-height: 1;'>⬢</span> Substrate & Kinase</span><span><span style='color: #10B981; font-size: 1.2rem; line-height: 1;'>●</span> Substrate</span><span style='margin-left: auto;'><span style='color: #CBD5E1; font-weight: bold; letter-spacing: 2px;'>---</span> Hallmark Association</span><span><span style='color: #9CA3AF; font-weight: bold;'>→</span> Phosphorylation</span></div>"), div(class = "plot-container", visNetworkOutput("kinNetwork", height = "750px")))), card(class = "shadow-sm border-0", card_header(class = "bg-transparent d-flex justify-content-between align-items-center", span(class="h4 text-primary fw-bold m-0", "Kinase-Substrate Breakdown"), downloadButton("dl_kin_csv", "📥 Download Table", class = "btn btn-sm btn-success fw-bold")), card_body(DTOutput("kinTable"))))
  })
  
  output$kinNetwork <- renderVisNetwork({
    req(kin_results()); res <- kin_results(); edges <- res$edges; hm_map <- res$hm_mapping
    node_ids <- unique(na.omit(c(edges$source, edges$target))); nodes <- data.frame(id = node_ids, stringsAsFactors = FALSE); nodes$label <- nodes$id
    is_hm <- nodes$id %in% hm_map$HallMarks; is_src <- nodes$id %in% res$bio_edges$source; is_tgt <- nodes$id %in% res$bio_edges$target
    nodes$group <- dplyr::case_when(is_hm ~ "Aging Hallmark", is_src & is_tgt ~ "Kinase & Substrate", is_src ~ "Kinase", TRUE ~ "Substrate")
    nodes$shape <- dplyr::case_when(nodes$group == "Aging Hallmark" ~ "star", nodes$group == "Kinase" ~ "diamond", nodes$group == "Kinase & Substrate" ~ "hexagon", TRUE ~ "dot")
    nodes$size <- ifelse(nodes$group == "Aging Hallmark", 45, ifelse(nodes$group == "Substrate", 20, 35))
    nodes$title <- sapply(nodes$id, function(n) { if(n %in% hm_map$HallMarks) return(paste0("<b>Hallmark:</b> ", n)); ent <- res$net_entrez_map[n]; paste0("<b>Gene:</b> ", n, "<br><b>Entrez ID:</b> ", ent) })
    
    vis_edges <- data.frame(from = edges$source, to = edges$target, dashes = (edges$type == "hallmark_association"), stringsAsFactors = FALSE) %>% drop_na(from, to)
    vis_edges$color <- ifelse(vis_edges$dashes, "#CBD5E1", "#9CA3AF"); vis_edges$arrows <- ifelse(vis_edges$dashes, "", "to")
    
    visNetwork(nodes, vis_edges) %>% visGroups(groupname = "Aging Hallmark", color = list(background = "#FBBF24", border = "#D97706")) %>% visGroups(groupname = "Kinase", color = list(background = "#8B5CF6", border = "#6D28D9")) %>% visGroups(groupname = "Kinase & Substrate", color = list(background = "#EC4899", border = "#A855F7")) %>% visGroups(groupname = "Substrate", color = list(background = "#10B981", border = "#047857")) %>% visExport(type = "png", name = "Kinase_Network", label = "Export PNG", style = "background-color: #0EA5E9; color: white; border: none; padding: 5px 10px; border-radius: 4px; cursor: pointer; float: right; margin-top: 10px;") %>% visOptions(highlightNearest = list(enabled = TRUE, degree = 1, hover = TRUE), nodesIdSelection = TRUE) %>% visEdges(smooth = FALSE) %>% visPhysics(solver = "forceAtlas2Based", forceAtlas2Based = list(gravitationalConstant = -100, springLength = 200), stabilization = list(enabled = TRUE, iterations = 150)) %>% visEvents(type = "once", stabilizationIterationsDone = "function() { this.setOptions({ physics: false }); }") %>% visInteraction(dragNodes = TRUE)
  })
  output$kinTable <- renderDT({ req(kin_results()); df <- kin_results()$table; names(df) <- c("Kinase Symbol", "Kinase Entrez", "Substrate Symbol", "Substrate Entrez", "Target Count"); datatable(df, class = 'table table-striped table-hover', options = list(pageLength = 10, scrollX = TRUE), rownames = FALSE) %>% formatStyle('Kinase Symbol', color = '#8B5CF6', fontWeight = 'bold') })
  
}

shinyApp(ui = ui, server = server)