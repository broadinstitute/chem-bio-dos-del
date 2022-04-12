# Import packages ####
packages <- c("data.table", "tidyverse", "plotly", "shiny", "shinyjs", "shinythemes", "DT", "jsonlite", "rcdk", "rclipboard")
package_check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = T)) {
      install.packages(x, dependencies = T)
      library(x, character.only = T)
    }
  }
)

# Read config file ####
config <- read_json("config_file.json")

# Import functions ####
source(file.path(config$project_dir, config$functions_script))

# Import metadata ####
lib_meta <- file.path(config$project_dir, config$library_metadata) %>% fread() %>% as_tibble()
comb_lib_meta <- file.path(config$project_dir, config$combined_library_metadata) %>% fread() %>% as_tibble()
bb_meta <- file.path(config$project_dir, config$building_block_metadata) %>% fread() %>% as_tibble()
run_meta <- file.path(config$project_dir, config$sequencing_run_metadata) %>% fread() %>% as_tibble()
ct_file_meta <- file.path(config$project_dir, config$count_filename_metadata) %>% fread() %>% as_tibble()
samp_meta <- file.path(config$project_dir, config$sample_metadata) %>% fread() %>% as_tibble()
targ_meta <- file.path(config$project_dir, config$target_metadata) %>% fread() %>% as_tibble()
resin_meta <- file.path(config$project_dir, config$resin_metadata) %>% fread() %>% as_tibble()
scr_meta <- file.path(config$project_dir, config$screen_metadata) %>% fread() %>% as_tibble()
agg_meta <- file.path(config$project_dir, config$aggregation_metadata) %>% fread() %>% as_tibble()
merged_meta <- ct_file_meta %>%
  left_join(run_meta) %>%
  left_join(samp_meta) %>%
  left_join(scr_meta) %>%
  left_join(targ_meta) %>%
  left_join(resin_meta)
res_lst <- list.files(file.path(config$project_dir, config$result_metadata_location), recursive = F, full.names = F) %>% str_remove_all(".csv")
if(length(res_lst) == 0) { res_lst <- NULL }
agg_lst <- as_tibble(fread(file.path(config$project_dir, config$aggregation_metadata))) %>%
  mutate(temp = paste0(agg_id, ": ", an_type_name)) %>%
  pull(temp)
lib_lst <- lib_meta %>%
  mutate(temp = paste0(lib_id, ": ", lib_name)) %>%
  pull(temp)

# Run Shiny app ####
for(i in 1) {
  ui <- fluidPage(
    useShinyjs(),
    theme = shinytheme("simplex"),
    titlePanel("DEL Analysis App 1.0"),
    sidebarLayout(
      sidebarPanel(
        conditionalPanel(
          'input.dataset === "lib_meta"',
          div("Metadata associated with the DNA-encoded libraries."),
          br(),
          div("Submit this form to add new entries:"),
          textInput("user_lib_seq", "lib_seq", ""),
          textInput("user_lib_name", "lib_name", ""),
          textInput("user_chemist", "chemist", ""),
          textInput("user_lib_completion_date", "lib_completion_date", ""),
          textInput("user_lib_loc", "lib_loc", ""),
          actionButton("append_lib_meta", "Add entry"),
          actionButton("save_lib_meta", "Save changes"),
          actionButton("reload_lib_meta", "Reload")
        ),
        conditionalPanel(
          'input.dataset === "comb_lib_meta"',
          div("Metadata associated with the pooled DNA-encoded libraries. Each unique combination of pools, including the ratios, is assigned a unique ID."),
          br(),
          div("Submit this form to add new entries (double-click to edit):"),
          DT::dataTableOutput("input_comb_lib"),
          actionButton("input_comb_lib_add_row", "Add row"),
          actionButton("input_comb_lib_remove_row", "Remove row"),
          br(),
          br(),
          actionButton("append_comb_lib_meta", "Add entry"),
          actionButton("save_comb_lib_meta", "Save changes"),
          actionButton("reload_comb_lib_meta", "Reload")
        ),
        conditionalPanel(
          'input.dataset === "targ_meta"',
          div("Metadata associated with the targets."),
          br(),
          div("Submit this form to add new entries:"),
          textInput("user_targ_name", "targ_name", ""),
          actionButton("append_targ_meta", "Add entry"),
          actionButton("save_targ_meta", "Save changes"),
          actionButton("reload_targ_meta", "Reload")
        ),
        conditionalPanel(
          'input.dataset === "resin_meta"',
          div("Metadata associated with the resins (i.e., beads)."),
          br(),
          div("Submit this form to add new entries:"),
          textInput("user_resin_name", "resin_name", ""),
          actionButton("append_resin_meta", "Add entry"),
          actionButton("save_resin_meta", "Save changes"),
          actionButton("reload_resin_meta", "Reload")
        ),
        conditionalPanel(
          'input.dataset === "scr_meta"',
          div("Metadata associated with the DEL screens."),
          br(),
          div("Submit this form to add new entries:"),
          textInput("user_scr_name", "scr_name", ""),
          textInput("user_scr_completion_date", "scr_completion_date", ""),
          actionButton("append_scr_meta", "Add entry"),
          actionButton("save_scr_meta", "Save changes"),
          actionButton("reload_scr_meta", "Reload")
        ),
        conditionalPanel(
          'input.dataset === "samp_meta"',
          div("Metadata associated with the samples. To add new rows, use the 'New samples' tab."),
          br(),
          actionButton("save_samp_meta", "Save changes"),
          actionButton("reload_samp_meta", "Reload")
        ),
        conditionalPanel(
          'input.dataset === "merged_meta"',
          div("Merged metadata to aid in defining analysis experiments. Click on rows to select them, then click 'Save selection' to export them as a csv file. To use this output for analyses, add another column named 'description,' which will be used to axis labels."),
          br(),
          checkboxGroupInput("show_vars", "Columns in merged_meta to show:",
                             names(merged_meta), selected = names(merged_meta)),
          fluidRow(
            actionButton("merged_select_all_cols", "Select all columns"),
            actionButton("merged_clear_all_cols", "Clear column selection")
          ),
          br(),
          actionButton("reload_merged_meta", "Reload")
        ),
        conditionalPanel(
          'input.dataset === "New samples"',
          div("Use an input csv table to update the samp_meta table."),
          br(),
          fileInput("new_samples", "Choose csv file",
                    multiple = FALSE,
                    accept = c("text/csv",
                               "text/comma-separated-values,text/plain",
                               ".csv")),
          actionButton("append_samp_meta", "Add to samp_meta")
        ),
        conditionalPanel(
          'input.dataset === "run_meta"',
          div("Metadata associated with the sequencing runs."),
          br(),
          div("Submit this form to add new entries:"),
          textInput("user_seq_type", "seq_type", ""),
          textInput("user_seq_location", "seq_location", ""),
          textInput("user_seq_completion_date", "seq_completion_date", ""),
          textInput("user_fastq_folder_name", "fastq_folder_name", ""),
          actionButton("append_run_meta", "Add entry"),
          actionButton("save_run_meta", "Save changes"),
          actionButton("reload_run_meta", "Reload"),
        ),
        conditionalPanel(
          'input.dataset === "Create config"',
          div("Create config file (json) from a csv table."),
          br(),
          div("Note: Before the new counts can be used for analyses, ct_file_meta must be updated (click below). The same csv file can be used for both tasks."),
          br(),
          fileInput("config_upload", "Choose csv file",
                    multiple = FALSE,
                    accept = c("text/csv",
                               "text/comma-separated-values,text/plain",
                               ".csv")),
          actionButton("update_ct_file_meta", "Update ct_file_meta"),
          br(),
          br(),
          downloadButton("download_json", "Create json")
        ),
        conditionalPanel(
          'input.dataset === "Get structure"',
          div("Get expected structures (SMILES) of library members. Metadata must be loaded (once per session) before structures can be generated."),
          actionButton("load_cpd_struct", "Load metadata"),
          br(),
          br(),
          div("Single compound:"),
          fluidRow(column(6, numericInput("input_cpd_struct_single_lib_id", "lib_id", NA)),
                   column(6, numericInput("input_cpd_struct_single_cy1", "cycle1", NA)),
                   column(6, numericInput("input_cpd_struct_single_cy2", "cycle2", NA)),
                   column(6, numericInput("input_cpd_struct_single_cy3", "cycle3", NA))),
          actionButton("get_cpd_struct", "Get structure"),
          br(),
          br(),
          fileInput("cpd_struct_upload", "Multiple compounds:",
                    multiple = FALSE,
                    accept = c("text/csv",
                               "text/comma-separated-values,text/plain",
                               ".csv")),
          actionButton("get_cpd_structs", "Get structures"),
          actionButton("clear_cpd_struct", "Clear"),
          downloadButton("save_cpd_struct", "Save")
        ),
        conditionalPanel(
          'input.dataset === "Run analysis (simple entry)"',
          div("Perform analysis, provided a list of controls and targets. Uploaded files must contain columns 'ct_filename' and 'description'."),
          br(),
          div("Note: Samples that are replicates should be given the same 'description'."),
          br(),
          fileInput("run_an_simp_ctrl_upload", "Choose controls",
                    multiple = FALSE,
                    accept = c("text/csv",
                               "text/comma-separated-values,text/plain",
                               ".csv")),
          fileInput("run_an_simp_targ_upload", "Choose targets",
                    multiple = FALSE,
                    accept = c("text/csv",
                               "text/comma-separated-values,text/plain",
                               ".csv")),
          actionButton("run_an_simp_assign_res_name", "Assign result ID"),
          br(),
          br(),
          actionButton("run_an_simp_run", "Run analysis")
        ),
        conditionalPanel(
          'input.dataset === "Run analysis (advanced entry)"',
          div("Perform analysis, provided a table containing the controls and targets. Uploaded file must contain columns 'targ_ct_lst', 'ctrl_ct_lst', 'targ_description', and 'ctrl_description'. See example input for guidance."),
          br(),
          fileInput("run_an_adv_upload", "Choose csv file",
                    multiple = FALSE,
                    accept = c("text/csv",
                               "text/comma-separated-values,text/plain",
                               ".csv")),
          numericInput("input_an_adv_ci", "Choose a confidence interval for enrichments (%). Decimals allowed.", 95),
          br(),
          actionButton("run_an_adv_assign_res_name", "Assign result ID"),
          br(),
          br(),
          actionButton("run_an_adv_run", "Run analysis")
        ),
        conditionalPanel(
          'input.dataset === "Completed analyses"',
          div("Here you can view completed analyses and check how analysis experiments were defined."),
          br(),
          selectInput("input_completed_analyses", "Choose a res ID", res_lst)
        ),
        conditionalPanel(
          'input.dataset === "Table viewer"',
          div("Find count and enrichment data for specific compounds."),
          selectInput("input_tbl_viewer_res_id", "Choose res_id", res_lst),
          selectInput("input_tbl_viewer_targ_desc", "Choose target", NULL),
          selectInput("input_tbl_viewer_ctrl_desc", "Choose control", NULL),
          selectInput("input_tbl_viewer_lib_id", "Choose library", NULL),
          actionButton("visualize_tbl_viewer", "View table"),
          br(),
          br(),
          fileInput("tbl_viewer_upload", "Show specific compounds:",
                    multiple = FALSE,
                    accept = c("text/csv",
                               "text/comma-separated-values,text/plain",
                               ".csv")),
          downloadButton("save_tbl_viewer", "Save"),
          br(),
          br(),
          div("Key:"),
          div("er = enrichment ratio (maximum likelihood)"),
          div("er_lb = enrichment ratio (lower bound of confidence interval)"),
          div("er_ub = enrichment ratio (upper bound of confidence interval)")
        ),
        conditionalPanel(
          'input.dataset === "3D enrichment scatterplots (interactive)"',
          selectInput("input_3d_er_interactive_res_id", "Choose res_id", res_lst),
          selectInput("input_3d_er_interactive_targ_desc", "Choose target", NULL),
          selectInput("input_3d_er_interactive_ctrl_desc", "Choose control", NULL),
          selectInput("input_3d_er_interactive_lib_id", "Choose library", NULL),
          br(),
          sliderInput("visualize_3d_er_sp_interactive_num_points", "Number of points to visualize", min = 100, max = 100000, value = 1000, step = 100),
          sliderInput("visualize_3d_er_sp_interactive_scale", "Adjust point size", min = 0, max = 5, value = 1, step = 0.01),
          radioButtons("visualize_3d_er_sp_interactive_scale_log_radio", "Scale point size by:", choices = c("enrichment", "log(enrichment)")),
          actionButton("visualize_3d_er_sp_interactive_load_cpd_struct", "Load structure metadata"),
          actionButton("visualize_3d_er_sp_interactive", "Visualize plot"),
          br(),
          br(),
          div("Each point is a trisynthon where the color and size both correspond to the enrichment ratio (lower bound). Click on points to render structures.")
        ),
        conditionalPanel(
          'input.dataset === "Interactive scatterplots"',
          selectInput("input_interactive_lib_id", "Choose library", lib_lst),
          selectInput("input_interactive_agg_id", "Choose level of aggregation", agg_lst),
          h4("y axis"),
          selectInput("input_interactive_y_type", "Choose type of value", c("count", "rank by counts", "enrichment")),
          selectInput("input_interactive_y_res_id", "Choose res_id", res_lst),
          conditionalPanel(
            'input.input_interactive_y_type == "count" | input.input_interactive_y_type == "rank by counts"',
            selectInput("input_interactive_y_samp_desc", "Choose sample", NULL)
          ),
          conditionalPanel(
            'input.input_interactive_y_type == "enrichment"',
            selectInput("input_interactive_y_targ_desc", "Choose target", NULL),
            selectInput("input_interactive_y_ctrl_desc", "Choose control", NULL)
          ),
          conditionalPanel(
            'input.input_interactive_y_type == "enrichment"',
            checkboxInput("input_interactive_y_log", "log scale y-axis", value = F)
          ),
          h4("x axis"),
          selectInput("input_interactive_x_type", "Choose type of value", c("count", "rank by counts", "enrichment")),
          selectInput("input_interactive_x_res_id", "Choose res_id", res_lst),
          conditionalPanel(
            'input.input_interactive_x_type == "count" | input.input_interactive_x_type == "rank by counts"',
            selectInput("input_interactive_x_samp_desc", "Choose sample", NULL)
          ),
          conditionalPanel(
            'input.input_interactive_x_type == "enrichment"',
            selectInput("input_interactive_x_targ_desc", "Choose target", NULL),
            selectInput("input_interactive_x_ctrl_desc", "Choose control", NULL)
          ),
          conditionalPanel(
            'input.input_interactive_x_type == "enrichment"',
            checkboxInput("input_interactive_x_log", "log scale x-axis", value = F)
          ),
          # br(),
          radioButtons("input_interactive_vis", "Plot visualization", choices = c("Off", "On"), inline = T)
        )
      ),
      
      mainPanel(
        tabsetPanel(
          id = 'dataset',
          tabPanel("lib_meta", DT::dataTableOutput("lib_meta_table")),
          tabPanel("comb_lib_meta", DT::dataTableOutput("comb_lib_meta_table")),
          tabPanel("targ_meta", DT::dataTableOutput("targ_meta_table")),
          tabPanel("resin_meta", DT::dataTableOutput("resin_meta_table")),
          tabPanel("scr_meta", DT::dataTableOutput("scr_meta_table")),
          tabPanel("samp_meta", DT::dataTableOutput("samp_meta_table")),
          tabPanel("merged_meta", fluidPage(fluidRow(DT::dataTableOutput("merged_meta_table")),
                                            br(),
                                            fluidRow(
                                              actionButton("merged_tbl_select_all", "Select all samples"),
                                              actionButton("merged_tbl_clear_all", "Clear sample selection"),
                                              downloadButton("save_merged_selection", "Save selection")
                                              )
                                            )),
          tabPanel("New samples", DT::dataTableOutput("new_samples_table")),
          tabPanel("run_meta", DT::dataTableOutput("run_meta_table")),
          tabPanel("Create config", DT::dataTableOutput("create_config_table")),
          tabPanel("Get structure", fluidPage(h4("Single compound"),
                                              verbatimTextOutput("single_cpd_struct"),
                                              h4("Multiple compounds"),
                                              DT::dataTableOutput("get_structure_table"))),
          tabPanel("Run analysis (simple entry)", fluidPage(fluidRow(DT::dataTableOutput("run_an_simp_ctrl_table"),
                                                                     DT::dataTableOutput("run_an_simp_targ_table")),
                                                    verbatimTextOutput("run_an_simp_result_name"))),
          tabPanel("Run analysis (advanced entry)", fluidPage(DT::dataTableOutput("run_an_adv_table"),
                                                              verbatimTextOutput("run_an_adv_result_name"))),
          tabPanel("Completed analyses", DT::dataTableOutput("res_meta_viewer")),
          tabPanel("Table viewer", DT::dataTableOutput("table_viewer")),
          tabPanel("3D enrichment scatterplots (interactive)", fluidRow(column(width = 8, 
                                                                               tags$head(tags$script('$(document).on("shiny:connected", function(e) {
                            Shiny.onInputChange("main_panel_innerWidth_3d_er_sp_interactive", window.innerWidth);
                            });
                            $(window).resize(function(e) {
                            Shiny.onInputChange("main_panel_innerWidth_3d_er_sp_interactive", window.innerWidth);
                            });
                            ')),
                                                                               plotlyOutput("output_3d_er_sp_interactive", height = "auto"),
                                                                               br(),
                                                                               downloadButton("output_3d_er_sp_interactive_save_plot", "Save plot"),
                                                                               downloadButton("output_3d_er_sp_interactive_save_tbl", "Save data as table")),
                                                                        column(width = 4, 
                                                                               tags$head(tags$script('$(document).on("shiny:connected", function(e) {
                            Shiny.onInputChange("side_panel_innerWidth_3d_er_sp_interactive", window.innerWidth);
                            });
                            $(window).resize(function(e) {
                            Shiny.onInputChange("side_panel_innerWidth_3d_er_sp_interactive", window.innerWidth);
                            });
                            ')),
                                                                               h4("Click on points to view structure (load metadata first)."),
                                                                               plotOutput("output_3d_er_sp_interactive_selected_struct", inline = T),
                                                                               br(),
                                                                               actionButton("input_3d_er_sp_interactive_smi_to_clip", "Copy SMILES to clipboard"),
                                                                               br(),
                                                                               actionButton("input_3d_er_sp_interactive_reorder_bbs", "Reorder by structural similarity"),
                                                                               actionButton("input_3d_er_sp_interactive_reset_order", "Reset order"),
                                                                               br(),
                                                                               verbatimTextOutput("output_3d_er_sp_interactive_selection"),
                                                                               plotlyOutput("output_3d_er_sp_interactive_selected_er")))),
          tabPanel("Interactive scatterplots", fluidPage(conditionalPanel('input.input_interactive_vis == "On"',
                                                                          fluidRow(
                                                                            column(width = 6,
                                                                                   sliderInput("input_interactive_max_points", "Max number of points",
                                                                                               100, 100000, value = 1000, step = 100)),
                                                                            column(width = 6,
                                                                                   sliderInput("visualize_interactive_sp_size", "Adjust height", min = 100, max = 1000, value = 400, step = 50))
                                                                          ),
                                                                          div("Color specific compounds. Separate BBs with semicolons and use hyphens to denote ranges (e.g. 1-5;7;9-12)."),
                                                                          fluidRow(
                                                                            column(width = 2,
                                                                                   textInput("input_interactive_color_cy1", "Cycle 1", NA)),
                                                                            column(width = 2, h4("AND", style="padding:20px;")),
                                                                            column(width = 2,
                                                                                   textInput("input_interactive_color_cy2", "Cycle 2", NA)),
                                                                            column(width = 2, h4("AND", style="padding:20px;")),
                                                                            column(width = 2,
                                                                                   textInput("input_interactive_color_cy3", "Cycle 3", NA)),
                                                                            column(width = 2,
                                                                                   radioButtons("input_interactive_color", "", choices = c("Off", "On")))
                                                                          ),
                                                                          div(""),
                                                                          fluidRow(
                                                                            column(width = 6,
                                                                                   fileInput("input_interactive_color_tbl_request", "Color specific compounds using input table",
                                                                                             multiple = FALSE,
                                                                                             accept = c("text/csv",
                                                                                                        "text/comma-separated-values,text/plain",
                                                                                                        ".csv"))),
                                                                            column(width = 4,
                                                                                   radioButtons("input_interactive_color_tbl", "", choices = c("Off", "On")))
                                                                          )
                                                                          ),
                                                         plotlyOutput("output_interactive_sp", height = "auto"),
                                                         br(),
                                                         conditionalPanel('input.input_interactive_vis == "On"', downloadButton("output_interactive_sp_save_plot", "Save plot"),
                                                                          downloadButton("output_interactive_sp_save_tbl", "Save data as table"))
                                                         ))
        )
      )
    )
  )
}

addResourcePath("tempdir", getwd())
server <- function(input, output, session) {
  rv <- reactiveValues(rv_lib_meta = lib_meta,
                       rv_comb_lib_meta = comb_lib_meta,
                       rv_input_comb_lib = tibble(lib_id = 0, frac = 0),
                       rv_run_meta = run_meta,
                       rv_targ_meta = targ_meta,
                       rv_resin_meta = resin_meta,
                       rv_scr_meta = scr_meta,
                       rv_samp_meta = samp_meta,
                       rv_ct_file_meta = ct_file_meta,
                       rv_merged_meta = ct_file_meta %>%
                         left_join(run_meta) %>%
                         left_join(samp_meta) %>%
                         left_join(scr_meta) %>%
                         left_join(targ_meta) %>%
                         left_join(resin_meta),
                       rv_single_cpd_struct = NULL,
                       rv_cpd_structure = tibble(lib_id = numeric(),
                                                 cycle1 = numeric(),
                                                 cycle2 = numeric(),
                                                 cycle3 = numeric(),
                                                 structure = character()),
                       rv_lib_enum_map = NULL,
                       rv_res_lst = res_lst,
                       rv_run_an_simp_ctrl_lst = c(),
                       rv_run_an_simp_targ_lst = c(),
                       rv_run_an_simp_res_name = NULL,
                       rv_run_an_simp_ctrl_tbl = c(),
                       rv_run_an_simp_targ_tbl = c(),
                       rv_run_an_adv_res_name = NULL,
                       rv_run_an_adv_in_tbl = c(),
                       rv_res_meta_viewer = NULL,
                       rv_tbl_viewer = NULL,
                       rv_3d_er_sp_interactive_tbl = NULL,
                       rv_3d_er_sp_interactive_struct_map = NULL,
                       rv_3d_er_sp_interactive_selected_struct = NULL,
                       rv_3d_er_sp_interactive_x_map = NULL,
                       rv_3d_er_sp_interactive_y_map = NULL,
                       rv_3d_er_sp_interactive_z_map = NULL,
                       rv_interactive_meta = c(),
                       rv_interactive_y_vals = NULL,
                       rv_interactive_x_vals = NULL,
                       rv_interactive_y_samp_ct = NULL,
                       rv_interactive_y_targ_ct = NULL,
                       rv_interactive_y_ctrl_ct = NULL,
                       rv_interactive_x_samp_ct = NULL,
                       rv_interactive_x_targ_ct = NULL,
                       rv_interactive_x_ctrl_ct = NULL,
                       rv_interactive_col_vals = NULL,
                       rv_interactive_tbl = NULL)
  
  shinyjs::hide("input_3d_er_sp_interactive_smi_to_clip")
  shinyjs::hide("input_3d_er_sp_interactive_reorder_bbs")
  shinyjs::hide("input_3d_er_sp_interactive_reset_order")
  shinyjs::hide("output_3d_er_sp_interactive_save_plot")
  shinyjs::hide("output_3d_er_sp_interactive_save_tbl")
  
  output$lib_meta_table <- DT::renderDataTable({
    DT::datatable(rv$rv_lib_meta, rownames = F, filter = "top", selection = "none")
  })
  
  observeEvent(input$append_lib_meta, {
    if(input$user_lib_seq == "") {
      showNotification("Specify the tag sequence associated with the library", type = "warning")
    }
    if(input$user_lib_name == "") {
      showNotification("Specify the name of the library", type = "warning")
    }
    if(input$user_chemist == "") {
      showNotification("Specify the chemist(s) who made the library", type = "warning")
    }
    if(input$user_lib_completion_date == "") {
      showNotification("Specify the date of completion", type = "warning")
    }
    if(input$user_lib_loc == "") {
      showNotification("Specify the current location of the library", type = "warning")
    }
    if(input$user_lib_seq != "" & input$user_lib_name != "" & input$user_chemist != "" &
       input$user_lib_completion_date != "" & input$user_lib_loc != "") {
      temp <- tibble(lib_id = min(setdiff(1:999, rv$rv_lib_meta$lib_id)),
                     lib_seq = input$user_lib_seq,
                     lib_name = input$user_lib_name,
                     chemist = input$user_chemist,
                     lib_completion_date = input$user_lib_completion_date,
                     lib_loc = input$user_lib_loc)
      rv$rv_lib_meta <- rv$rv_lib_meta %>%
        bind_rows(temp)
      showNotification("Entry added")
    }
  })
  
  observeEvent(input$save_lib_meta, {
    write.csv(rv$rv_lib_meta, file.path(config$project_dir, config$library_metadata), row.names = F)
    showNotification("Changes saved")
  })
  
  observeEvent(input$reload_lib_meta, {
    rv$rv_lib_meta <- file.path(config$project_dir, config$library_metadata) %>% fread() %>% as_tibble()
    showNotification("Reloaded")
  })
  
  output$comb_lib_meta_table <- DT::renderDataTable({
    DT::datatable(rv$rv_comb_lib_meta, rownames = F, filter = "top", selection = "none")
  })
  
  output$input_comb_lib <- DT::renderDataTable({
    DT::datatable(rv$rv_input_comb_lib, editable = "cell",
                  rownames = F, selection = "none", options = list(dom = "t"))
  })
  
  observeEvent(input$input_comb_lib_cell_edit, {
    info <- input$input_comb_lib_cell_edit
    edit_row <-  info$row
    edit_col <-  info$col + 1
    edit_value <-  as.numeric(info$value)
    
    rv$rv_input_comb_lib[edit_row, edit_col] <- edit_value
  })
  
  observeEvent(input$input_comb_lib_add_row, {
    rv$rv_input_comb_lib <- rv$rv_input_comb_lib %>%
      bind_rows(tibble(lib_id = 0, frac = 0))
  })
  
  observeEvent(input$input_comb_lib_remove_row, {
    rv$rv_input_comb_lib <- rv$rv_input_comb_lib %>%
      head(-1)
  })
  
  observeEvent(input$append_comb_lib_meta, {
    if(0 %in% rv$rv_input_comb_lib$lib_id) {
      showNotification("Specify all library IDs", type = "warning")
    }
    if(abs(sum(rv$rv_input_comb_lib$frac) - 1) > 0.01) {
      showNotification("Fractional compositions don't add up to 1", type = "warning")
    }
    if(!0 %in% rv$rv_input_comb_lib$lib_id & abs(sum(rv$rv_input_comb_lib$frac) - 1) < 0.01) {
      temp <- rv$rv_input_comb_lib %>%
        mutate(comb_lib_id = min(setdiff(1:99999, rv$rv_comb_lib_meta$comb_lib_id)))
      rv$rv_comb_lib_meta <- rv$rv_comb_lib_meta %>%
        bind_rows(temp)
      showNotification("Entry added")
    }
  })
  
  observeEvent(input$save_comb_lib_meta, {
    write.csv(rv$rv_comb_lib_meta, file.path(config$project_dir, config$combined_library_metadata), row.names = F)
    showNotification("Changes saved")
  })
  
  observeEvent(input$reload_comb_lib_meta, {
    rv$rv_comb_lib_meta <- file.path(config$project_dir, config$combined_library_metadata) %>% fread() %>% as_tibble()
    showNotification("Reloaded")
  })
  
  output$run_meta_table <- DT::renderDataTable({
    DT::datatable(rv$rv_run_meta, rownames = F, filter = "top", selection = "none")
  })

  output$targ_meta_table <- DT::renderDataTable({
    DT::datatable(rv$rv_targ_meta, rownames = F, filter = "top", selection = "none")
  })
  
  observeEvent(input$append_targ_meta, {
    if(input$user_targ_name == "") {
      showNotification("Specify the name of the target", type = "warning")
    }
    if(input$user_targ_name != "") {
      temp <- tibble(targ_id = min(setdiff(1:99999, rv$rv_targ_meta$targ_id)),
                     targ_name = input$user_targ_name)
      rv$rv_targ_meta <- rv$rv_targ_meta %>%
        bind_rows(temp)
      showNotification("Entry added")
    }
  })
  
  observeEvent(input$save_targ_meta, {
    write.csv(rv$rv_targ_meta, file.path(config$project_dir, config$target_metadata), row.names = F)
    showNotification("Changes saved")
  })
  
  observeEvent(input$reload_targ_meta, {
    rv$rv_targ_meta <- file.path(config$project_dir, config$target_metadata) %>% fread() %>% as_tibble()
    showNotification("Reloaded")
  })
  
  output$resin_meta_table <- DT::renderDataTable({
    DT::datatable(rv$rv_resin_meta, rownames = F, filter = "top", selection = "none")
  })
  
  observeEvent(input$append_resin_meta, {
    if(input$user_resin_name == "") {
      showNotification("Specify the name of the resin", type = "warning")
    }
    if(input$user_resin_name != "") {
      temp <- tibble(resin_id = min(setdiff(1:99999, rv$rv_resin_meta$resin_id)),
                     resin_name = input$user_resin_name)
      rv$rv_resin_meta <- rv$rv_resin_meta %>%
        bind_rows(temp)
      showNotification("Entry added")
    }
  })
  
  observeEvent(input$save_resin_meta, {
    write.csv(rv$rv_resin_meta, file.path(config$project_dir, config$resin_metadata), row.names = F)
    showNotification("Changes saved")
  })
  
  observeEvent(input$reload_resin_meta, {
    rv$rv_resin_meta <- file.path(config$project_dir, config$resin_metadata) %>% fread() %>% as_tibble()
    showNotification("Reloaded")
  })
  
  output$scr_meta_table <- DT::renderDataTable({
    DT::datatable(rv$rv_scr_meta, rownames = F, filter = "top", selection = "none")
  })
  
  observeEvent(input$append_scr_meta, {
    if(input$user_scr_name == "") {
      showNotification("Specify the name of the DEL screen", type = "warning")
    }
    if(input$user_scr_completion_date == "") {
      showNotification("Specify the completion date for the DEL screen", type = "warning")
    }
    if(input$user_scr_name != "" & input$user_scr_completion_date != "") {
      temp <- tibble(scr_id = min(setdiff(1:99999, rv$rv_scr_meta$scr_id)),
                     scr_name = input$user_scr_name,
                     scr_completion_date = input$user_scr_completion_date)
      rv$rv_scr_meta <- rv$rv_scr_meta %>%
        bind_rows(temp)
      showNotification("Entry added")
    }
  })
  
  observeEvent(input$save_scr_meta, {
    write.csv(rv$rv_scr_meta, file.path(config$project_dir, config$screen_metadata), row.names = F)
    showNotification("Changes saved")
  })
  
  observeEvent(input$reload_scr_meta, {
    rv$rv_scr_meta <- file.path(config$project_dir, config$screen_metadata) %>% fread() %>% as_tibble()
    showNotification("Reloaded")
  })
  
  output$samp_meta_table <- DT::renderDataTable({
    DT::datatable(rv$rv_samp_meta, rownames = F, filter = "top", selection = "none")
  })
  
  observeEvent(input$save_samp_meta, {
    write.csv(rv$rv_samp_meta, file.path(config$project_dir, config$sample_metadata), row.names = F)
    showNotification("Changes saved")
  })
  
  observeEvent(input$reload_samp_meta, {
    rv$rv_samp_meta <- file.path(config$project_dir, config$sample_metadata) %>% fread() %>% as_tibble()
    showNotification("Reloaded")
  })
  
  observeEvent(input$append_samp_meta, {
    unused_samp_ids <- setdiff(1:999999, rv$rv_samp_meta$samp_id)
    temp <- as_tibble(fread(input$new_samples$datapath)) %>%
      mutate(samp_id = unused_samp_ids[1:nrow(.)])
    rv$rv_samp_meta <- rv$rv_samp_meta %>%
      bind_rows(temp)
    showNotification("Entries added")
  })
  
  observeEvent(input$append_run_meta, {
    if(input$user_seq_type == "") {
      showNotification("Specify the type of sequencing (e.g. MiSeq, HiSeq)", type = "warning")
    }
    if(input$user_seq_location == "") {
      showNotification("Specify the location of the sequencing", type = "warning")
    }
    if(input$user_seq_completion_date == "") {
      showNotification("Specify the sequencing was completed", type = "warning")
    }
    if(input$user_fastq_folder_name == "") {
      showNotification("Specify the name of the folder holding the fastq files", type = "warning")
    }
    if(input$user_seq_type != "" & input$user_seq_location != "" & input$user_fastq_folder_name != "" & input$user_fastq_folder_name != "") {
      temp <- tibble(run_id = min(setdiff(1:10000, rv$rv_run_meta$run_id)),
                     seq_type = input$user_seq_type,
                     seq_location = input$user_seq_location,
                     seq_completion_date = input$user_seq_completion_date,
                     fastq_folder_name = input$user_fastq_folder_name)
      rv$rv_run_meta <- rv$rv_run_meta %>%
        bind_rows(temp)
      showNotification("Entry added")
    }
  })
  
  observeEvent(input$save_run_meta, {
    write.csv(rv$rv_run_meta, file.path(config$project_dir, config$sequencing_run_metadata), row.names = F)
    showNotification("Changes saved")
  })
  
  observeEvent(input$reload_run_meta, {
    rv$rv_run_meta <- file.path(config$project_dir, config$sequencing_run_metadata) %>% fread() %>% as_tibble()
    showNotification("Reloaded")
  })
  
  output$merged_meta_table <- DT::renderDataTable({
    DT::datatable(rv$rv_merged_meta[, input$show_vars, drop = FALSE] %>%
                    mutate(across(everything(), as.character)),
                  rownames = F, filter = "top")
  })
  
  observeEvent(input$merged_select_all_cols, {
    updateCheckboxGroupInput(session, "show_vars", choices = names(rv$rv_merged_meta), selected = names(rv$rv_merged_meta))
  })
  
  observeEvent(input$merged_clear_all_cols, {
    updateCheckboxGroupInput(session, "show_vars", choices = names(rv$rv_merged_meta), selected = NULL)
  })
  
  merged_proxy <- dataTableProxy("merged_meta_table")
  
  observeEvent(input$reload_merged_meta, {
    rv$rv_merged_meta <- rv$rv_ct_file_meta %>%
      left_join(rv$rv_run_meta) %>%
      left_join(rv$rv_samp_meta) %>%
      left_join(rv$rv_targ_meta) %>%
      left_join(rv$rv_resin_meta)
    showNotification("Reloaded")
  })
  
  observeEvent(input$merged_tbl_select_all, {
    merged_proxy %>% selectRows(input$merged_meta_table_rows_all)
    showNotification("All samples selected")
  })
  
  observeEvent(input$merged_tbl_clear_all, {
    merged_proxy %>% selectRows(NULL)
    showNotification("Selection cleared")
  })
  
  output$save_merged_selection <- downloadHandler("csv_with_selection.csv", 
                                                  content = function(file) {
                                                    write.csv(rv$rv_merged_meta[input$merged_meta_table_rows_selected,], file, row.names = F)
                                                  },
                                                  contentType = "text/csv")
  
  output$new_samples_table <- DT::renderDataTable({
    if(is.null(input$new_samples)) {
      return(NULL)
    }
    DT::datatable(fread(input$new_samples$datapath),
                  rownames = F, selection = "none")
  })
  
  output$create_config_table <- DT::renderDataTable({
    if(is.null(input$config_upload)) {
      return(NULL)
    }
    DT::datatable(fread(input$config_upload$datapath),
                  rownames = F, selection = "none")
  })
  
  output$download_json <- downloadHandler("config_name.json", 
                                          content = function(file) {
                                            temp_out <- as_tibble(fread(input$config_upload$datapath)) %>%
                                              mutate(id = paste0("run", sprintf("%003d", run_id), "_samp", sprintf("%006d", samp_id)))
                                            exp_list <- list()
                                            for(i in 1:nrow(temp_out)) {
                                              temp_exp_lib_ids <- rv$rv_comb_lib_meta %>%
                                                filter(comb_lib_id == pull(temp_out, comb_lib_id)[i]) %>%
                                                pull(lib_id)
                                              temp_exp <- list(
                                                id = pull(temp_out, id)[i] %>% unbox(),
                                                fname = pull(temp_out, seq_filename)[i] %>% str_split(", ") %>% unlist(),
                                                library_ids = temp_exp_lib_ids
                                              )
                                              exp_list[[as.character(i-1)]] <- temp_exp
                                            }
                                            lib_list <- rv$rv_comb_lib_meta %>%
                                              filter(comb_lib_id %in% temp_out$comb_lib_id) %>%
                                              pull(lib_id) %>%
                                              unique()
                                            meta_list <- list(
                                              library = '../lib_meta.csv' %>% unbox(),
                                              encoding = '../bb_meta.csv' %>% unbox(),
                                              encoding_enumerated = '../lib_enum_struct.csv' %>% unbox(),
                                              library_ids = lib_list
                                            )
                                            seq_slices <- list(
                                              c(0,9),
                                              c(9,18),
                                              c(18,27),
                                              c(27,37),
                                              c(37,50)
                                            )
                                            settings_list <- list(
                                              fuzzy = unbox(T),
                                              umi = "directional" %>% unbox(),
                                              minread = 1 %>% unbox(),
                                              seq_slices = seq_slices,
                                              read_subs = {}
                                            )
                                            
                                            cfg <- list(
                                              project_dir = NA %>% unbox(),
                                              datasets = exp_list,
                                              enrichments = {},
                                              metadata = meta_list,
                                              settings = settings_list,
                                              save_often = 5 %>% unbox(),
                                              # outfile = str_remove_all(input$config_upload$name, ".csv") %>% unbox(),
                                              outfile = "save_file" %>% unbox(),
                                              reprocess_reads = F %>% unbox(),
                                              reprocess_counts = F %>% unbox(),
                                              reprocess_enrichments = T %>% unbox(),
                                              njobs = 4 %>% unbox()
                                            )
                                            
                                            write(toJSON(cfg, pretty = T), file)
                                          },
                                          contentType = "json")
  
  observeEvent(input$update_ct_file_meta, {
    temp_out <- as_tibble(fread(input$config_upload$datapath)) %>%
      left_join(select(rv$rv_comb_lib_meta, comb_lib_id, lib_id)) %>%
      mutate(ct_filename = paste0("run", sprintf("%003d", run_id), 
                                  "_samp", sprintf("%006d", samp_id), 
                                  "_lib", sprintf("%003d", lib_id),
                                  ".csv")) %>%
      select(run_id, samp_id, lib_id, seq_filename, ct_filename)
    rv$rv_ct_file_meta <- rv$rv_ct_file_meta %>%
      bind_rows(temp_out) %>%
      distinct()
    
    write.csv(rv$rv_ct_file_meta, file.path(config$project_dir, config$count_filename_metadata), row.names = F)
    
    showNotification("ct_file_meta updated")
  })
  
  output$single_cpd_struct <- renderText({
    if(is.null(rv$rv_single_cpd_struct)) {
      "NULL"
    } else {
      unname(rv$rv_single_cpd_struct)
    }
  })

  output$get_structure_table <- DT::renderDataTable({
    DT::datatable(rv$rv_cpd_structure, rownames = F, selection = "none",
                  options = list(dom = "tp"))
  })
  
  observeEvent(input$cpd_struct_upload, {
    temp <- as_tibble(fread(input$cpd_struct_upload$datapath))
    rv$rv_cpd_structure <- rv$rv_cpd_structure %>%
      bind_rows(temp)
  })
  
  observeEvent(input$load_cpd_struct, {
    if(length(rv$rv_lib_enum_map) > 0) {
      showNotification("Metadata already loaded")
    }
    if(length(rv$rv_lib_enum_map) == 0) {
      load_cpd_notif <- showNotification("Loading metadata...", duration = NULL)
      temp <- file.path(config$project_dir, config$enumerated_structure_metadata) %>% fread() %>% as_tibble()
      rv$rv_lib_enum_map <- structure(temp$structure, names = paste(temp$lib_id,
                                                                    temp$cycle1,
                                                                    temp$cycle2,
                                                                    temp$cycle3))
      removeNotification(load_cpd_notif)
      showNotification("Finished loading metadata")
    }
  })
  
  observeEvent(input$get_cpd_struct, {
    if(length(rv$rv_lib_enum_map) == 0) {
      showNotification("Load metadata first", type = "warning")
    }
    if(length(rv$rv_lib_enum_map) > 0) {
      rv$rv_single_cpd_struct <- rv$rv_lib_enum_map[paste(input$input_cpd_struct_single_lib_id,
                                                          input$input_cpd_struct_single_cy1,
                                                          input$input_cpd_struct_single_cy2,
                                                          input$input_cpd_struct_single_cy3)]
    }
  })
  
  observeEvent(input$get_cpd_structs, {
    if(length(rv$rv_lib_enum_map) == 0) {
      showNotification("Load metadata first", type = "warning")
    }
    if(length(rv$rv_lib_enum_map) > 0) {
      rv$rv_cpd_structure$structure <- rv$rv_lib_enum_map[paste(rv$rv_cpd_structure$lib_id,
                                                                rv$rv_cpd_structure$cycle1,
                                                                rv$rv_cpd_structure$cycle2,
                                                                rv$rv_cpd_structure$cycle3)]
    }
  })
  
  observeEvent(input$clear_cpd_struct, {
    rv$rv_cpd_structure <- tibble(lib_id = numeric(),
                                  cycle1 = numeric(),
                                  cycle2 = numeric(),
                                  cycle3 = numeric(),
                                  structure = character())
  })
  
  output$save_cpd_struct <- downloadHandler("csv_with_structures.csv", 
                                            content = function(file) {
                                              write.csv(rv$rv_cpd_structure, file, row.names = F)
                                            },
                                            contentType = "text/csv")
  
  output$run_an_simp_result_name <- renderPrint({
    list(
      res_name = rv$rv_run_an_simp_res_name
    )
  })
  
  output$run_an_simp_ctrl_table <- DT::renderDataTable({
    DT::datatable(rv$rv_run_an_simp_ctrl_tbl, rownames = F, width = "50%", selection = "none",
                  options = list(dom = 'tp',
                                 pageLength = 5,
                                 lengthMenu = c(5, 10, 15, 20)))
  })
  
  output$run_an_simp_targ_table <- DT::renderDataTable({
    DT::datatable(rv$rv_run_an_simp_targ_tbl, rownames = F, width = "50%", selection = "none",
                  options = list(dom = 'tp',
                                 pageLength = 5,
                                 lengthMenu = c(5, 10, 15, 20)))
  })
  
  observeEvent(input$run_an_simp_ctrl_upload, {
    rv$rv_run_an_simp_ctrl_tbl <- as_tibble(fread(input$run_an_simp_ctrl_upload$datapath)) %>%
      select(ct_filename, description)
    rv$rv_run_an_simp_ctrl_lst <- rv$rv_run_an_simp_ctrl_tbl %>%
      pull(ct_filename)
  })
  
  observeEvent(input$run_an_simp_targ_upload, {
    rv$rv_run_an_simp_targ_tbl <- as_tibble(fread(input$run_an_simp_targ_upload$datapath)) %>%
      select(ct_filename, description)
    rv$rv_run_an_simp_targ_lst <- rv$rv_run_an_simp_targ_tbl %>%
      pull(ct_filename)
  })
  
  observeEvent(input$run_an_simp_assign_res_name, {
    temp <- list.files("./meta/res/", recursive = F, full.names = F) %>% 
      str_remove_all(".csv") %>%
      str_remove_all("res") %>% 
      as.numeric()
    unused_res_ids <- setdiff(1:9999, temp)
    rv$rv_run_an_simp_res_name <- paste0("res", sprintf("%04d", min(unused_res_ids)))
  })
  
  observeEvent(input$run_an_simp_run, {
    if(is.null(rv$rv_run_an_simp_ctrl_lst)) {
      showNotification("Specify control samples", type = "warning")
    }
    if(is.null(rv$rv_run_an_simp_targ_lst)) {
      showNotification("Specify target samples", type = "warning")
    }
    if(is.null(rv$rv_run_an_simp_res_name)) {
      showNotification("Assign a result ID", type = "warning")
    }
    if(!is.null(rv$rv_run_an_simp_targ_lst) &
       !is.null(rv$rv_run_an_simp_ctrl_lst) &
       !is.null(rv$rv_run_an_simp_res_name)) {
      run_an_simp_start_notif <- showNotification("Running analysis...", duration = NULL)
      
      ptm <- proc.time()
      report_name <- format(Sys.time(), "%Y%m%d_%H%M%S")
      report_lines <- c()
      
      # Import metadata
      source(file.path(config$project_dir, config$metadata_script))
      
      # input from user - lists for tgts and ctrls
      targ_ct_lst <- rv$rv_run_an_simp_targ_lst
      ctrl_ct_lst <- rv$rv_run_an_simp_ctrl_lst
      res_name <- rv$rv_run_an_simp_res_name
      
      submeta <- run_meta %>%
        left_join(ct_file_meta) %>%
        left_join(samp_meta) %>%
        left_join(targ_meta) %>%
        right_join(bind_rows(rv$rv_run_an_simp_targ_tbl, rv$rv_run_an_simp_ctrl_tbl)) %>%
        group_by_at(setdiff(names(.), c("ct_filename", "samp_id"))) %>%
        mutate(group_id = cur_group_id()) %>%
        ungroup()
      
      # generate replicate scatter plots
      temp_ptm <- proc.time()
      gids <- unique(submeta$group_id)
      for(gid in gids) {
        ids <- submeta %>%
          filter(group_id == gid)
        if(nrow(ids) > 1) {
          for(i in 1:(nrow(ids)-1)) {
            for(j in (i+1):nrow(ids)) {
              generate_replicate_scatterplots(
                id1 = ids$ct_filename[i],
                id2 = ids$ct_filename[j],
                xlab = ids$description[i],
                ylab = ids$description[j],
                plot_title = "replicate scatterplot", 
                width = 6, # in inches
                height = 6, # in inches
                out_loc = file.path(config$project_dir, config$result_location, res_name, "plot")
              )
            }
          }
        }
      }
      temp_dur <- round((proc.time() - temp_ptm)[[3]],2)
      report_lines <- c(report_lines, paste0("Time to generate replicate scatter plots: ", temp_dur, " seconds"))
      
      # parse through input from user and automatically generate relevant analysis experiments
      gids <- submeta %>%
        filter(ct_filename %in% targ_ct_lst) %>%
        pull(group_id) %>%
        unique()
      
      ctrl_desc_lst <- submeta %>%
        filter(ct_filename %in% ctrl_ct_lst) %>%
        pull(description) %>%
        unique()
      
      curr_targ_ct_lst <- c()
      targ_desc_lst <- c()
      for(gid in gids) {
        temp <- submeta %>%
          filter(group_id == gid)
        tids <- temp %>%
          pull(ct_filename)
        tdesc <- temp %>%
          pull(description)
        curr_targ_ct_lst <- c(curr_targ_ct_lst, paste0(tids, collapse = ", "))
        targ_desc_lst <- c(targ_desc_lst, tdesc[1])
      }
      
      in_tbl <- tibble(
        targ_ct_lst = curr_targ_ct_lst,
        ctrl_ct_lst = paste0(ctrl_ct_lst, collapse = ", "),
        targ_description = targ_desc_lst,
        ctrl_description = ctrl_desc_lst,
        agg_id = paste0(agg_meta$agg_id, collapse = ", "),
        an_type_id = paste0(an_type_meta$an_type_id[an_type_meta$an_type_name == "enrichment"], collapse = ", "),
        ci = 0.99
      )
      
      # generate/assign analysis ids and update an_meta
      in_tbl <- assign_an_ids(in_tbl)
      
      # perform analyses to get lists of values for plotting
      if(nrow(in_tbl) > 0) {
        temp_ptm <- proc.time()
        prfm_an_expts(in_tbl)
        temp_dur <- round((proc.time() - temp_ptm)[[3]],2)
        report_lines <- c(report_lines, paste0("Time to perform analysis experiments: ", temp_dur, " seconds"))
      }
      
      temp_dur <- round((proc.time() - ptm)[[3]],2)
      report_lines <- c(report_lines, paste0("Time to complete script: ", temp_dur, " seconds"))
      
      write_lines(report_lines, file.path(config$project_dir, config$result_location, res_name, paste0("report_", report_name, ".txt")))
      
      in_tbl %>%
        write.csv(file.path(config$project_dir, config$result_metadata_location, paste0(res_name, ".csv")), row.names = F)
      
      rv$rv_res_lst <- list.files(file.path(config$project_dir, config$result_metadata_location), 
                                  recursive = F, full.names = F) %>% str_remove_all(".csv")
          
      removeNotification(run_an_simp_start_notif)
      showNotification("Finished running analysis", duration = NULL)
    }
  })
  
  output$run_an_adv_result_name <- renderPrint({
    list(
      res_name = rv$rv_run_an_adv_res_name
    )
  })
  
  output$run_an_adv_table <- DT::renderDataTable({
    DT::datatable(rv$rv_run_an_adv_tbl, rownames = F, width = "50%", selection = "none",
                  options = list(dom = 'tp',
                                 pageLength = 5,
                                 lengthMenu = c(5, 10, 15, 20)))
  })
  
  observeEvent(input$run_an_adv_upload, {
    rv$rv_run_an_adv_tbl <- as_tibble(fread(input$run_an_adv_upload$datapath)) %>%
      select(targ_ct_lst, ctrl_ct_lst, targ_description, ctrl_description)
  })
  
  observeEvent(input$run_an_adv_assign_res_name, {
    temp <- list.files(file.path(config$project_dir, config$result_metadata_location), recursive = F, full.names = F) %>% 
      str_remove_all(".csv") %>%
      str_remove_all("res") %>% 
      as.numeric()
    unused_res_ids <- setdiff(1:9999, temp)
    rv$rv_run_an_adv_res_name <- paste0("res", sprintf("%04d", min(unused_res_ids)))
  })
  
  observe({
    if(is.numeric(input$input_an_adv_ci)) {
      if(input$input_an_adv_ci >= 100) {
        updateNumericInput(session, "input_an_adv_ci", value = 99.99)
      }
      else if(input$input_an_adv_ci < 0) {
        updateNumericInput(session, "input_an_adv_ci", value = 0)
      }
    }
  })
  
  observeEvent(input$run_an_adv_run, {
    if(!"ctrl_ct_lst" %in% names(rv$rv_run_an_adv_tbl)) {
      showNotification("Specify control samples", type = "warning")
    }
    if(!"targ_ct_lst" %in% names(rv$rv_run_an_adv_tbl)) {
      showNotification("Specify target samples", type = "warning")
    }
    if(is.null(rv$rv_run_an_adv_res_name)) {
      showNotification("Assign a result ID", type = "warning")
    }
    if("ctrl_ct_lst" %in% names(rv$rv_run_an_adv_tbl) &
       "targ_ct_lst" %in% names(rv$rv_run_an_adv_tbl) &
       !is.null(rv$rv_run_an_adv_res_name)) {
      run_an_adv_start_notif <- showNotification("Running analysis...", duration = NULL)
      
      ptm <- proc.time()
      report_name <- format(Sys.time(), "%Y%m%d_%H%M%S")
      report_lines <- c()
      
      # Import metadata
      source(file.path(config$project_dir, config$metadata_script))
      
      res_name <- rv$rv_run_an_adv_res_name
      
      # generate replicate scatter plots
      temp_ptm <- proc.time()
      for(i in 1:nrow(rv$rv_run_an_adv_tbl)) {
        if(grepl(", ", rv$rv_run_an_adv_tbl$ctrl_ct_lst[i])) {
          ids <- unlist(str_split(rv$rv_run_an_adv_tbl$ctrl_ct_lst[i], ", "))
          for(i in 1:(length(ids)-1)) {
            for(j in (i+1):length(ids)) {
              generate_replicate_scatterplots(
                id1 = ids[i],
                id2 = ids[j],
                xlab = rv$rv_run_an_adv_tbl$ctrl_description[i],
                ylab = rv$rv_run_an_adv_tbl$ctrl_description[i],
                plot_title = "replicate scatterplot", 
                width = 6, # in inches
                height = 6, # in inches
                out_loc = file.path(config$project_dir, config$result_location, res_name, "plot")
              )
            }
          }
        }
        if(grepl(", ", rv$rv_run_an_adv_tbl$targ_ct_lst[i])) {
          ids <- unlist(str_split(rv$rv_run_an_adv_tbl$targ_ct_lst[i], ", "))
          for(i in 1:(length(ids)-1)) {
            for(j in (i+1):length(ids)) {
              generate_replicate_scatterplots(
                id1 = ids[i],
                id2 = ids[j],
                xlab = rv$rv_run_an_adv_tbl$targ_description[i],
                ylab = rv$rv_run_an_adv_tbl$targ_description[i],
                plot_title = "replicate scatterplot", 
                width = 6, # in inches
                height = 6, # in inches
                out_loc = file.path(config$project_dir, config$result_location, res_name, "plot")
              )
            }
          }
        }
      }
      temp_dur <- round((proc.time() - temp_ptm)[[3]],2)
      report_lines <- c(report_lines, paste0("Time to generate replicate scatter plots: ", temp_dur, " seconds"))
      
      # parse through analysis experiments
      in_tbl <- rv$rv_run_an_adv_tbl %>%
        mutate(agg_id = paste0(agg_meta$agg_id, collapse = ", "),
               an_type_id = paste0(an_type_meta$an_type_id[an_type_meta$an_type_name == "enrichment"], collapse = ", "),
               ci = input$input_an_adv_ci/100)
      
      # generate/assign analysis ids and update an_meta
      in_tbl <- assign_an_ids(in_tbl)
      
      # perform analyses to get lists of values for plotting
      if(nrow(in_tbl) > 0) {
        temp_ptm <- proc.time()
        prfm_an_expts(in_tbl)
        temp_dur <- round((proc.time() - temp_ptm)[[3]],2)
        report_lines <- c(report_lines, paste0("Time to perform analysis experiments: ", temp_dur, " seconds"))
      }
      
      temp_dur <- round((proc.time() - ptm)[[3]],2)
      report_lines <- c(report_lines, paste0("Time to complete script: ", temp_dur, " seconds"))
      
      write_lines(report_lines, file.path(config$project_dir, config$result_location, res_name, paste0("report_", report_name, ".txt")))
      
      in_tbl %>%
        write.csv(file.path(config$project_dir, config$result_metadata_location, paste0(res_name, ".csv")), row.names = F)
      
      rv$rv_res_lst <- list.files(file.path(config$project_dir, config$result_metadata_location), 
                                  recursive = F, full.names = F) %>% str_remove_all(".csv")
      
      removeNotification(run_an_adv_start_notif)
      showNotification("Finished running analysis", duration = NULL)
    }
  })
  
  
  
  
  
  
  
  
  
  
  
  
  observeEvent(input$input_completed_analyses, {
    req(input$input_completed_analyses)
    rv$rv_res_meta_viewer <- file.path(config$project_dir, config$result_metadata_location, paste0(input$input_completed_analyses, ".csv")) %>%
      fread() %>%
      as_tibble()
  })
  
  observe({
    updateSelectInput(session, "input_completed_analyses", choices = rv$rv_res_lst)
  })
  
  output$res_meta_viewer <- DT::renderDataTable({
    req(rv$rv_res_meta_viewer)
    DT::datatable(rv$rv_res_meta_viewer %>%
                    select(-agg_id, -an_type_id, -an_id), 
                  rownames = F, filter = "top", selection = "none")
  })
  
  
  
  
  
  
  
  
  
  
  
  
  
  input_tbl_viewer_targ_choices <- reactive({
    req(input$input_tbl_viewer_res_id)
    as_tibble(fread(file.path(config$project_dir, config$result_metadata_location, paste0(input$input_tbl_viewer_res_id, ".csv")))) %>%
      pull(targ_description)
  })
  
  observe({
    updateSelectInput(session, "input_tbl_viewer_targ_desc", choices = input_tbl_viewer_targ_choices())
  })
  
  input_tbl_viewer_ctrl_choices <- reactive({
    req(input$input_tbl_viewer_res_id)
    as_tibble(fread(file.path(config$project_dir, config$result_metadata_location, paste0(input$input_tbl_viewer_res_id, ".csv")))) %>%
      filter(targ_description == input$input_tbl_viewer_targ_desc) %>%
      pull(ctrl_description)
  })
  
  observe({
    updateSelectInput(session, "input_tbl_viewer_ctrl_desc", choices = input_tbl_viewer_ctrl_choices())
  })
  
  input_tbl_viewer_lib_choices <- reactive({
    req(input$input_tbl_viewer_targ_desc)
    req(input$input_tbl_viewer_ctrl_desc)
    req(input$input_tbl_viewer_res_id)
    new_choices <- as_tibble(fread(file.path(config$project_dir, config$result_metadata_location, paste0(input$input_tbl_viewer_res_id, ".csv")))) %>%
      filter(targ_description == input$input_tbl_viewer_targ_desc & ctrl_description == input$input_tbl_viewer_ctrl_desc)
    if(nrow(new_choices) > 0) {
      new_choices <- new_choices %>%
        select(ctrl_ct_lst, targ_ct_lst) %>%
        unlist() %>%
        str_split(", ") %>%
        unlist() %>%
        unique()
      new_choices <- sapply(strsplit(new_choices, "_"), `[`, 3) %>%
        str_remove_all("lib|\\.csv") %>%
        as.numeric() %>%
        unique()
      tibble(lib_id = new_choices) %>%
        left_join(rv$rv_lib_meta) %>%
        mutate(temp = paste0(lib_id, ": ", lib_name)) %>%
        pull(temp)
    }
  })
  
  observe({
    updateSelectInput(session, "input_tbl_viewer_lib_id", choices = input_tbl_viewer_lib_choices())
  })
  
  observe({
    updateSelectInput(session, "input_tbl_viewer_res_id", "Choose res_id", rv$rv_res_lst)
  })
  
  output$table_viewer <- DT::renderDataTable({
    req(rv$rv_tbl_viewer)
    DT::datatable(rv$rv_tbl_viewer %>%
                    select(-lib_id) %>%
                    mutate(across(matches("cycle"), as.character)), 
                  rownames = F, filter = "top", selection = "none")
  })
  
  observeEvent(input$tbl_viewer_upload, {
    temp_notif <- showNotification("Selecting compounds...", duration = NULL)
    temp <- as_tibble(fread(input$tbl_viewer_upload$datapath))
    rv$rv_tbl_viewer <- rv$rv_tbl_viewer %>%
      right_join(temp)
    removeNotification(temp_notif)
  })
  
  observeEvent(input$visualize_tbl_viewer, {
    vis_tbl_viewer_notif <- showNotification("Gathering table...", duration = NULL)
    temp_tbl <- as_tibble(fread(file.path(config$project_dir, config$result_metadata_location, paste0(input$input_tbl_viewer_res_id, ".csv")))) %>%
      filter(targ_description == input$input_tbl_viewer_targ_desc & ctrl_description == input$input_tbl_viewer_ctrl_desc)
    
    temp_ctrl_lst <- temp_tbl %>% pull(ctrl_ct_lst) %>% str_split(", ") %>% unlist()
    temp_ctrl_ct <- as_tibble(fread(file.path(config$project_dir, config$data_count_location, temp_ctrl_lst[1])))
    if(length(temp_ctrl_lst) > 1) {
      for(i in 2:length(temp_ctrl_lst)) {
        temp_ctrl_ct <- as_tibble(temp_ctrl_ct + as_tibble(fread(file.path(config$project_dir, config$data_count_location, temp_ctrl_lst[i]))))
      }
    }
    temp_targ_lst <- temp_tbl %>% pull(targ_ct_lst) %>% str_split(", ") %>% unlist()
    temp_targ_ct <- as_tibble(fread(file.path(config$project_dir, config$data_count_location, temp_targ_lst[1])))
    if(length(temp_targ_lst) > 1) {
      for(i in 2:length(temp_targ_lst)) {
        temp_targ_ct <- as_tibble(temp_targ_ct + as_tibble(fread(file.path(config$project_dir, config$data_count_location, temp_targ_lst[i]))))
      }
    }
    temp_an_id <- temp_tbl %>%
      pull(an_id)
    temp_ci <- temp_tbl %>%
      pull(ci)
    temp_lib_id <- as.numeric(unlist(str_split(input$input_tbl_viewer_lib_id, ": "))[1])
    temp_er <- as_tibble(fread(file.path(config$project_dir, config$result_analysis_location, 
                                         paste0("an", sprintf("%06d", temp_an_id),
                                                "_lib", sprintf("%03d", temp_lib_id),
                                                "_agg07_type01.csv")))) %>%
      pull(value) %>%
      round(2)
    temp_er_lb <- as_tibble(fread(file.path(config$project_dir, config$result_analysis_location, 
                                            paste0("an", sprintf("%06d", temp_an_id),
                                                   "_lib", sprintf("%03d", temp_lib_id),
                                                   "_agg07_type02.csv")))) %>%
      pull(value) %>%
      round(2)
    temp_er_ub <- as_tibble(fread(file.path(config$project_dir, config$result_analysis_location, 
                                            paste0("an", sprintf("%06d", temp_an_id),
                                                   "_lib", sprintf("%03d", temp_lib_id),
                                                   "_agg07_type03.csv")))) %>%
      pull(value) %>%
      round(2)
    rv$rv_tbl_viewer <- as_tibble(fread(file.path(config$project_dir, config$library_metadata_column_location,
                                               sprintf("lib%03d_agg07_meta_cols.csv", temp_lib_id)))) %>%
      mutate(ctrl_ct = temp_ctrl_ct %>% pull(value),
             targ_ct = temp_targ_ct %>% pull(value),
             er = temp_er,
             er_lb = temp_er_lb,
             er_ub = temp_er_ub)
    removeNotification(vis_tbl_viewer_notif)
  })
  
  output$save_tbl_viewer <- downloadHandler("table_name.csv", 
                                            content = function(file) {
                                              write.csv(rv$rv_tbl_viewer, file, row.names = F)
                                            },
                                            contentType = "text/csv")
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  

  input_3d_er_interactive_targ_choices <- reactive({
    req(input$input_3d_er_interactive_res_id)
    as_tibble(fread(file.path(config$project_dir, config$result_metadata_location, paste0(input$input_3d_er_interactive_res_id, ".csv")))) %>%
      pull(targ_description)
  })
  
  observe({
    updateSelectInput(session, "input_3d_er_interactive_targ_desc", choices = input_3d_er_interactive_targ_choices())
  })
  
  input_3d_er_interactive_ctrl_choices <- reactive({
    req(input$input_3d_er_interactive_res_id)
    as_tibble(fread(file.path(config$project_dir, config$result_metadata_location, paste0(input$input_3d_er_interactive_res_id, ".csv")))) %>%
      filter(targ_description == input$input_3d_er_interactive_targ_desc) %>%
      pull(ctrl_description)
  })
  
  observe({
    updateSelectInput(session, "input_3d_er_interactive_ctrl_desc", choices = input_3d_er_interactive_ctrl_choices())
  })
  
  input_3d_er_interactive_lib_choices <- reactive({
    req(input$input_3d_er_interactive_targ_desc)
    req(input$input_3d_er_interactive_ctrl_desc)
    req(input$input_3d_er_interactive_res_id)
    new_choices <- as_tibble(fread(file.path(config$project_dir, config$result_metadata_location, paste0(input$input_3d_er_interactive_res_id, ".csv")))) %>%
      filter(targ_description == input$input_3d_er_interactive_targ_desc & ctrl_description == input$input_3d_er_interactive_ctrl_desc)
    if(nrow(new_choices) > 0) {
      new_choices <- new_choices %>%
        select(ctrl_ct_lst, targ_ct_lst) %>%
        unlist() %>%
        str_split(", ") %>%
        unlist() %>%
        unique()
      new_choices <- sapply(strsplit(new_choices, "_"), `[`, 3) %>%
        str_remove_all("lib|\\.csv") %>%
        as.numeric() %>%
        unique()
      tibble(lib_id = new_choices) %>%
        left_join(rv$rv_lib_meta) %>%
        mutate(temp = paste0(lib_id, ": ", lib_name)) %>%
        pull(temp)
    }
  })
  
  observe({
    updateSelectInput(session, "input_3d_er_interactive_lib_id", choices = input_3d_er_interactive_lib_choices())
  })
  
  observe({
    updateSelectInput(session, "input_3d_er_interactive_res_id", "Choose res_id", rv$rv_res_lst)
  })
  
  observeEvent(input$visualize_3d_er_sp_interactive_load_cpd_struct, {
    temp_lib_id <- as.numeric(unlist(str_split(input$input_3d_er_interactive_lib_id, ": "))[1])
    if(length(rv$rv_3d_er_sp_interactive_struct_map) > 0) {
      showNotification("Metadata already loaded")
    }
    else if(length(rv$rv_lib_enum_map) > 0) {
      load_cpd_notif <- showNotification("Loading metadata...", duration = NULL)
      rv$rv_3d_er_sp_interactive_struct_map <- rv$rv_lib_enum_map[grep(paste0(temp_lib_id, " "), names(rv$rv_lib_enum_map))]
      removeNotification(load_cpd_notif)
      showNotification("Finished loading metadata")
    }
    else if(length(rv$rv_lib_enum_map) == 0) {
      load_cpd_notif <- showNotification("Loading metadata...", duration = NULL)
      temp <- as_tibble(fread(file.path(config$project_dir, config$enumerated_structure_metadata)))
      rv$rv_lib_enum_map <- structure(temp$structure, names = paste(temp$lib_id,
                                                                    temp$cycle1,
                                                                    temp$cycle2,
                                                                    temp$cycle3))
      rv$rv_3d_er_sp_interactive_struct_map <- rv$rv_lib_enum_map[grep(paste0(temp_lib_id, " "), names(rv$rv_lib_enum_map))]
      removeNotification(load_cpd_notif)
      showNotification("Finished loading metadata")
    }
  })
  
  observeEvent(input$input_3d_er_interactive_lib_id, {
    if(length(rv$rv_lib_enum_map) > 0) {
      temp_lib_id <- as.numeric(unlist(str_split(input$input_3d_er_interactive_lib_id, ": "))[1])
      load_cpd_notif <- showNotification("Loading metadata...", duration = NULL)
      rv$rv_3d_er_sp_interactive_struct_map <- rv$rv_lib_enum_map[grep(paste0(temp_lib_id, " "), names(rv$rv_lib_enum_map))]
      removeNotification(load_cpd_notif)
      showNotification("Finished loading metadata")
    }
  })
  
  observeEvent(input$visualize_3d_er_sp_interactive, {
    vis_3d_er_sp_interactive_notif <- showNotification("Rendering plot...", duration = NULL)
    
    temp_tbl <- as_tibble(fread(file.path(config$project_dir, config$result_metadata_location, paste0(input$input_3d_er_interactive_res_id, ".csv")))) %>%
      filter(targ_description == input$input_3d_er_interactive_targ_desc & ctrl_description == input$input_3d_er_interactive_ctrl_desc)
    
    temp_ctrl_lst <- temp_tbl %>% pull(ctrl_ct_lst) %>% str_split(", ") %>% unlist()
    temp_ctrl_ct <- as_tibble(fread(file.path(config$project_dir, config$data_count_location, temp_ctrl_lst[1])))
    if(length(temp_ctrl_lst) > 1) {
      for(i in 2:length(temp_ctrl_lst)) {
        temp_ctrl_ct <- as_tibble(temp_ctrl_ct + as_tibble(fread(file.path(config$project_dir, config$data_count_location, temp_ctrl_lst[i]))))
      }
    }
    temp_targ_lst <- temp_tbl %>% pull(targ_ct_lst) %>% str_split(", ") %>% unlist()
    temp_targ_ct <- as_tibble(fread(file.path(config$project_dir, config$data_count_location, temp_targ_lst[1])))
    if(length(temp_targ_lst) > 1) {
      for(i in 2:length(temp_targ_lst)) {
        temp_targ_ct <- as_tibble(temp_targ_ct + as_tibble(fread(file.path(config$project_dir, config$data_count_location, temp_targ_lst[i]))))
      }
    }
    temp_an_id <- temp_tbl %>%
      pull(an_id)
    temp_lib_id <- as.numeric(unlist(str_split(input$input_3d_er_interactive_lib_id, ": "))[1])
    temp <- as_tibble(fread(file.path(config$project_dir, config$library_metadata_column_location,
                                      sprintf("lib%03d_agg07_meta_cols.csv", temp_lib_id)))) %>%
      bind_cols(as_tibble(fread(file.path(config$project_dir, config$result_analysis_location,
                                          sprintf("an%06d_lib%03d_agg07_type02.csv", temp_an_id, temp_lib_id))))) %>%
      mutate(ctrl_ct = temp_ctrl_ct %>% pull(value),
             targ_ct = temp_targ_ct %>% pull(value))
    rv$rv_3d_er_sp_interactive_x_map <- structure(unique(temp$cycle1), names = unique(temp$cycle1))
    rv$rv_3d_er_sp_interactive_y_map <- structure(unique(temp$cycle2), names = unique(temp$cycle2))
    rv$rv_3d_er_sp_interactive_z_map <- structure(unique(temp$cycle3), names = unique(temp$cycle3))
    temp_labels <- temp[,order(colnames(temp))] %>%
      select(matches("cycle")) %>%
      apply(1, function(x) {paste0(x, collapse = ", ")})
    rv$rv_3d_er_sp_interactive_tbl <- temp %>%
      mutate(scaled_size = value*input$visualize_3d_er_sp_interactive_scale,
             label = paste0("(", cycle1, ", ", cycle2, ", ", cycle3 ,")\n", "er_lb = ", round(value, 2))) %>%
      arrange(-value)
    shinyjs::show("input_3d_er_sp_interactive_smi_to_clip")
    shinyjs::show("input_3d_er_sp_interactive_reorder_bbs")
    shinyjs::show("input_3d_er_sp_interactive_reset_order")
    shinyjs::show("output_3d_er_sp_interactive_save_plot")
    shinyjs::show("output_3d_er_sp_interactive_save_tbl")
    
    temp_scale <- max(25 / max(rv$rv_3d_er_sp_interactive_tbl$value, na.rm = T), 0.01)
    updateRadioButtons(session, "visualize_3d_er_sp_interactive_scale_log_radio", selected = "enrichment")
    updateSliderInput(session, "visualize_3d_er_sp_interactive_scale", value = temp_scale)
    removeNotification(vis_3d_er_sp_interactive_notif)
  })
  
  observe({
    req(rv$rv_3d_er_sp_interactive_tbl)
    input$visualize_3d_er_sp_interactive_scale
    input$visualize_3d_er_sp_interactive_scale_log_radio
    if(input$visualize_3d_er_sp_interactive_scale_log_radio == "enrichment") {
      rv$rv_3d_er_sp_interactive_tbl <- rv$rv_3d_er_sp_interactive_tbl %>%
        mutate(scaled_size = value*input$visualize_3d_er_sp_interactive_scale)
    }
    else if(input$visualize_3d_er_sp_interactive_scale_log_radio == "log(enrichment)") {
      rv$rv_3d_er_sp_interactive_tbl <- rv$rv_3d_er_sp_interactive_tbl %>%
        mutate(scaled_size = (value*input$visualize_3d_er_sp_interactive_scale) %>% log() * 10)
    }
  })
  
  output$output_3d_er_sp_interactive <- renderPlotly({
    req(rv$rv_3d_er_sp_interactive_tbl)
    plot_ly(rv$rv_3d_er_sp_interactive_tbl[1:input$visualize_3d_er_sp_interactive_num_points,] %>%
              mutate(x = rv$rv_3d_er_sp_interactive_x_map[as.character(cycle1)],
                     y = rv$rv_3d_er_sp_interactive_y_map[as.character(cycle2)],
                     z = rv$rv_3d_er_sp_interactive_z_map[as.character(cycle3)]), 
            x = ~x, y = ~y, z = ~z, text = ~label,
            hoverinfo = "text",
            # height = input$visualize_3d_er_sp_interactive_size,
            height = input$main_panel_innerWidth_3d_er_sp_interactive*2/5,
            marker = list(color = ~value, size = ~scaled_size, colorscale = c("#FFE1A1", "683531"), showscale = T),
            source = "plotly_3d_er_interactive") %>%
      add_markers() %>%
      layout(scene = list(xaxis = list(nticks = 10,
                                       range = c(1,length(unique(rv$rv_3d_er_sp_interactive_tbl %>% pull(cycle1)))),
                                       title = "cycle 1"),
                          yaxis = list(nticks = 10,
                                       range = c(1,length(unique(rv$rv_3d_er_sp_interactive_tbl %>% pull(cycle2)))),
                                       title = "cycle 2"),
                          zaxis = list(nticks = 10,
                                       range = c(1,length(unique(rv$rv_3d_er_sp_interactive_tbl %>% pull(cycle3)))),
                                       title = "cycle 3"),
                          aspectmode = "cube")) %>%
      event_register("plotly_hover") %>%
      event_register("plotly_unhover") %>%
      event_register("plotly_click")
  })
  
  plotly_3d_er_interactive_hover_event <- reactive({
    req(rv$rv_3d_er_sp_interactive_tbl)
    event_data(event = "plotly_hover", source = "plotly_3d_er_interactive")
  })
  
  plotly_3d_er_interactive_unhover_event <- reactive({
    req(rv$rv_3d_er_sp_interactive_tbl)
    event_data(event = "plotly_unhover", source = "plotly_3d_er_interactive")
  })
  
  plotly_3d_er_interactive_click_event <- reactive({
    req(rv$rv_3d_er_sp_interactive_tbl)
    event_data(event = "plotly_click", source = "plotly_3d_er_interactive")
  })
  
  observeEvent(plotly_3d_er_interactive_click_event(), {
    req(rv$rv_3d_er_sp_interactive_struct_map)
    temp_lib_id <- as.numeric(unlist(str_split(input$input_3d_er_interactive_lib_id, ": "))[1])
    cy1 <- which(rv$rv_3d_er_sp_interactive_x_map == plotly_3d_er_interactive_click_event()$x) %>% names()
    cy2 <- which(rv$rv_3d_er_sp_interactive_y_map == plotly_3d_er_interactive_click_event()$y) %>% names()
    cy3 <- which(rv$rv_3d_er_sp_interactive_z_map == plotly_3d_er_interactive_click_event()$z) %>% names()
    rv$rv_3d_er_sp_interactive_selected_struct <- rv$rv_3d_er_sp_interactive_struct_map[paste(temp_lib_id, cy1, cy2, cy3)]
  })
  
  output$output_3d_er_sp_interactive_selected_struct <- renderPlot({
    req(rv$rv_3d_er_sp_interactive_selected_struct)
    dep <- get.depictor(width=400, height=400)
    img <- parse.smiles(rv$rv_3d_er_sp_interactive_selected_struct)[[1]] %>%
      view.image.2d(dep)
    par(mar = c(0, 0, 0, 0))
    plot("", xlim = c(0, 1), ylim = c(0,1), xlab = "", ylab = "", xaxt = "n", yaxt = "n")
    rasterImage(img, 0, 0, 1, 1)
  },
  height = reactive(ifelse(!is.null(input$side_panel_innerWidth_3d_er_sp_interactive),input$side_panel_innerWidth_3d_er_sp_interactive*1/5,0)),
  width = reactive(ifelse(!is.null(input$side_panel_innerWidth_3d_er_sp_interactive),input$side_panel_innerWidth_3d_er_sp_interactive*1/5,0))
  )
  
  observeEvent(input$input_3d_er_sp_interactive_smi_to_clip, {
    req(rv$rv_3d_er_sp_interactive_selected_struct)
    clipr::write_clip(rv$rv_3d_er_sp_interactive_selected_struct)
    showNotification("SMILES copied to clipboard")
  })
  
  observeEvent(input$input_3d_er_sp_interactive_reorder_bbs, {
    temp_lib_id <- as.numeric(unlist(str_split(input$input_3d_er_interactive_lib_id, ": "))[1])
    cy1 <- which(rv$rv_3d_er_sp_interactive_x_map == plotly_3d_er_interactive_click_event()$x) %>% names() %>% as.numeric()
    cy2 <- which(rv$rv_3d_er_sp_interactive_y_map == plotly_3d_er_interactive_click_event()$y) %>% names() %>% as.numeric()
    cy3 <- which(rv$rv_3d_er_sp_interactive_z_map == plotly_3d_er_interactive_click_event()$z) %>% names() %>% as.numeric()
    cy1_mols <- bb_meta %>%
      filter(lib_id == temp_lib_id) %>%
      filter(cyc == 1) %>%
      pull(bb_smiles) %>%
      parse.smiles()
    cy1_mols[sapply(cy1_mols, is.null)] <- parse.smiles("C")
    fps <- lapply(cy1_mols, get.fingerprint, type = "circular")
    rv$rv_3d_er_sp_interactive_x_map <- structure(1:length(cy1_mols), 
                                                  names = fingerprint::fp.sim.matrix(fps[cy1], fps, method='tanimoto') %>% order())
    cy2_mols <- bb_meta %>%
      filter(lib_id == temp_lib_id) %>%
      filter(cyc == 2) %>%
      pull(bb_smiles) %>%
      parse.smiles()
    cy2_mols[sapply(cy2_mols, is.null)] <- parse.smiles("C")
    fps <- lapply(cy2_mols, get.fingerprint, type = "circular")
    rv$rv_3d_er_sp_interactive_y_map <- structure(1:length(cy2_mols), 
                                                  names = fingerprint::fp.sim.matrix(fps[cy2], fps, method='tanimoto') %>% order())
    cy3_mols <- bb_meta %>%
      filter(lib_id == temp_lib_id) %>%
      filter(cyc == 3) %>%
      pull(bb_smiles) %>%
      parse.smiles()
    cy3_mols[sapply(cy3_mols, is.null)] <- parse.smiles("C")
    fps <- lapply(cy3_mols, get.fingerprint, type = "circular")
    rv$rv_3d_er_sp_interactive_z_map <- structure(1:length(cy3_mols), 
                                                  names = fingerprint::fp.sim.matrix(fps[cy3], fps, method='tanimoto') %>% order())
  })
  
  observeEvent(input$input_3d_er_sp_interactive_reset_order, {
    rv$rv_3d_er_sp_interactive_x_map <- structure(1:length(rv$rv_3d_er_sp_interactive_x_map), names = 1:length(rv$rv_3d_er_sp_interactive_x_map))
    rv$rv_3d_er_sp_interactive_y_map <- structure(1:length(rv$rv_3d_er_sp_interactive_y_map), names = 1:length(rv$rv_3d_er_sp_interactive_y_map))
    rv$rv_3d_er_sp_interactive_z_map <- structure(1:length(rv$rv_3d_er_sp_interactive_z_map), names = 1:length(rv$rv_3d_er_sp_interactive_z_map))
  })
  
  output$output_3d_er_sp_interactive_selection <- renderText({
    cy1 <- which(rv$rv_3d_er_sp_interactive_x_map == plotly_3d_er_interactive_click_event()$x) %>% names() %>% as.numeric()
    cy2 <- which(rv$rv_3d_er_sp_interactive_y_map == plotly_3d_er_interactive_click_event()$y) %>% names() %>% as.numeric()
    cy3 <- which(rv$rv_3d_er_sp_interactive_z_map == plotly_3d_er_interactive_click_event()$z) %>% names() %>% as.numeric()
    if(length(cy1) == 1 & length(cy2) == 1 & length(cy3) == 1) {
      temp <- rv$rv_3d_er_sp_interactive_tbl %>%
        filter(cycle1 == cy1 & cycle2 == cy2 & cycle3 == cy3)
      paste0("Selected trisynthon: (", paste(c(cy1, cy2, cy3), collapse = ", "), ")\n",
                         "Control counts: ", temp %>% pull(ctrl_ct), " (", 
                         (temp %>% pull(ctrl_ct) / sum(rv$rv_3d_er_sp_interactive_tbl$ctrl_ct) * 1E6) %>% round(3), " ppm)\n",
                         "Target counts: ", temp %>% pull(targ_ct), " (", 
                         (temp %>% pull(targ_ct) / sum(rv$rv_3d_er_sp_interactive_tbl$targ_ct) * 1E6) %>% round(3), " ppm)")
    }
  })
  
  output$output_3d_er_sp_interactive_selected_er <- renderPlotly({
    req(rv$rv_3d_er_sp_interactive_tbl)
    cy1 <- which(rv$rv_3d_er_sp_interactive_x_map == plotly_3d_er_interactive_click_event()$x) %>% names() %>% as.numeric()
    cy2 <- which(rv$rv_3d_er_sp_interactive_y_map == plotly_3d_er_interactive_click_event()$y) %>% names() %>% as.numeric()
    cy3 <- which(rv$rv_3d_er_sp_interactive_z_map == plotly_3d_er_interactive_click_event()$z) %>% names() %>% as.numeric()
    if(length(cy1) == 1 & length(cy2) == 1 & length(cy3) == 1) {
      temp <- rv$rv_3d_er_sp_interactive_tbl %>%
        filter(cycle1 == cy1 & cycle2 == cy2 & cycle3 == cy3)
      k1 <- temp %>% pull(targ_ct) %>% as.numeric()
      k2 <- temp %>% pull(ctrl_ct) %>% as.numeric()
      n1 <- sum(rv$rv_3d_er_sp_interactive_tbl$targ_ct)
      n2 <- sum(rv$rv_3d_er_sp_interactive_tbl$ctrl_ct)
      er_val <- temp %>% pull(value)
      temp_tbl <- tibble(x = seq(0.1, er_val*5, length.out = 5000),
                         y = z_func_prob(k1,k2,n1,n2,x)*2.5) %>%
        filter(y > 0.00001)
      if(nrow(temp_tbl) > 1000) {
        every_n <- floor(nrow(temp_tbl) / 500)
        temp_tbl <- temp_tbl[seq(1, nrow(temp_tbl), by = every_n),]
      }
      plot_ly(temp_tbl, x = ~x, y = ~y,
              type = "scattergl", mode = "lines", height = 250) %>%
      layout(xaxis = list(title = "enrichment ratio",
                          showspikes = T,
                          spikemode = "across+toaxis",
                          spikesnap = "cursor",
                          spikedash = "solid",
                          showline = T),
             yaxis = list(title = "relative probability",
                          showticklabels = F),
             hovermode = "x",
             spikedistance = -1)
    }
  })
  
  output$output_3d_er_sp_interactive_save_plot <- downloadHandler("plot_name.html", 
                                                                  content = function(file) {
                                                                    req(rv$rv_3d_er_sp_interactive_tbl)
                                                                    temp <- plot_ly(rv$rv_3d_er_sp_interactive_tbl[1:input$visualize_3d_er_sp_interactive_num_points,] %>%
                                                                                      mutate(x = rv$rv_3d_er_sp_interactive_x_map[as.character(cycle1)],
                                                                                             y = rv$rv_3d_er_sp_interactive_y_map[as.character(cycle2)],
                                                                                             z = rv$rv_3d_er_sp_interactive_z_map[as.character(cycle3)]), 
                                                                                    x = ~x, y = ~y, z = ~z, text = ~label,
                                                                                    hoverinfo = "text",
                                                                                    marker = list(color = ~value, size = ~scaled_size, colorscale = c("#FFE1A1", "683531"), showscale = T)) %>%
                                                                      add_markers() %>%
                                                                      layout(scene = list(xaxis = list(nticks = 10,
                                                                                                       range = c(1,length(unique(rv$rv_3d_er_sp_interactive_tbl %>% pull(cycle1)))),
                                                                                                       title = "Axis 1"),
                                                                                          yaxis = list(nticks = 10,
                                                                                                       range = c(1,length(unique(rv$rv_3d_er_sp_interactive_tbl %>% pull(cycle2)))),
                                                                                                       title = "Axis 2"),
                                                                                          zaxis = list(nticks = 10,
                                                                                                       range = c(1,length(unique(rv$rv_3d_er_sp_interactive_tbl %>% pull(cycle3)))),
                                                                                                       title = "Axis 3"),
                                                                                          aspectmode = "cube"))
                                                                    htmlwidgets::saveWidget(partial_bundle(temp), file)
                                                                  },
                                                                  contentType = "html")
  
  output$output_3d_er_sp_interactive_save_tbl <- downloadHandler("table_name.csv", 
                                                                 content = function(file) {
                                                                   req(rv$rv_3d_er_sp_interactive_tbl)
                                                                   rv$rv_3d_er_sp_interactive_tbl[1:input$visualize_3d_er_sp_interactive_num_points,] %>%
                                                                     mutate(x = rv$rv_3d_er_sp_interactive_x_map[as.character(cycle1)],
                                                                            y = rv$rv_3d_er_sp_interactive_y_map[as.character(cycle2)],
                                                                            z = rv$rv_3d_er_sp_interactive_z_map[as.character(cycle3)]) %>%
                                                                     rename(er_lb = value) %>%
                                                                     select(-scaled_size, -x, -y, -z, -label) %>%
                                                                     write.csv(file, row.names = F)
                                                                 },
                                                                 contentType = "text/csv")
  
  
  
  
  
  
  

  input_interactive_y_samp_choices <- reactive({
    req(input$input_interactive_y_res_id)
    as_tibble(fread(file.path(config$project_dir, config$result_metadata_location, paste0(input$input_interactive_y_res_id, ".csv")))) %>%
      filter(grepl(sprintf("lib%03d", as.numeric(unlist(str_split(input$input_interactive_lib_id, ": "))[1])), targ_ct_lst)) %>%
      select(targ_description, ctrl_description) %>%
      unlist() %>%
      unique()
  })
  
  observe({
    updateSelectInput(session, "input_interactive_y_samp_desc", choices = input_interactive_y_samp_choices())
  })
  
  input_interactive_y_targ_choices <- reactive({
    req(input$input_interactive_y_res_id)
    as_tibble(fread(file.path(config$project_dir, config$result_metadata_location, paste0(input$input_interactive_y_res_id, ".csv")))) %>%
      filter(grepl(sprintf("lib%03d", as.numeric(unlist(str_split(input$input_interactive_lib_id, ": "))[1])), targ_ct_lst)) %>%
      pull(targ_description)
  })
  
  observe({
    updateSelectInput(session, "input_interactive_y_targ_desc", choices = input_interactive_y_targ_choices())
  })
  
  input_interactive_y_ctrl_choices <- reactive({
    req(input$input_interactive_y_res_id)
    as_tibble(fread(file.path(config$project_dir, config$result_metadata_location, paste0(input$input_interactive_y_res_id, ".csv")))) %>%
      filter(grepl(sprintf("lib%03d", as.numeric(unlist(str_split(input$input_interactive_lib_id, ": "))[1])), ctrl_ct_lst)) %>%
      filter(targ_description == input$input_interactive_y_targ_desc) %>%
      pull(ctrl_description)
  })
  
  observe({
    updateSelectInput(session, "input_interactive_y_ctrl_desc", choices = input_interactive_y_ctrl_choices())
  })
  
  input_interactive_x_samp_choices <- reactive({
    req(input$input_interactive_x_res_id)
    as_tibble(fread(file.path(config$project_dir, config$result_metadata_location, paste0(input$input_interactive_x_res_id, ".csv")))) %>%
      filter(grepl(sprintf("lib%03d", as.numeric(unlist(str_split(input$input_interactive_lib_id, ": "))[1])), targ_ct_lst)) %>%
      select(targ_description, ctrl_description) %>%
      unlist() %>%
      unique()
  })
  
  observe({
    updateSelectInput(session, "input_interactive_x_samp_desc", choices = input_interactive_x_samp_choices())
  })
  
  input_interactive_x_targ_choices <- reactive({
    req(input$input_interactive_x_res_id)
    as_tibble(fread(file.path(config$project_dir, config$result_metadata_location, paste0(input$input_interactive_x_res_id, ".csv")))) %>%
      filter(grepl(sprintf("lib%03d", as.numeric(unlist(str_split(input$input_interactive_lib_id, ": "))[1])), targ_ct_lst)) %>%
      pull(targ_description)
  })
  
  observe({
    updateSelectInput(session, "input_interactive_x_targ_desc", choices = input_interactive_x_targ_choices())
  })
  
  input_interactive_x_ctrl_choices <- reactive({
    req(input$input_interactive_x_res_id)
    as_tibble(fread(file.path(config$project_dir, config$result_metadata_location, paste0(input$input_interactive_x_res_id, ".csv")))) %>%
      filter(grepl(sprintf("lib%03d", as.numeric(unlist(str_split(input$input_interactive_lib_id, ": "))[1])), ctrl_ct_lst)) %>%
      filter(targ_description == input$input_interactive_x_targ_desc) %>%
      pull(ctrl_description)
  })
  
  observe({
    updateSelectInput(session, "input_interactive_x_ctrl_desc", choices = input_interactive_x_ctrl_choices())
  })
  
  observe({
    req(input$input_interactive_lib_id)
    req(input$input_interactive_agg_id)
    rv$rv_interactive_meta <- as_tibble(fread(file.path(config$project_dir, config$library_metadata_column_location,
                                                        sprintf("lib%03d_agg%02d_meta_cols.csv", as.numeric(unlist(str_split(input$input_interactive_lib_id, ": "))[1]),
                                                                as.numeric(unlist(str_split(input$input_interactive_agg_id, ": "))[1])))))
  })
  
  observe({
    req(input$input_interactive_agg_id)
    req(input$input_interactive_y_samp_desc)
    req(input$input_interactive_y_type)
    
    temp <- as_tibble(fread(file.path(config$project_dir, config$result_metadata_location, paste0(input$input_interactive_y_res_id, ".csv"))))
    temp <- tibble(ct_filenames = c(temp$targ_ct_lst, temp$ctrl_ct_lst),
                   description = c(temp$targ_description, temp$ctrl_description)) %>%
      unique() %>%
      filter(description == input$input_interactive_y_samp_desc)
    if(nrow(temp) > 0) {
      temp <- temp %>%
        pull(ct_filenames) %>%
        str_split(", ") %>%
        unlist()
      temp_lib_id <- unlist(str_split(temp[1], "_"))[3] %>% str_remove_all("lib|.csv") %>% as.numeric()
      y_vals <- as_tibble(fread(file.path(config$project_dir, config$data_count_location, temp[1])))
      if(length(temp) > 1) {
        for(i in 2:length(temp)) {
          y_vals <- as_tibble(y_vals + as_tibble(fread(file.path(config$project_dir, config$data_count_location, temp[i]))))
        }
      }
      if(nrow(rv$rv_interactive_meta) < nrow(y_vals)) {
        y_vals_agg <- as_tibble(fread(file.path(config$project_dir, config$library_metadata_column_location,
                                                sprintf("lib%03d_agg%02d_meta_cols.csv", temp_lib_id, 7)))) %>%
          bind_cols(y_vals)
        cycs <- paste0("cycle", unlist(str_split(agg_meta %>% filter(agg_id == as.numeric(unlist(str_split(input$input_interactive_agg_id, ": "))[1])) %>%
                                                   pull(cyc), ", ")))
        y_vals_agg <- y_vals_agg %>%
          group_by_at(c("lib_id", cycs)) %>%
          summarise(value = sum(value),
                    .groups = "keep") %>%
          pull(value)
        
        rv$rv_interactive_y_samp_ct <- y_vals_agg
        
        if(input$input_interactive_y_type == "count") {
          rv$rv_interactive_y_vals <- y_vals_agg
        }
        else if(input$input_interactive_y_type == "rank by counts") {
          rv$rv_interactive_y_vals <- rank(y_vals_agg)
        }
      }
      else {
        rv$rv_interactive_y_samp_ct <- y_vals %>% pull(value)
        
        if(input$input_interactive_y_type == "count") {
          rv$rv_interactive_y_vals <- y_vals %>% pull(value)
        }
        else if(input$input_interactive_y_type == "rank by counts") {
          rv$rv_interactive_y_vals <- y_vals %>% pull(value) %>% rank()
        }
      }
    }
  })
  
  
  observe({
    input$input_interactive_lib_id
    input$input_interactive_x_res_id
    input$input_interactive_y_res_id
    updateRadioButtons(session, "input_interactive_vis", selected = "Off")
  })
  
  
  observe({
    req(input$input_interactive_agg_id)
    req(input$input_interactive_y_targ_desc)
    req(input$input_interactive_y_ctrl_desc)
    req(rv$rv_interactive_meta)
    
    if(input$input_interactive_y_type == "enrichment") {
      temp <- as_tibble(fread(file.path(config$project_dir, config$result_metadata_location, paste0(input$input_interactive_y_res_id, ".csv"))))
      temp <- tibble(ct_filenames = c(temp$targ_ct_lst, temp$ctrl_ct_lst),
                     description = c(temp$targ_description, temp$ctrl_description)) %>%
        unique() %>%
        filter(description == input$input_interactive_y_targ_desc)
      if(nrow(temp) > 0) {
        temp <- temp %>%
          pull(ct_filenames) %>%
          str_split(", ") %>%
          unlist()
        temp_lib_id <- unlist(str_split(temp[1], "_"))[3] %>% str_remove_all("lib|.csv") %>% as.numeric()
        y_vals <- as_tibble(fread(file.path(config$project_dir, config$data_count_location, temp[1])))
        if(length(temp) > 1) {
          for(i in 2:length(temp)) {
            y_vals <- as_tibble(y_vals + as_tibble(fread(file.path(config$project_dir, config$data_count_location, temp[i]))))
          }
        }
        if(nrow(rv$rv_interactive_meta) < nrow(y_vals)) {
          y_vals_agg <- as_tibble(fread(file.path(config$project_dir, config$library_metadata_column_location,
                                                  sprintf("lib%03d_agg%02d_meta_cols.csv", temp_lib_id, 7)))) %>%
            bind_cols(y_vals)
          cycs <- paste0("cycle", unlist(str_split(agg_meta %>% filter(agg_id == as.numeric(unlist(str_split(input$input_interactive_agg_id, ": "))[1])) %>%
                                                     pull(cyc), ", ")))
          y_vals_agg <- y_vals_agg %>%
            group_by_at(c("lib_id", cycs)) %>%
            summarise(value = sum(value),
                      .groups = "keep") %>%
            pull(value)
          
          rv$rv_interactive_y_targ_ct <- y_vals_agg
        }
        else {
          rv$rv_interactive_y_targ_ct <- y_vals %>% pull(value)
        }
      }
      
      temp <- as_tibble(fread(file.path(config$project_dir, config$result_metadata_location, paste0(input$input_interactive_y_res_id, ".csv"))))
      temp <- tibble(ct_filenames = c(temp$targ_ct_lst, temp$ctrl_ct_lst),
                     description = c(temp$targ_description, temp$ctrl_description)) %>%
        unique() %>%
        filter(description == input$input_interactive_y_ctrl_desc)
      if(nrow(temp) > 0) {
        temp <- temp %>%
          pull(ct_filenames) %>%
          str_split(", ") %>%
          unlist()
        temp_lib_id <- unlist(str_split(temp[1], "_"))[3] %>% str_remove_all("lib|.csv") %>% as.numeric()
        y_vals <- as_tibble(fread(file.path(config$project_dir, config$data_count_location, temp[1])))
        if(length(temp) > 1) {
          for(i in 2:length(temp)) {
            y_vals <- as_tibble(y_vals + as_tibble(fread(file.path(config$project_dir, config$data_count_location, temp[i]))))
          }
        }
        if(nrow(rv$rv_interactive_meta) < nrow(y_vals)) {
          y_vals_agg <- as_tibble(fread(file.path(config$project_dir, config$library_metadata_column_location,
                                                  sprintf("lib%03d_agg%02d_meta_cols.csv", temp_lib_id, 7)))) %>%
            bind_cols(y_vals)
          cycs <- paste0("cycle", unlist(str_split(agg_meta %>% filter(agg_id == as.numeric(unlist(str_split(input$input_interactive_agg_id, ": "))[1])) %>%
                                                     pull(cyc), ", ")))
          y_vals_agg <- y_vals_agg %>%
            group_by_at(c("lib_id", cycs)) %>%
            summarise(value = sum(value),
                      .groups = "keep") %>%
            pull(value)
          
          rv$rv_interactive_y_ctrl_ct <- y_vals_agg
        }
        else {
          rv$rv_interactive_y_ctrl_ct <- y_vals %>% pull(value)
        }
      }
      
      temp <- as_tibble(fread(file.path(config$project_dir, config$result_metadata_location, paste0(input$input_interactive_y_res_id, ".csv")))) %>%
        filter(targ_description == input$input_interactive_y_targ_desc & ctrl_description == input$input_interactive_y_ctrl_desc)
      if(nrow(temp) > 0) {
        temp_an_id <- temp %>% pull(an_id)
        temp_lib_id <- temp %>% pull(targ_ct_lst) %>% str_remove_all(".csv")
        temp_lib_id <- unlist(str_split(temp_lib_id, ", "))[1]
        temp_lib_id <- unlist(str_split(temp_lib_id, "_"))[3] %>% str_remove("lib") %>% as.numeric()
        y_vals <- as_tibble(fread(file.path(config$project_dir, config$result_analysis_location,
                                            sprintf("an%06d_lib%03d_agg%02d_type02.csv", 
                                                    temp_an_id, 
                                                    temp_lib_id,
                                                    as.numeric(unlist(str_split(input$input_interactive_agg_id, ": "))[1])))))
        rv$rv_interactive_y_vals <- y_vals %>% pull(value)
      }
    }
  })
  
  observe({
    req(input$input_interactive_agg_id)
    req(input$input_interactive_x_samp_desc)
    req(input$input_interactive_x_type)
    temp <- as_tibble(fread(file.path(config$project_dir, config$result_metadata_location, paste0(input$input_interactive_x_res_id, ".csv"))))
    temp <- tibble(ct_filenames = c(temp$targ_ct_lst, temp$ctrl_ct_lst),
                   description = c(temp$targ_description, temp$ctrl_description)) %>%
      unique() %>%
      filter(description == input$input_interactive_x_samp_desc)
    if(nrow(temp) > 0) {
      temp <- temp %>%
        pull(ct_filenames) %>%
        str_split(", ") %>%
        unlist()
      temp_lib_id <- unlist(str_split(temp[1], "_"))[3] %>% str_remove_all("lib|.csv") %>% as.numeric()
      x_vals <- as_tibble(fread(file.path(config$project_dir, config$data_count_location, temp[1])))
      if(length(temp) > 1) {
        for(i in 2:length(temp)) {
          x_vals <- as_tibble(x_vals + as_tibble(fread(file.path(config$project_dir, config$data_count_location, temp[i]))))
        }
      }
      if(nrow(rv$rv_interactive_meta) < nrow(x_vals)) {
        x_vals_agg <- as_tibble(fread(file.path(config$project_dir, config$library_metadata_column_location,
                                                sprintf("lib%03d_agg%02d_meta_cols.csv", temp_lib_id, 7)))) %>%
          bind_cols(x_vals)
        cycs <- paste0("cycle", unlist(str_split(agg_meta %>% filter(agg_id == as.numeric(unlist(str_split(input$input_interactive_agg_id, ": "))[1])) %>%
                                                   pull(cyc), ", ")))
        x_vals_agg <- x_vals_agg %>%
          group_by_at(c("lib_id", cycs)) %>%
          summarise(value = sum(value),
                    .groups = "keep") %>%
          pull(value)
        
        rv$rv_interactive_x_samp_ct <- x_vals_agg
        
        if(input$input_interactive_x_type == "count") {
          rv$rv_interactive_x_vals <- x_vals_agg
        }
        else if(input$input_interactive_x_type == "rank by counts") {
          rv$rv_interactive_x_vals <- rank(x_vals_agg)
        }
      }
      else {
        rv$rv_interactive_x_samp_ct <- x_vals %>% pull(value)
        
        if(input$input_interactive_x_type == "count") {
          rv$rv_interactive_x_vals <- x_vals %>% pull(value)
        }
        else if(input$input_interactive_x_type == "rank by counts") {
          rv$rv_interactive_x_vals <- x_vals %>% pull(value) %>% rank()
        }
      }
    }
  })
  
  observe({
    req(input$input_interactive_agg_id)
    req(input$input_interactive_x_targ_desc)
    req(input$input_interactive_x_ctrl_desc)
    req(rv$rv_interactive_meta)
    
    if(input$input_interactive_x_type == "enrichment") {
      temp <- as_tibble(fread(file.path(config$project_dir, config$result_metadata_location, paste0(input$input_interactive_x_res_id, ".csv"))))
      temp <- tibble(ct_filenames = c(temp$targ_ct_lst, temp$ctrl_ct_lst),
                     description = c(temp$targ_description, temp$ctrl_description)) %>%
        unique() %>%
        filter(description == input$input_interactive_x_targ_desc)
      if(nrow(temp) > 0) {
        temp <- temp %>%
          pull(ct_filenames) %>%
          str_split(", ") %>%
          unlist()
        temp_lib_id <- unlist(str_split(temp[1], "_"))[3] %>% str_remove_all("lib|.csv") %>% as.numeric()
        x_vals <- as_tibble(fread(file.path(config$project_dir, config$data_count_location, temp[1])))
        if(length(temp) > 1) {
          for(i in 2:length(temp)) {
            x_vals <- as_tibble(x_vals + as_tibble(fread(file.path(config$project_dir, config$data_count_location, temp[i]))))
          }
        }
        if(nrow(rv$rv_interactive_meta) < nrow(x_vals)) {
          x_vals_agg <- as_tibble(fread(file.path(config$project_dir, config$library_metadata_column_location,
                                                  sprintf("lib%03d_agg%02d_meta_cols.csv", temp_lib_id, 7)))) %>%
            bind_cols(x_vals)
          cycs <- paste0("cycle", unlist(str_split(agg_meta %>% filter(agg_id == as.numeric(unlist(str_split(input$input_interactive_agg_id, ": "))[1])) %>%
                                                     pull(cyc), ", ")))
          x_vals_agg <- x_vals_agg %>%
            group_by_at(c("lib_id", cycs)) %>%
            summarise(value = sum(value),
                      .groups = "keep") %>%
            pull(value)
          
          rv$rv_interactive_x_targ_ct <- x_vals_agg
        }
        else {
          rv$rv_interactive_x_targ_ct <- x_vals %>% pull(value)
        }
      }
      
      temp <- as_tibble(fread(file.path(config$project_dir, config$result_metadata_location, paste0(input$input_interactive_x_res_id, ".csv"))))
      temp <- tibble(ct_filenames = c(temp$targ_ct_lst, temp$ctrl_ct_lst),
                     description = c(temp$targ_description, temp$ctrl_description)) %>%
        unique() %>%
        filter(description == input$input_interactive_x_ctrl_desc)
      if(nrow(temp) > 0) {
        temp <- temp %>%
          pull(ct_filenames) %>%
          str_split(", ") %>%
          unlist()
        temp_lib_id <- unlist(str_split(temp[1], "_"))[3] %>% str_remove_all("lib|.csv") %>% as.numeric()
        x_vals <- as_tibble(fread(file.path(config$project_dir, config$data_count_location, temp[1])))
        if(length(temp) > 1) {
          for(i in 2:length(temp)) {
            x_vals <- as_tibble(x_vals + as_tibble(fread(file.path(config$project_dir, config$data_count_location, temp[i]))))
          }
        }
        if(nrow(rv$rv_interactive_meta) < nrow(x_vals)) {
          x_vals_agg <- as_tibble(fread(file.path(config$project_dir, config$library_metadata_column_location,
                                                  sprintf("lib%03d_agg%02d_meta_cols.csv", temp_lib_id, 7)))) %>%
            bind_cols(x_vals)
          cycs <- paste0("cycle", unlist(str_split(agg_meta %>% filter(agg_id == as.numeric(unlist(str_split(input$input_interactive_agg_id, ": "))[1])) %>%
                                                     pull(cyc), ", ")))
          x_vals_agg <- x_vals_agg %>%
            group_by_at(c("lib_id", cycs)) %>%
            summarise(value = sum(value),
                      .groups = "keep") %>%
            pull(value)
          
          rv$rv_interactive_x_ctrl_ct <- x_vals_agg
        }
        else {
          rv$rv_interactive_x_ctrl_ct <- x_vals %>% pull(value)
        }
      }
      
      temp <- as_tibble(fread(file.path(config$project_dir, config$result_metadata_location, paste0(input$input_interactive_x_res_id, ".csv")))) %>%
        filter(targ_description == input$input_interactive_x_targ_desc & ctrl_description == input$input_interactive_x_ctrl_desc)
      if(nrow(temp) > 0) {
        temp_an_id <- temp %>% pull(an_id)
        temp_lib_id <- temp %>% pull(targ_ct_lst) %>% str_remove_all(".csv")
        temp_lib_id <- unlist(str_split(temp_lib_id, ", "))[1]
        temp_lib_id <- unlist(str_split(temp_lib_id, "_"))[3] %>% str_remove("lib") %>% as.numeric()
        x_vals <- as_tibble(fread(file.path(config$project_dir, config$result_analysis_location,
                                            sprintf("an%06d_lib%03d_agg%02d_type02.csv", 
                                                    temp_an_id,
                                                    temp_lib_id,
                                                    as.numeric(unlist(str_split(input$input_interactive_agg_id, ": "))[1])))))
        rv$rv_interactive_x_vals <- x_vals %>% pull(value)
      }
    }
  })
  
  observeEvent(input$input_interactive_agg_id, {
    cycs <- paste0("cycle", unlist(str_split(agg_meta %>% filter(agg_id == as.numeric(unlist(str_split(input$input_interactive_agg_id, ": "))[1])) %>%
                                               pull(cyc), ", ")))
    if("cycle1" %in% cycs) {
      enable("input_interactive_color_cy1")
    }
    else {
      disable("input_interactive_color_cy1")
    }
    if("cycle2" %in% cycs) {
      enable("input_interactive_color_cy2")
    }
    else {
      disable("input_interactive_color_cy2")
    }
    if("cycle3" %in% cycs) {
      enable("input_interactive_color_cy3")
    }
    else {
      disable("input_interactive_color_cy3")
    }
  })
  
  observeEvent(input$input_interactive_color_tbl, {
    if(input$input_interactive_color_tbl == "On") {
      updateRadioButtons(session, "input_interactive_color", selected = "Off")
    }
  })
  
  observeEvent(input$input_interactive_color, {
    if(input$input_interactive_color == "On") {
      updateRadioButtons(session, "input_interactive_color_tbl", selected = "Off")
    }
  })
  
  observe({
    req(input$input_interactive_color)
    req(input$input_interactive_color_tbl)
    req(rv$rv_interactive_meta)
    temp_cols <- rep("Other", nrow(rv$rv_interactive_meta))
    if(input$input_interactive_color == "On") {
      cy1 <- input$input_interactive_color_cy1 %>% str_split(";") %>% unlist() %>% lapply(get_range) %>% unlist() %>% as.numeric()
      cy2 <- input$input_interactive_color_cy2 %>% str_split(";") %>% unlist() %>% lapply(get_range) %>% unlist() %>% as.numeric()
      cy3 <- input$input_interactive_color_cy3 %>% str_split(";") %>% unlist() %>% lapply(get_range) %>% unlist() %>% as.numeric()
      curr_selection <- rep(T, nrow(rv$rv_interactive_meta))
      if(!is.na(cy1)) {
        curr_selection[!rv$rv_interactive_meta$cycle1 %in% cy1] <- F
      }
      if(!is.na(cy2)) {
        curr_selection[!rv$rv_interactive_meta$cycle2 %in% cy2] <- F
      }
      if(!is.na(cy3)) {
        curr_selection[!rv$rv_interactive_meta$cycle3 %in% cy3] <- F
      }
      temp_cols[curr_selection] <- "Selection"
    }
    else if(input$input_interactive_color_tbl == "On") {
      req(input$input_interactive_color_tbl_request)
      temp_tbl <- as_tibble(fread(input$input_interactive_color_tbl_request$datapath)) 
      if(length(setdiff(colnames(temp_tbl), colnames(rv$rv_interactive_meta))) == 0) {
        temp_tbl <- temp_tbl %>%
          mutate(color = "Selection")
        temp_cols <- rv$rv_interactive_meta %>%
          left_join(temp_tbl) %>%
          pull(color)
        temp_cols[is.na(temp_cols)] <- "Other"
      }
    }
    rv$rv_interactive_col_vals <- temp_cols
  })

  observeEvent(input$update_interactive_res_id, {
    updateSelectInput(session, "input_interactive_y_res_id", "Choose res_id", rv$rv_res_lst)
    updateSelectInput(session, "input_interactive_x_res_id", "Choose res_id", rv$rv_res_lst)
  })
  
  observe({
    req(rv$rv_interactive_x_vals)
    req(rv$rv_interactive_y_vals)
    req(rv$rv_interactive_meta)
    req(isTruthy(rv$rv_interactive_x_samp_ct) || (isTruthy(rv$rv_interactive_x_ctrl_ct) & isTruthy(rv$rv_interactive_x_targ_ct)))
    req(isTruthy(rv$rv_interactive_y_samp_ct) || (isTruthy(rv$rv_interactive_y_ctrl_ct) & isTruthy(rv$rv_interactive_y_targ_ct)))
    req(input$input_interactive_max_points)
    
    if(nrow(rv$rv_interactive_meta) == length(rv$rv_interactive_x_vals) &
       nrow(rv$rv_interactive_meta) == length(rv$rv_interactive_y_vals) &
       ((input$input_interactive_x_type == "enrichment" & 
        nrow(rv$rv_interactive_meta) == length(rv$rv_interactive_x_ctrl_ct) &
        nrow(rv$rv_interactive_meta) == length(rv$rv_interactive_x_targ_ct)) |
        (input$input_interactive_x_type %in% c("count", "rank by counts") & 
         nrow(rv$rv_interactive_meta) == length(rv$rv_interactive_x_samp_ct))) &
       ((input$input_interactive_y_type == "enrichment" & 
        nrow(rv$rv_interactive_meta) == length(rv$rv_interactive_y_ctrl_ct) &
        nrow(rv$rv_interactive_meta) == length(rv$rv_interactive_y_targ_ct)) | 
        (input$input_interactive_y_type %in% c("count", "rank by counts") & 
         nrow(rv$rv_interactive_meta) == length(rv$rv_interactive_y_samp_ct)))) {
      temp <- rv$rv_interactive_meta %>%
        mutate(x_value = rv$rv_interactive_x_vals,
               y_value = rv$rv_interactive_y_vals)
      
      if(input$input_interactive_x_type == "enrichment") {
        temp <- temp %>%
          mutate(x_ctrl_ct = rv$rv_interactive_x_ctrl_ct,
                 x_targ_ct = rv$rv_interactive_x_targ_ct,
                 x_ctrl_ct_ppm = round(x_ctrl_ct / sum(x_ctrl_ct) * 1E6, 3),
                 x_targ_ct_ppm = round(x_targ_ct / sum(x_targ_ct) * 1E6, 3))
      }
      else {
        temp <- temp %>%
          mutate(x_samp_ct = rv$rv_interactive_x_samp_ct,
                 x_samp_ct_ppm = round(x_samp_ct / sum(x_samp_ct) * 1E6, 3))
      }
      
      if(input$input_interactive_y_type == "enrichment") {
        temp <- temp %>%
          mutate(y_ctrl_ct = rv$rv_interactive_y_ctrl_ct,
                 y_targ_ct = rv$rv_interactive_y_targ_ct,
                 y_ctrl_ct_ppm = round(y_ctrl_ct / sum(y_ctrl_ct) * 1E6, 3),
                 y_targ_ct_ppm = round(y_targ_ct / sum(y_targ_ct) * 1E6, 3))
      }
      else {
        temp <- temp %>%
          mutate(y_samp_ct = rv$rv_interactive_y_samp_ct,
                 y_samp_ct_ppm = round(y_samp_ct / sum(y_samp_ct) * 1E6, 3))
      }
      
      if(nrow(temp) == length(rv$rv_interactive_col_vals)) {
        temp <- temp %>%
          mutate(color = rv$rv_interactive_col_vals) %>%
          top_n(input$input_interactive_max_points, y_value)
        temp_agg_id <- as.numeric(unlist(str_split(input$input_interactive_agg_id, ": "))[1])
        missing_cycs <- paste0("cycle", unlist(str_split(agg_meta %>% filter(agg_id == max(agg_id)) %>% pull(cyc), ", "))) %>%
          setdiff(colnames(temp))
        for(cyc in missing_cycs) {
          temp[cyc] <- "?"
        }
        temp_labels <- temp[,order(colnames(temp))] %>%
          select(matches("cycle")) %>%
          apply(1, function(x) {paste0(x, collapse = ", ")})
        temp_labels <- paste0("(", temp_labels, ")\n\ncounts")
        temp_samp_lst <- c()
        if(input$input_interactive_x_type == "enrichment") {
          temp_labels <- paste0(temp_labels, "\n",
                                input$input_interactive_x_ctrl_desc, ": ", temp$x_ctrl_ct, "\n",
                                input$input_interactive_x_targ_desc, ": ", temp$x_targ_ct)
          temp_samp_lst <- c(temp_samp_lst, input$input_interactive_x_ctrl_desc, input$input_interactive_x_targ_desc)
        }
        else {
          temp_labels <- paste0(temp_labels, "\n",
                                input$input_interactive_x_samp_desc, ": ", temp$x_samp_ct)
          temp_samp_lst <- c(temp_samp_lst, input$input_interactive_x_samp_desc)
        }
        
        if(input$input_interactive_y_type == "enrichment") {
          if(!input$input_interactive_y_ctrl_desc %in% temp_samp_lst) {
            temp_labels <- paste0(temp_labels, "\n",
                                  input$input_interactive_y_ctrl_desc, ": ", temp$y_ctrl_ct)
            temp_samp_lst <- c(temp_samp_lst, input$input_interactive_y_ctrl_desc)
          }
          if(!input$input_interactive_y_targ_desc %in% temp_samp_lst) {
            temp_labels <- paste0(temp_labels, "\n",
                                  input$input_interactive_y_targ_desc, ": ", temp$y_targ_ct)
            temp_samp_lst <- c(temp_samp_lst, input$input_interactive_y_targ_desc)
          }
        }
        else {
          if(!input$input_interactive_y_samp_desc %in% temp_samp_lst) {
            temp_labels <- paste0(temp_labels, "\n",
                                  input$input_interactive_y_samp_desc, ": ", temp$y_samp_ct)
            temp_samp_lst <- c(temp_samp_lst, input$input_interactive_y_samp_desc)
          }
        }
        
        temp_labels <- paste0(temp_labels, "\n\ncounts (ppm)")
        
        temp_samp_lst <- c()
        if(input$input_interactive_x_type == "enrichment") {
          if(!input$input_interactive_x_ctrl_desc %in% temp_samp_lst) {
            temp_labels <- paste0(temp_labels, "\n",
                                  input$input_interactive_x_ctrl_desc, ": ", temp$x_ctrl_ct_ppm)
            temp_samp_lst <- c(temp_samp_lst, input$input_interactive_x_ctrl_desc)
          }
          if(!input$input_interactive_x_targ_desc %in% temp_samp_lst) {
            temp_labels <- paste0(temp_labels, "\n",
                                  input$input_interactive_x_targ_desc, ": ", temp$x_targ_ct_ppm)
            temp_samp_lst <- c(temp_samp_lst, input$input_interactive_x_targ_desc)
          }
        }
        else {
          if(!input$input_interactive_x_samp_desc %in% temp_samp_lst) {
            temp_labels <- paste0(temp_labels, "\n",
                                  input$input_interactive_x_samp_desc, ": ", temp$x_samp_ct_ppm)
            temp_samp_lst <- c(temp_samp_lst, input$input_interactive_x_samp_desc)
          }
        }
        
        if(input$input_interactive_y_type == "enrichment") {
          if(!input$input_interactive_y_ctrl_desc %in% temp_samp_lst) {
            temp_labels <- paste0(temp_labels, "\n",
                                  input$input_interactive_y_ctrl_desc, ": ", temp$y_ctrl_ct_ppm)
            temp_samp_lst <- c(temp_samp_lst, input$input_interactive_y_ctrl_desc)
          }
          if(!input$input_interactive_y_targ_desc %in% temp_samp_lst) {
            temp_labels <- paste0(temp_labels, "\n",
                                  input$input_interactive_y_targ_desc, ": ", temp$y_targ_ct_ppm)
            temp_samp_lst <- c(temp_samp_lst, input$input_interactive_y_targ_desc)
          }
        }
        else {
          if(!input$input_interactive_y_samp_desc %in% temp_samp_lst) {
            temp_labels <- paste0(temp_labels, "\n",
                                  input$input_interactive_y_samp_desc, ": ", temp$y_samp_ct_ppm)
            temp_samp_lst <- c(temp_samp_lst, input$input_interactive_y_samp_desc)
          }
        }
        
        rv$rv_interactive_tbl <- temp %>%
          mutate(label = temp_labels)
      }
    }
  })

  output$output_interactive_sp <- renderPlotly({
    req(rv$rv_interactive_tbl)
    if(input$input_interactive_vis == "On") {
      if(input$input_interactive_x_type %in% c("count", "rank by counts")) {
        temp_x_lab <- paste(input$input_interactive_x_samp_desc, "-", input$input_interactive_x_type)
      }
      else if(input$input_interactive_x_type == "enrichment") {
        temp_x_lab <- paste(input$input_interactive_x_targ_desc, "vs", 
                            input$input_interactive_x_ctrl_desc, 
                            "-", input$input_interactive_x_type, "(lb)")
      }
      if(input$input_interactive_y_type %in% c("count", "rank by counts")) {
        temp_y_lab <- paste(input$input_interactive_y_samp_desc, "-", input$input_interactive_y_type)
      }
      else if(input$input_interactive_y_type == "enrichment") {
        temp_y_lab <- paste(input$input_interactive_y_targ_desc, "vs", 
                            input$input_interactive_y_ctrl_desc, 
                            "-", input$input_interactive_y_type, "(lb)")
      }
      
      plotly_empty()
      if(!input$input_interactive_x_log & !input$input_interactive_y_log) {
        plot_ly(rv$rv_interactive_tbl, source = "interactive", x = ~x_value, y = ~y_value, color = ~color, colors = c("#132B43", "#56B1F7"),
                type = "scattergl", mode = "markers", text = ~label,
                hoverinfo = "text",
                height = input$visualize_interactive_sp_size) %>%
          layout(xaxis = list(title = temp_x_lab),
                 yaxis = list(title = temp_y_lab),
                 hoverlabel = list(align = "left"))
      }
      else if(input$input_interactive_x_log & !input$input_interactive_y_log) {
        if(input$input_interactive_x_type == "enrichment") {
          plot_ly(rv$rv_interactive_tbl, source = "interactive", x = ~x_value, y = ~y_value, color = ~color, colors = c("#132B43", "#56B1F7"),
                  type = "scattergl", mode = "markers", text = ~label,
                  hoverinfo = "text",
                  height = input$visualize_interactive_sp_size) %>%
            layout(xaxis = list(title = temp_x_lab, type = "log"),
                   yaxis = list(title = temp_y_lab),
                   hoverlabel = list(align = "left"))
        }
        else {
          plot_ly(rv$rv_interactive_tbl, source = "interactive", x = ~x_value, y = ~y_value, color = ~color, colors = c("#132B43", "#56B1F7"),
                  type = "scattergl", mode = "markers", text = ~label,
                  hoverinfo = "text",
                  height = input$visualize_interactive_sp_size) %>%
            layout(xaxis = list(title = temp_x_lab),
                   yaxis = list(title = temp_y_lab),
                   hoverlabel = list(align = "left"))
        }
      }
      else if(!input$input_interactive_x_log & input$input_interactive_y_log) {
        if(input$input_interactive_y_type == "enrichment") {
          plot_ly(rv$rv_interactive_tbl, source = "interactive", x = ~x_value, y = ~y_value, color = ~color, colors = c("#132B43", "#56B1F7"),
                  type = "scattergl", mode = "markers", text = ~label,
                  hoverinfo = "text",
                  height = input$visualize_interactive_sp_size) %>%
            layout(xaxis = list(title = temp_x_lab),
                   yaxis = list(title = temp_y_lab, type = "log"),
                   hoverlabel = list(align = "left"))
        }
        else {
          plot_ly(rv$rv_interactive_tbl, source = "interactive", x = ~x_value, y = ~y_value, color = ~color, colors = c("#132B43", "#56B1F7"),
                  type = "scattergl", mode = "markers", text = ~label,
                  hoverinfo = "text",
                  height = input$visualize_interactive_sp_size) %>%
            layout(xaxis = list(title = temp_x_lab),
                   yaxis = list(title = temp_y_lab),
                   hoverlabel = list(align = "left"))
        }
      }
      else {
        if(input$input_interactive_y_type != "enrichment") {
          plot_ly(rv$rv_interactive_tbl, source = "interactive", x = ~x_value, y = ~y_value, color = ~color, colors = c("#132B43", "#56B1F7"),
                  type = "scattergl", mode = "markers", text = ~label,
                  hoverinfo = "text",
                  height = input$visualize_interactive_sp_size) %>%
            layout(xaxis = list(title = temp_x_lab, type = "log"),
                   yaxis = list(title = temp_y_lab),
                   hoverlabel = list(align = "left"))
        }
        else if(input$input_interactive_x_type != "enrichment") {
          plot_ly(rv$rv_interactive_tbl, source = "interactive", x = ~x_value, y = ~y_value, color = ~color, colors = c("#132B43", "#56B1F7"),
                  type = "scattergl", mode = "markers", text = ~label,
                  hoverinfo = "text",
                  height = input$visualize_interactive_sp_size) %>%
            layout(xaxis = list(title = temp_x_lab),
                   yaxis = list(title = temp_y_lab, type = "log"),
                   hoverlabel = list(align = "left"))
        }
        else {
          plot_ly(rv$rv_interactive_tbl, source = "interactive", x = ~x_value, y = ~y_value, color = ~color, colors = c("#132B43", "#56B1F7"),
                  type = "scattergl", mode = "markers", text = ~label,
                  hoverinfo = "text",
                  height = input$visualize_interactive_sp_size) %>%
            layout(xaxis = list(title = temp_x_lab, type = "log"),
                   yaxis = list(title = temp_y_lab, type = "log"),
                   hoverlabel = list(align = "left"))
        }
      }
    }
  })
  
  output$output_interactive_selection <- DT::renderDataTable({
    req(rv$rv_interactive_tbl)
    if(input$input_interactive_vis == "On") {
      temp <- rv$rv_interactive_tbl %>%
        select(-label) %>%
        select(-lib_id)
      temp <- temp[which(rownames(temp) %in% (event_data("plotly_selected", source = "interactive")$pointNumber+1)),
                   order(colnames(temp))] %>%
        mutate(x_value = x_value %>% round(2),
               y_value = y_value %>% round(2)) %>%
        select(-color)
      DT::datatable(temp, rownames = F, filter = "top", selection = "none", options = list(dom = "tp"))
    }
  })
  
  output$output_interactive_sp_save_plot <- downloadHandler("plot_name.html", 
                                                              content = function(file) {
                                                                req(rv$rv_interactive_tbl)
                                                                if(input$input_interactive_vis == "On") {
                                                                  if(input$input_interactive_x_type %in% c("count", "rank by counts")) {
                                                                    temp_x_lab <- paste(input$input_interactive_x_samp_desc, "-", input$input_interactive_x_type)
                                                                  }
                                                                  else if(input$input_interactive_x_type == "enrichment") {
                                                                    temp_x_lab <- paste(input$input_interactive_x_targ_desc, "vs", 
                                                                                        input$input_interactive_x_ctrl_desc, 
                                                                                        "-", input$input_interactive_x_type, "(lb)")
                                                                  }
                                                                  if(input$input_interactive_y_type %in% c("count", "rank by counts")) {
                                                                    temp_y_lab <- paste(input$input_interactive_y_samp_desc, "-", input$input_interactive_y_type)
                                                                  }
                                                                  else if(input$input_interactive_y_type == "enrichment") {
                                                                    temp_y_lab <- paste(input$input_interactive_y_targ_desc, "vs", 
                                                                                        input$input_interactive_y_ctrl_desc, 
                                                                                        "-", input$input_interactive_y_type, "(lb)")
                                                                  }
                                                                  
                                                                  plotly_empty()
                                                                  if(!input$input_interactive_x_log & !input$input_interactive_y_log) {
                                                                    temp <- plot_ly(rv$rv_interactive_tbl, source = "interactive", x = ~x_value, y = ~y_value, color = ~color, colors = c("#132B43", "#56B1F7"),
                                                                                    type = "scattergl", mode = "markers", text = ~label,
                                                                                    hoverinfo = "text") %>%
                                                                      layout(xaxis = list(title = temp_x_lab),
                                                                             yaxis = list(title = temp_y_lab))
                                                                  }
                                                                  else if(input$input_interactive_x_log & !input$input_interactive_y_log) {
                                                                    if(input$input_interactive_x_type == "enrichment") {
                                                                      temp <- plot_ly(rv$rv_interactive_tbl, source = "interactive", x = ~x_value, y = ~y_value, color = ~color, colors = c("#132B43", "#56B1F7"),
                                                                                      type = "scattergl", mode = "markers", text = ~label,
                                                                                      hoverinfo = "text") %>%
                                                                        layout(xaxis = list(title = temp_x_lab, type = "log"),
                                                                               yaxis = list(title = temp_y_lab))
                                                                    }
                                                                    else {
                                                                      temp <- plot_ly(rv$rv_interactive_tbl, source = "interactive", x = ~x_value, y = ~y_value, color = ~color, colors = c("#132B43", "#56B1F7"),
                                                                                      type = "scattergl", mode = "markers", text = ~label,
                                                                                      hoverinfo = "text") %>%
                                                                        layout(xaxis = list(title = temp_x_lab),
                                                                               yaxis = list(title = temp_y_lab))
                                                                    }
                                                                  }
                                                                  else if(!input$input_interactive_x_log & input$input_interactive_y_log) {
                                                                    if(input$input_interactive_y_type == "enrichment") {
                                                                      temp <- plot_ly(rv$rv_interactive_tbl, source = "interactive", x = ~x_value, y = ~y_value, color = ~color, colors = c("#132B43", "#56B1F7"),
                                                                                      type = "scattergl", mode = "markers", text = ~label,
                                                                                      hoverinfo = "text") %>%
                                                                        layout(xaxis = list(title = temp_x_lab),
                                                                               yaxis = list(title = temp_y_lab, type = "log"))
                                                                    }
                                                                    else {
                                                                      temp <- plot_ly(rv$rv_interactive_tbl, source = "interactive", x = ~x_value, y = ~y_value, color = ~color, colors = c("#132B43", "#56B1F7"),
                                                                                      type = "scattergl", mode = "markers", text = ~label,
                                                                                      hoverinfo = "text") %>%
                                                                        layout(xaxis = list(title = temp_x_lab),
                                                                               yaxis = list(title = temp_y_lab))
                                                                    }
                                                                  }
                                                                  else {
                                                                    if(input$input_interactive_y_type != "enrichment") {
                                                                      temp <- plot_ly(rv$rv_interactive_tbl, source = "interactive", x = ~x_value, y = ~y_value, color = ~color, colors = c("#132B43", "#56B1F7"),
                                                                                      type = "scattergl", mode = "markers", text = ~label,
                                                                                      hoverinfo = "text") %>%
                                                                        layout(xaxis = list(title = temp_x_lab, type = "log"),
                                                                               yaxis = list(title = temp_y_lab))
                                                                    }
                                                                    else if(input$input_interactive_x_type != "enrichment") {
                                                                      temp <- plot_ly(rv$rv_interactive_tbl, source = "interactive", x = ~x_value, y = ~y_value, color = ~color, colors = c("#132B43", "#56B1F7"),
                                                                                      type = "scattergl", mode = "markers", text = ~label,
                                                                                      hoverinfo = "text") %>%
                                                                        layout(xaxis = list(title = temp_x_lab),
                                                                               yaxis = list(title = temp_y_lab, type = "log"))
                                                                    }
                                                                    else {
                                                                      temp <- plot_ly(rv$rv_interactive_tbl, source = "interactive", x = ~x_value, y = ~y_value, color = ~color, colors = c("#132B43", "#56B1F7"),
                                                                                      type = "scattergl", mode = "markers", text = ~label,
                                                                                      hoverinfo = "text") %>%
                                                                        layout(xaxis = list(title = temp_x_lab, type = "log"),
                                                                               yaxis = list(title = temp_y_lab, type = "log"))
                                                                    }
                                                                  }
                                                                }
                                                                htmlwidgets::saveWidget(partial_bundle(temp), file)
                                                              },
                                                              contentType = "html")
  
  output$output_interactive_sp_save_tbl <- downloadHandler("table_name.csv", 
                                                            content = function(file) {
                                                              req(rv$rv_interactive_tbl)
                                                              if(input$input_interactive_x_type %in% c("count", "rank by counts")) {
                                                                temp_x_lab <- paste("x: ", input$input_interactive_x_samp_desc, "-", input$input_interactive_x_type)
                                                              }
                                                              else if(input$input_interactive_x_type == "enrichment") {
                                                                temp_x_lab <- paste("x: ", input$input_interactive_x_targ_desc, "vs", 
                                                                                    input$input_interactive_x_ctrl_desc, 
                                                                                    "-", input$input_interactive_x_type, "(lb)")
                                                              }
                                                              if(input$input_interactive_y_type %in% c("count", "rank by counts")) {
                                                                temp_y_lab <- paste("y: ", input$input_interactive_y_samp_desc, "-", input$input_interactive_y_type)
                                                              }
                                                              else if(input$input_interactive_y_type == "enrichment") {
                                                                temp_y_lab <- paste("y: ", input$input_interactive_y_targ_desc, "vs", 
                                                                                    input$input_interactive_y_ctrl_desc, 
                                                                                    "-", input$input_interactive_y_type, "(lb)")
                                                              }
                                                              
                                                              rv$rv_interactive_tbl %>%
                                                                rename(!!temp_x_lab := x_value,
                                                                       !!temp_y_lab := y_value) %>%
                                                                select(-label) %>%
                                                                write.csv(file, row.names = F)
                                                            },
                                                            contentType = "text/csv")
  
}

shinyApp(ui, server)
  