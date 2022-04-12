# Import packages ####
packages <- c("data.table", "tidyverse", "plotly", "shiny", "shinydashboard", "shinyjs", "DT", "jsonlite", "rcdk", "rclipboard")
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
  left_join(run_meta, by = "run_id") %>%
  left_join(samp_meta, by = "samp_id") %>%
  left_join(scr_meta, by = "scr_id") %>%
  left_join(targ_meta, by = "targ_id") %>%
  left_join(resin_meta, by = "resin_id")
res_lst <- list.files(file.path(config$project_dir, config$result_metadata_location), recursive = F, full.names = F) %>% str_remove_all(".csv")
if(length(res_lst) == 0) { res_lst <- NULL }
agg_lst <- as_tibble(fread(file.path(config$project_dir, config$aggregation_metadata))) %>%
  mutate(temp = paste0(agg_id, ": ", an_type_name)) %>%
  pull(temp)
lib_lst <- lib_meta %>%
  mutate(temp = paste0(lib_id, ": ", lib_name)) %>%
  pull(temp)
input_3d_er_sp_interactive_cycle_selections <- list()
for(lib in unique(bb_meta$lib_id)) {
  temp <- bb_meta %>%
    filter(lib_id == lib)
  for(cycle_id in unique(temp$cyc)) {
    input_3d_er_sp_interactive_cycle_selections[[as.character(lib)]][[cycle_id]] <- temp %>%
      filter(cyc == cycle_id) %>%
      pull(tag_id)
  }
}

# Run Shiny app ####
for(i in 1) {
  ui <- dashboardPage(
    dashboardHeader(title = "DEL Analysis App 2.1"),
    dashboardSidebar(
      sidebarMenu(
        id = "menu1",
        menuItem(
          "Dashboard", tabName = "dashboard", icon = icon("dashboard")
        ),
        menuItem(
          "Metadata", tabName = "metadata", icon = icon("database"),
          menuSubItem("Libraries", tabName = "lib_meta"),
          menuSubItem("Combined libraries", tabName = "comb_lib_meta"),
          menuSubItem("Targets", tabName = "targ_meta"),
          menuSubItem("Resins", tabName = "resin_meta"),
          menuSubItem("Screens", tabName = "scr_meta"),
          menuSubItem("Samples", tabName = "samp_meta"),
          menuSubItem("Sequencing runs", tabName = "run_meta")
        ),
        menuItem(
          "Create config file", tabName = "create_config", icon = icon("file-code")
        ),
        menuItem(
          "Analyses", tabName = "del_analyses", icon = icon("laptop-code"),
          menuSubItem("Create analysis config file", tabName = "create_an_input"),
          menuSubItem("Run analysis (simple entry)", tabName = "run_an_simp"),
          menuSubItem("Run analysis (advanced entry)", tabName = "run_an_adv"),
          menuSubItem("Completed analyses", tabName = "completed_an")
        ),
        menuItem(
          "Data retrieval", tabName = "retrieve_data", icon = icon("search"),
          menuSubItem("Get compound structure(s)", tabName = "get_structure"),
          menuSubItem("View table", tabName = "table_viewer")
        ),
        menuItem(
          "Visualizations", tabName = "visualize_data", icon = icon("chart-bar"),
          menuSubItem("Jitter plot", tabName = "vis_jitterplot"),
          menuSubItem("3D box plot", tabName = "interactive_3d_sp"),
          menuSubItem("2D scatter plot", tabName = "interactive_2d_sp"),
          menuSubItem("Monosynthon heatmap", tabName = "mono_heatmap")
        )
      )
    ),
      
    dashboardBody(
      tabItems(
        tabItem(tabName = "dashboard",
                fluidPage(
                  fluidRow(column(width = 4,
                                  valueBoxOutput("dashboard_num_libs", width = NULL),
                                  valueBoxOutput("dashboard_num_screens", width = NULL),
                                  valueBoxOutput("dashboard_num_targs", width = NULL),
                                  valueBoxOutput("dashboard_num_an_exps", width = NULL)),
                           box(title = "Resin usage", width = 8, status = "primary", solidHeader = T,
                               plotlyOutput("dashboard_resin_usage", height = "auto"))
                  ),
                  fluidRow(box(title = "Timeline of DEL screens", width = 12, status = "primary", solidHeader = T,
                               plotOutput("dashboard_screens_timeline"))
                  )
                )
        ),
        tabItem(tabName = "lib_meta",
                fluidPage(
                  tags$head(tags$style('
                          td[data-type="factor"] input {
                            width: 120px !important;
                          }
                          ')),
                  fluidRow(
                    box(title = "Add new entries", width = 4, status = "primary", solidHeader = T, collapsible = T, collapsed = T,
                        textInput("user_lib_seq", "lib_seq", ""),
                        textInput("user_lib_name", "lib_name", ""),
                        textInput("user_chemist", "chemist", ""),
                        textInput("user_lib_completion_date", "lib_completion_date", ""),
                        textInput("user_lib_loc", "lib_loc", ""),
                        actionButton("append_lib_meta", "Preview"),
                        actionButton("save_lib_meta", "Add entry"),
                        actionButton("reload_lib_meta", "Reload table"))
                  ),
                  fluidRow(
                    box(title = "Library metadata", width = 12, status = "primary", solidHeader = T, 
                        DT::dataTableOutput("lib_meta_table"), style = "overflow-x: scroll;")
                  ),
                  downloadButton("output_lib_meta_download", label = "Save a copy")
                )
        ),
        tabItem(tabName = "comb_lib_meta",
                fluidPage(
                  fluidRow(
                    box(title = "Add new entries", width = 4, status = "primary", solidHeader = T, collapsible = T, collapsed = T,
                        DT::dataTableOutput("input_comb_lib"),
                        actionButton("input_comb_lib_add_row", "Add row"),
                        actionButton("input_comb_lib_remove_row", "Remove row"),
                        br(),
                        br(),
                        actionButton("append_comb_lib_meta", "Preview"),
                        actionButton("save_comb_lib_meta", "Add entry"),
                        actionButton("reload_comb_lib_meta", "Reload table"))
                  ),
                  fluidRow(
                    box(title = "Combined library metadata", width = 12, status = "primary", solidHeader = T, 
                        DT::dataTableOutput("comb_lib_meta_table"), style = "overflow-x: scroll;")
                  ),
                  downloadButton("output_comb_lib_meta_download", label = "Save a copy")
                )
        ),
        tabItem(tabName = "targ_meta",
                fluidPage(
                  fluidRow(
                    box(title = "Add new entries", width = 4, status = "primary", solidHeader = T, collapsible = T, collapsed = T,
                        textInput("user_targ_name", "targ_name", ""),
                        actionButton("append_targ_meta", "Preview"),
                        actionButton("save_targ_meta", "Add entry"),
                        actionButton("reload_targ_meta", "Reload table"))
                  ),
                  fluidRow(
                    box(title = "Target metadata", width = 12, status = "primary", solidHeader = T, 
                        DT::dataTableOutput("targ_meta_table"), style = "overflow-x: scroll;")
                  ),
                  downloadButton("output_targ_meta_download", label = "Save a copy")
                )
        ),
        tabItem(tabName = "resin_meta",
                fluidPage(
                  fluidRow(
                    box(title = "Add new entries", width = 4, status = "primary", solidHeader = T, collapsible = T, collapsed = T,
                        textInput("user_resin_name", "resin_name", ""),
                        actionButton("append_resin_meta", "Preview"),
                        actionButton("save_resin_meta", "Add entry"),
                        actionButton("reload_resin_meta", "Reload table"))
                  ),
                  fluidRow(
                    box(title = "Resin metadata", width = 12, status = "primary", solidHeader = T, 
                        DT::dataTableOutput("resin_meta_table"), style = "overflow-x: scroll;")
                  ),
                  downloadButton("output_resin_meta_download", label = "Save a copy")
                )
        ),
        tabItem(tabName = "scr_meta",
                fluidPage(
                  fluidRow(
                    box(title = "Add new entries", width = 4, status = "primary", solidHeader = T, collapsible = T, collapsed = T,
                        textInput("user_scr_name", "scr_name", ""),
                        textInput("user_scr_completion_date", "scr_completion_date", ""),
                        actionButton("append_scr_meta", "Preview"),
                        actionButton("save_scr_meta", "Add entry"),
                        actionButton("reload_scr_meta", "Reload table"))
                  ),
                  fluidRow(
                    box(title = "DEL screen metadata", width = 12, status = "primary", solidHeader = T, 
                        DT::dataTableOutput("scr_meta_table"), style = "overflow-x: scroll;")
                  ),
                  downloadButton("output_scr_meta_download", label = "Save a copy")
                )
        ),
        tabItem(tabName = "samp_meta",
                fluidPage(
                  fluidRow(
                    box(title = "Add new entries", width = 12, status = "primary", solidHeader = T, collapsible = T, collapsed = T,
                        fileInput("new_samples", "Upload table",
                                           multiple = FALSE,
                                           accept = c("text/csv",
                                                      "text/comma-separated-values,text/plain",
                                                      ".csv")),
                        DT::dataTableOutput("new_samples_table"), style = "overflow-x: scroll;",
                        br(),
                        actionButton("append_samp_meta", "Preview"),
                        actionButton("save_samp_meta", "Add entries"),
                        actionButton("reload_samp_meta", "Reload table"))
                    ),
                  fluidRow(
                    box(title = "Sample metadata", width = 12, status = "primary", solidHeader = T, 
                        DT::dataTableOutput("samp_meta_table"), style = "overflow-x: scroll;")
                  ),
                  downloadButton("output_samp_meta_download", label = "Save a copy")
                )
        ),
        tabItem(tabName = "run_meta",
                fluidPage(
                  fluidRow(
                    box(title = "Add new entries", width = 4, status = "primary", solidHeader = T, collapsible = T, collapsed = T,
                        textInput("user_seq_type", "seq_type", ""),
                        textInput("user_seq_location", "seq_location", ""),
                        textInput("user_seq_completion_date", "seq_completion_date", ""),
                        textInput("user_fastq_folder_name", "fastq_folder_name", ""),
                        actionButton("append_run_meta", "Preview"),
                        actionButton("save_run_meta", "Add entry"),
                        actionButton("reload_run_meta", "Reload table"))
                  ),
                  fluidRow(
                    box(title = "Sequencing run metadata", width = 12, status = "primary", solidHeader = T, 
                        DT::dataTableOutput("run_meta_table"), style = "overflow-x: scroll;")
                  ),
                  downloadButton("output_run_meta_download", label = "Save a copy")
                )
        ),
        tabItem(tabName = "create_config",
                fluidPage(
                  div("Create config file (json) from a csv table."),
                  br(),
                  fluidRow(
                    box(title = "Upload table", width = 6, status = "primary", solidHeader = T, collapsible = T,
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
                        downloadButton("download_json", "Download json"))
                  ),
                  fluidRow(
                    box(title = "Uploaded table", width = 12, status = "primary", solidHeader = T, 
                        DT::dataTableOutput("create_config_table"), style = "overflow-x: scroll;")
                  )
                )
        ),
        tabItem(tabName = "create_an_input",
                fluidPage(
                  list(tags$head(tags$style(HTML("
                                 .multicol { 
                                   height: auto;
                                   -webkit-column-count: 4; /* Chrome, Safari, Opera */ 
                                   -moz-column-count: 4;    /* Firefox */ 
                                   column-count: 4; 
                                   -moz-column-fill: auto;
                                   -column-fill: auto;
                                 } 
                                 ")) 
                  )),
                  fluidRow(
                    box(title = "Column selection", width = 12, status = "primary", solidHeader = T, collapsible = T, collapsed = T,
                        tags$div(align = 'left', 
                                 class = 'multicol',
                                 checkboxGroupInput("input_create_an_show_vars", "Columns to show:",
                                                    names(merged_meta), selected = names(merged_meta), inline = F)),
                        br(),
                        actionButton("input_create_an_select_all_cols", "Select all columns"),
                        actionButton("input_create_an_clear_all_cols", "Clear column selection"),
                        actionButton("reload_merged_meta_create_an", "Reload table"))
                  ),
                  fluidRow(
                    box(title = "Define samples for comparison", width = 12, status = "primary", solidHeader = T, collapsible = T, 
                        DT::dataTableOutput("input_create_an_tbl"), style = "overflow-x: scroll;",
                        textOutput("create_an_input_n_samp_selected"),
                        actionButton("input_create_an_tbl_select_all", "Select all samples"),
                        actionButton("input_create_an_tbl_clear_all", "Clear sample selection"),
                        br(),
                        br(),
                        textInput("input_create_an_tbl_assign_samp_name", "Nickname", "example: (BRD4, heat elution, His beads)"),
                        actionButton("input_create_an_tbl_assign_samp", "Assign selected samples to nickname")
                    )
                  ),
                  fluidRow(
                    box(title = "Define target(s) and control(s)", width = 12, status = "primary", solidHeader = T, collapsible = T, 
                      fluidRow(
                        column(6, selectInput("input_create_an_tbl_targ_choices", label = "Choose target", choices = NULL)),
                        column(6, selectInput("input_create_an_tbl_ctrl_choices", label = "Choose control", choices = NULL))
                      ),
                      actionButton("input_create_an_tbl_add_row", "Add row"),
                      actionButton("input_create_an_tbl_remove_row", "Remove last row")
                    )
                  ),
                  fluidRow(
                    box(title = "Analysis config table", width = 12, status = "primary", solidHeader = T, 
                        DT::dataTableOutput("create_an_input_tbl"),
                        downloadButton("create_an_input_tbl_save", "Save config table"), 
                        style = "overflow-x: scroll;")
                  )
                )
        ),
        tabItem(tabName = "run_an_simp",
                fluidPage(
                  div("Perform analysis, provided a list of controls and targets. Uploaded files must contain columns 'ct_filename' and 'description'."),
                  div("Note: Samples that are replicates should be given the same 'description'."),
                  br(),
                  fluidRow(
                    box(title = "Upload tables", width = 6, status = "primary", solidHeader = T, collapsible = T,
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
                        actionButton("run_an_simp_run", "Run analysis")),
                    valueBoxOutput("run_an_simp_result_name")
                  ),
                  fluidRow(
                    box(title = "Control", width = 12, status = "primary", solidHeader = T, 
                        DT::dataTableOutput("run_an_simp_ctrl_table"), style = "overflow-x: scroll;")
                  ),
                  fluidRow(
                    box(title = "Target(s)", width = 12, status = "primary", solidHeader = T, 
                        DT::dataTableOutput("run_an_simp_targ_table"), style = "overflow-x: scroll;")
                  )
                )
        ),
        tabItem(tabName = "run_an_adv",
                fluidPage(
                  div("Perform analysis, provided a table containing the controls and targets. Uploaded file must contain columns 'targ_ct_lst', 'ctrl_ct_lst', 'targ_description', and 'ctrl_description'. See example input for guidance."),
                  br(),
                  fluidRow(
                    box(title = "Upload table", width = 6, status = "primary", solidHeader = T, collapsible = T,
                        fileInput("run_an_adv_upload", "Choose csv file",
                                  multiple = FALSE,
                                  accept = c("text/csv",
                                             "text/comma-separated-values,text/plain",
                                             ".csv")),
                        numericInput("input_an_adv_ci", "Choose a confidence interval for enrichments (%). Decimals allowed.", 95),
                        br(),
                        actionButton("run_an_adv_assign_res_name", "Assign result ID"),
                        actionButton("run_an_adv_run", "Run analysis")),
                    valueBoxOutput("run_an_adv_result_name")
                  ),
                  fluidRow(
                    box(title = "Uploaded table", width = 12, status = "primary", solidHeader = T, 
                        DT::dataTableOutput("run_an_adv_table"), style = "overflow-x: scroll;")
                  )
                )
        ),
        tabItem(tabName = "completed_an",
                fluidPage(
                  div("Here you can view completed analyses and check how analysis experiments were defined."),
                  br(),
                  selectInput("input_completed_analyses", "Choose a res ID", c("All", res_lst)),
                  fluidRow(
                    box(title = "Completed analyses", width = 12, status = "primary", solidHeader = T, 
                        DT::dataTableOutput("res_meta_viewer"), style = "overflow-x: scroll;")
                  )
                )
        ),
        tabItem(tabName = "get_structure",
                fluidPage(
                  fluidRow(
                    box(title = "Single compound", width = 12, status = "primary", solidHeader = T, collapsible = T,
                        fluidRow(
                          column(
                            width = 6, 
                            numericInput("input_cpd_struct_single_lib_id", "lib_id", NA),
                            numericInput("input_cpd_struct_single_cy1", "cycle1", NA),
                            numericInput("input_cpd_struct_single_cy2", "cycle2", NA),
                            numericInput("input_cpd_struct_single_cy3", "cycle3", NA)),
                          column(
                            width = 6,
                            tags$head(tags$script('$(document).on("shiny:connected", function(e) {
                          Shiny.onInputChange("get_structure_single_cpd_panel_width", window.innerWidth);
                          });
                          $(window).resize(function(e) {
                          Shiny.onInputChange("get_structure_single_cpd_panel_width", window.innerWidth);
                          });
                          ')),
                            plotOutput("output_get_structure_single_cpd_struct", inline = T))
                          ),
                        actionButton("get_cpd_struct", "Get structure"),
                        actionButton("input_get_structure_single_cpd_smi_to_clip", "Copy SMILES to clipboard"))
                  ),
                  fluidRow(
                    box(title = "Multiple compounds", width = 12, status = "primary", solidHeader = T, collapsible = T,
                        fileInput("cpd_struct_upload", "Multiple compounds:",
                                  multiple = FALSE,
                                  accept = c("text/csv",
                                             "text/comma-separated-values,text/plain",
                                             ".csv")),
                        actionButton("get_cpd_structs", "Get structures"),
                        actionButton("clear_cpd_struct", "Clear"),
                        downloadButton("save_cpd_struct", "Save"),
                        DT::dataTableOutput("get_structure_table"), style = "overflow-x: scroll;")
                  )
                )
        ),
        tabItem(tabName = "table_viewer",
                fluidPage(
                  fluidRow(
                    box(title = "Select dataset", width = 6, status = "primary", solidHeader = T, collapsible = T,
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
                        div("Key:"),
                        div("er = enrichment ratio (maximum likelihood)"),
                        div("er_lb = enrichment ratio (lower bound of confidence interval)"),
                        div("er_ub = enrichment ratio (upper bound of confidence interval)"))
                  ),
                  fluidRow(
                    box(title = "Table", width = 12, status = "primary", solidHeader = T, 
                        DT::dataTableOutput("table_viewer"), style = "overflow-x: scroll;",
                        downloadButton("save_tbl_viewer", "Save"))
                  )
                )
        ),
        tabItem(tabName = "vis_jitterplot",
                box(title = "Select samples", status = "primary", solidHeader = T, collapsible = T, width = NULL,
                    selectInput("input_vis_jitterplot_completed_analyses", "Choose a res ID", c("All", res_lst)),
                    DT::dataTableOutput("vis_jitterplot_completed_an_tbl"), style = "overflow-x: scroll;",
                    textOutput("vis_jitterplot_n_samp_selected"),
                    actionButton("input_vis_jitterplot_completed_an_tbl_select_all", "Select all"),
                    actionButton("input_vis_jitterplot_completed_an_tbl_clear_all", "Clear selection"),
                    actionButton("vis_jitterplot_add_sample_to_list", "Add selected analyses to list")),
                box(title = "Add entries to plot", status = "primary", solidHeader = T, collapsible = T, width = NULL,
                    fluidRow(
                      column(6, selectInput("input_vis_jitterplot_sample_choices", label = "Choose sample", choices = NULL)),
                      column(6, selectInput("input_vis_jitterplot_color_choices", label = "Color by:", choices = NULL))
                    ),
                    textInput("input_vis_jitterplot_nickname", "Nickname (axis label):"),
                    actionButton("input_vis_jitterplot_add_entry", "Add to plot")),
                box(title = "Remove entries from plot", status = "primary", solidHeader = T, collapsible = T, collapsed = T, width = NULL,
                    fluidRow(
                      column(8, selectInput("input_vis_jitterplot_remove_entry_choice", NULL, choices = NULL),
                             actionButton("input_vis_jitterplot_remove_all", "Remove all")),
                      column(4, actionButton("input_vis_jitterplot_remove_entry", "Remove entry"))
                    )),
                box(title = "Jitter plot", status = "primary", solidHeader = T, width = NULL,
                    fluidRow(
                      column(width = 6,
                             selectInput("vis_jitterplot_agg_id", "Choose level of aggregation", agg_lst)),
                      column(width = 6,
                             sliderInput("vis_jitterplot_num_points", "Number of points to visualize", min = 100, max = 100000, value = 10000, step = 100))
                    ),
                    sliderInput("vis_jitterplot_selectivity_cutoff", "Selectivity cutoff", min = 0, max = 20, value = 0, step = 0.05),
                    plotlyOutput("output_vis_jitterplot", height = "auto"),
                    br(),
                    downloadButton("output_vis_jitterplot_save_plot", "Save plot"),
                    downloadButton("output_vis_jitterplot_save_tbl", "Save data as table"))
        ),
        tabItem(tabName = "interactive_3d_sp",
                fluidRow(
                  column(width = 6,
                         tags$head(tags$script('
                          var main_panel_innerWidth_3d_er_sp_interactive = 0;
                          $(document).on("shiny:connected", function(e) {
                          main_panel_innerWidth_3d_er_sp_interactive = window.innerWidth;
                          Shiny.onInputChange("main_panel_innerWidth_3d_er_sp_interactive", window.innerWidth);
                          });
                          $(window).resize(function(e) {
                          Shiny.onInputChange("main_panel_innerWidth_3d_er_sp_interactive", window.innerWidth);
                          });   
                          ')),
                         box(title = "Select dataset", status = "primary", solidHeader = T, collapsible = T, width = NULL,
                             selectInput("input_3d_er_interactive_res_id", "Choose res_id", res_lst),
                             selectInput("input_3d_er_interactive_targ_desc", "Choose target", NULL),
                             selectInput("input_3d_er_interactive_ctrl_desc", "Choose control", NULL),
                             selectInput("input_3d_er_interactive_lib_id", "Choose library", NULL),
                             actionButton("visualize_3d_er_sp_interactive", "Visualize plot"),
                             br(),
                             br(),
                             div("Each point is a trisynthon where the color and size both correspond to the enrichment ratio (lower bound). Click on points to render structures.")),
                         box(title = "Scatter plot", status = "primary", solidHeader = T, width = NULL,
                             fluidRow(
                               column(width = 4,
                                      sliderInput("visualize_3d_er_sp_interactive_num_points", "Number of points to visualize", min = 100, max = 100000, value = 1000, step = 100)),
                               column(width = 4,
                                      sliderInput("visualize_3d_er_sp_interactive_scale", "Adjust point size", min = 0, max = 5, value = 1, step = 0.01)),
                               column(width = 4,
                                      radioButtons("visualize_3d_er_sp_interactive_scale_log_radio", "Scale point size by:", choices = c("enrichment", "log(enrichment)")))
                             ),
                             plotlyOutput("output_3d_er_sp_interactive", height = "auto"),
                             br(),
                             downloadButton("output_3d_er_sp_interactive_save_plot", "Save plot"),
                             downloadButton("output_3d_er_sp_interactive_save_tbl", "Save data as table"))
                  ),
                  column(width = 3,
                         box(title = "Filter by cycle", status = "primary", solidHeader = T, collapsible = T, collapsed = T, width = NULL,
                             actionButton("input_3d_er_sp_interactive_filter_reset", "Reset"),
                             br(),
                             selectInput("input_3d_er_sp_interactive_filter_cycle_id", "Cycle", choices = 1:3),
                             uiOutput("output_3d_er_sp_interactive_filter_ui")),
                         box(title = "Filter by MW", status = "primary", solidHeader = T, collapsible = T, collapsed = T, width = NULL,
                             plotOutput("input_3d_er_sp_interactive_filter_mw_plot", height = "auto"),
                             sliderInput("input_3d_er_sp_interactive_filter_mw_slider", NULL, min = 0, max = 1000, value = c(0,1000), step = 10)),
                         box(title = "Filter by fsp3", status = "primary", solidHeader = T, collapsible = T, collapsed = T, width = NULL,
                             plotOutput("input_3d_er_sp_interactive_filter_fsp3_plot", height = "auto"),
                             sliderInput("input_3d_er_sp_interactive_filter_fsp3_slider", NULL, min = 0, max = 1, value = c(0,1), step = 0.01)),
                         box(title = "Filter by nrb", status = "primary", solidHeader = T, collapsible = T, collapsed = T, width = NULL,
                             plotOutput("input_3d_er_sp_interactive_filter_nrb_plot", height = "auto"),
                             sliderInput("input_3d_er_sp_interactive_filter_nrb_slider", NULL, min = 0, max = 20, value = c(0,20), step = 1)),
                         box(title = "Filter by logP", status = "primary", solidHeader = T, collapsible = T, collapsed = T, width = NULL,
                             plotOutput("input_3d_er_sp_interactive_filter_slogp_plot", height = "auto"),
                             sliderInput("input_3d_er_sp_interactive_filter_slogp_slider", NULL, min = -5, max = 10, value = c(-5,10), step = 0.01))
                  ),
                  column(width = 3,
                         tags$head(tags$script('$(document).on("shiny:connected", function(e) {
                          Shiny.onInputChange("side_panel_innerWidth_3d_er_sp_interactive", window.innerWidth);
                          });
                          $(window).resize(function(e) {
                          Shiny.onInputChange("side_panel_innerWidth_3d_er_sp_interactive", window.innerWidth);
                          });
                          ')),
                         box(title = "Structure viewer", status = "primary", solidHeader = T, collapsible = T, width = NULL,
                             plotOutput("output_3d_er_sp_interactive_selected_struct", inline = T),
                             verbatimTextOutput("output_3d_er_sp_interactive_selected_struct_props"),
                             br(),
                             actionButton("input_3d_er_sp_interactive_smi_to_clip", "Copy SMILES to clipboard"),
                             actionButton("input_3d_er_sp_interactive_reorder_bbs", "Reorder by structural similarity"),
                             actionButton("input_3d_er_sp_interactive_reset_order", "Reset order")),
                         box(title = "Enrichment viewer", status = "primary", solidHeader = T, collapsible = T, width = NULL,
                             verbatimTextOutput("output_3d_er_sp_interactive_selection"),
                             plotlyOutput("output_3d_er_sp_interactive_selected_er", height = "auto"))
                  )
                )
        ),
        tabItem(tabName = "interactive_2d_sp",
                fluidRow(
                  column(width = 4,
                         box(title = "Select dataset", status = "primary", solidHeader = T, width = NULL,
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
                             radioButtons("input_interactive_vis", "Plot visualization", choices = c("Off", "On"), inline = T))
                  ),
                  column(width = 8,
                         box(title = "Highlight specific points", status = "primary", solidHeader = T, collapsible = T, width = NULL,
                             selectInput("input_interactive_color_dropdown", "Filter:", c("Off", "Manually", "By table")),
                             conditionalPanel(
                               'input.input_interactive_color_dropdown == "Off"'
                             ),
                             conditionalPanel(
                               'input.input_interactive_color_dropdown == "Manually"',
                               div("Separate BBs with semicolons and use hyphens to denote ranges (e.g. 1-5;7;9-12)."),
                               fluidRow(
                                 column(width = 4,
                                        textInput("input_interactive_color_cy1", "Cycle 1", NA)),
                                 column(width = 4,
                                        textInput("input_interactive_color_cy2", "Cycle 2", NA)),
                                 column(width = 4,
                                        textInput("input_interactive_color_cy3", "Cycle 3", NA))
                               )
                             ),
                             conditionalPanel(
                               'input.input_interactive_color_dropdown == "By table"',
                               fluidRow(
                                 column(width = 6,
                                        fileInput("input_interactive_color_tbl_request", "Color specific compounds using input table",
                                                  multiple = FALSE,
                                                  accept = c("text/csv",
                                                             "text/comma-separated-values,text/plain",
                                                             ".csv")))
                               )
                             )),
                         box(title = "Scatter plot", status = "primary", solidHeader = T, width = NULL,
                             fluidPage(tags$head(tags$script('
                                                  var interactive_2d_sp_box_height = 0;
                                                  $(document).on("shiny:connected", function(e) {
                                                  interactive_2d_sp_box_height = window.innerHeight;
                                                  Shiny.onInputChange("interactive_2d_sp_box_height", interactive_2d_sp_box_height);
                                                  });
                                                  $(window).resize(function(e) {
                                                  interactive_2d_sp_box_height = window.innerHeight
                                                  Shiny.onInputChange("interactive_2d_sp_box_height", interactive_2d_sp_box_height);
                                                  });
                                                  ')),
                                       conditionalPanel('input.input_interactive_vis == "On"',
                                                        sliderInput("input_interactive_max_points", "Max number of points",
                                                                             100, 100000, value = 1000, step = 100),
                                                        plotlyOutput("output_interactive_sp", height = "auto"),
                                                        br(),
                                                        downloadButton("output_interactive_sp_save_plot_html", "Save plot (html)"),
                                                        downloadButton("output_interactive_sp_save_tbl", "Save data as table"),
                                                        style = "height: calc(100vh - 280px) !important;"
                                                        )))
                  )
                )
        ),
        tabItem(tabName = "mono_heatmap",
                box(title = "Add entries to plot", status = "primary", solidHeader = T, collapsible = T, width = NULL,
                    selectInput("input_vis_mono_heatmap_completed_analyses", "Choose a res ID", c("All", res_lst)),
                    DT::dataTableOutput("vis_mono_heatmap_completed_an_tbl"), style = "overflow-x: scroll;",
                    textInput("input_vis_mono_heatmap_nickname", "Nickname (axis label):"),
                    actionButton("input_vis_mono_heatmap_add_entry", "Add to plot")),
                box(title = "Remove entries from plot", status = "primary", solidHeader = T, collapsible = T, collapsed = T, width = NULL,
                    fluidRow(
                      column(8, selectInput("input_vis_mono_heatmap_remove_entry_choice", NULL, choices = NULL),
                             actionButton("input_vis_mono_heatmap_remove_all", "Remove all")),
                      column(4, actionButton("input_vis_mono_heatmap_remove_entry", "Remove entry"))
                    )),
                box(title = "Heatmap", status = "primary", solidHeader = T, width = NULL,
                    fluidRow(
                      column(8, selectInput("vis_mono_heatmap_agg_id", "Choose monosynthon", agg_lst[1:3])),
                      column(4, sliderInput("input_vis_mono_heatmap_col_scale", "Color scale", min = 0.05, max = 10, value = 1, step = 0.05))
                    ),
                    plotlyOutput("output_vis_mono_heatmap", height = "auto"),
                    br(),
                    downloadButton("output_vis_mono_heatmap_save_plot", "Save plot"),
                    downloadButton("output_vis_mono_heatmap_save_tbl", "Save data as table"))
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
                       rv_merged_meta = merged_meta,
                       rv_create_an_tbl = tibble(targ_ct_lst = character(),
                                                 ctrl_ct_lst = character(),
                                                 targ_description = character(),
                                                 ctrl_description = character()),
                       rv_create_an_tbl_search_columns = c(),
                       rv_single_cpd_struct = NULL,
                       rv_cpd_structure = tibble(lib_id = numeric(),
                                                 cycle1 = numeric(),
                                                 cycle2 = numeric(),
                                                 cycle3 = numeric(),
                                                 structure = character()),
                       rv_lib_enum_struct_map = list(),
                       rv_lib_enum_prop_map = list(),
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
                       rv_vis_jitterplot_completed_an_tbl = NULL,
                       rv_vis_jitterplot_selected_samples_tbl = tibble(),
                       rv_vis_jitterplot_vis_tbl = rep(list(tibble()), length(agg_lst)),
                       rv_3d_er_sp_interactive_cycle_selections = input_3d_er_sp_interactive_cycle_selections,
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
                       rv_interactive_tbl = NULL,
                       rv_vis_mono_heatmap_completed_an_tbl = NULL,
                       # rv_vis_jitterplot_selected_samples_tbl = tibble(),
                       rv_vis_mono_heatmap_vis_tbl = rep(list(tibble()), length(agg_lst)))

  
  shinyjs::hide("input_3d_er_sp_interactive_smi_to_clip")
  shinyjs::hide("input_3d_er_sp_interactive_reorder_bbs")
  shinyjs::hide("input_3d_er_sp_interactive_reset_order")
  shinyjs::hide("output_3d_er_sp_interactive_save_plot")
  shinyjs::hide("output_3d_er_sp_interactive_save_tbl")
  
  output$dashboard_num_libs <- renderValueBox({
    valueBox(
      value = rv$rv_merged_meta %>% pull(lib_id) %>% unique() %>% length(), subtitle = "Libraries", color = "teal", icon = icon("flask")
    )
  })
  
  output$dashboard_num_screens <- renderValueBox({
    valueBox(
      value = rv$rv_samp_meta %>% nrow(), subtitle = "Screens performed", color = "teal", icon = icon("filter")
    )
  })
  
  output$dashboard_num_targs <- renderValueBox({
    valueBox(
      value = rv$rv_targ_meta %>% nrow(), subtitle = "Targets screened", color = "teal", icon = icon("splotch")
    )
  })
  
  output$dashboard_num_an_exps <- renderValueBox({
    valueBox(
      value = length(rv$rv_res_lst), subtitle = "Analysis experiments", color = "teal", icon = icon("laptop-code")
    )
  })
  
  output$dashboard_resin_usage <- renderPlotly({
    temp <- rv$rv_samp_meta %>%
      left_join(rv$rv_resin_meta, by = "resin_id") %>%
      filter(resin_name != "None") %>%
      group_by(resin_name) %>%
      summarise(n = n(), .groups = "keep") %>%
      plot_ly(labels = ~resin_name, values = ~n, type = "pie", height = 405) %>%
      layout(xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
             yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
             legend = list(bgcolor = "rgba(205,255,255,0.3)"))
             # plot_bgcolor = "white", paper_bgcolor = "rgba(0,0,0,0.1)")
  })
  
  output$dashboard_screens_timeline <- renderPlot({
    temp <- rv$rv_samp_meta %>%
      left_join(rv$rv_scr_meta, by = "scr_id") %>%
      left_join(rv$rv_comb_lib_meta, by = "comb_lib_id") %>%
      left_join(rv$rv_lib_meta, by = "lib_id")
    temp$yr <- sapply(temp$scr_completion_date, function(x) {(str_split(x, "/") %>% unlist())[3]})
    temp$quart <- sapply(temp$scr_completion_date, date_to_quarter)
    
    start_year <- min(temp$yr)
    end_year <- max(temp$yr)
    x_vals <- expand.grid(paste0("/Q", 1:4), seq(start_year, end_year)) %>% as_tibble() %>%
      rename(quart = 1, yr = 2) %>%
      mutate(yr = substr(yr, nchar(yr)-1, nchar(yr)))
    x_vals <- paste0("'", x_vals$yr, x_vals$quart)

    temp %>%
      group_by(lib_name, quart) %>%
      summarise(n = n(), .groups = "keep") %>%
      mutate(quart = factor(quart, levels = x_vals),
             lib_name = factor(lib_name, levels = rev(rv$rv_lib_meta$lib_name))) %>%
      ggplot(aes(x = quart, y = lib_name, size = n, group = lib_name)) +
      geom_point(col = "#225ea8") +
      scale_x_discrete(breaks = x_vals, labels = x_vals) +
      theme_bw() +
      xlab("") +
      ylab("") +
      theme(axis.text.x = element_text(face = "plain", size = 12, angle = 30, hjust = 1, color = "#444444"),
            axis.text.y = element_text(face = "plain", size = 12, color = "#444444"),
            axis.title = element_blank(),
            title = element_text(face = "bold", size = 12),
            legend.title = element_text(face = "plain", size = 12, color = "#444444"),
            legend.text = element_text(size = 12, color = "#444444"),
            legend.key = element_blank(),
            # panel.border = element_rect(color = "black", size = 1),
            panel.border = element_blank(),
            axis.line = element_blank(),
            axis.ticks = element_blank(),
            axis.ticks.length = unit(1, "mm"),
            # panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            plot.margin=unit(c(1,1,1,1),"mm")) +
      scale_size_continuous(name = "Screens performed", range = c(0.5,10))
  })

  output$lib_meta_table <- DT::renderDataTable({
    DT::datatable(rv$rv_lib_meta %>%
                    mutate(lib_id = as.factor(lib_id)), rownames = F, filter = "top", selection = "none",
                  class = "nowrap display")
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
      showNotification("Entry added (preview)")
    }
  })
  
  observeEvent(input$save_lib_meta, {
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
      temp_current <- file.path(config$project_dir, config$library_metadata) %>% fread() %>% as_tibble()
      temp <- tibble(lib_id = min(setdiff(1:999, temp_current$lib_id)),
                     lib_seq = input$user_lib_seq,
                     lib_name = input$user_lib_name,
                     chemist = input$user_chemist,
                     lib_completion_date = input$user_lib_completion_date,
                     lib_loc = input$user_lib_loc)
      temp_current %>%
        bind_rows(temp) %>%
        write.csv(file.path(config$project_dir, config$library_metadata), row.names = F)
      showNotification("Entry added")
    }
  })
  
  observeEvent(input$reload_lib_meta, {
    rv$rv_lib_meta <- file.path(config$project_dir, config$library_metadata) %>% fread() %>% as_tibble()
    showNotification("Reloaded")
  })

  output$output_lib_meta_download <- downloadHandler("lib_meta.csv", 
                                                     content = function(file) {
                                                       write.csv(file.path(config$project_dir, config$library_metadata) %>% 
                                                                   fread() %>% 
                                                                   as_tibble(), 
                                                                 file, row.names = F)
                                                     },
                                                     contentType = "text/csv")
  
  output$comb_lib_meta_table <- DT::renderDataTable({
    DT::datatable(rv$rv_comb_lib_meta %>%
                    mutate(across(matches("lib_id"), as.factor)), 
                  rownames = F, filter = "top", selection = "none")
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
      showNotification("Entry added (preview)")
    }
  })
  
  observeEvent(input$save_comb_lib_meta, {
    if(0 %in% rv$rv_input_comb_lib$lib_id) {
      showNotification("Specify all library IDs", type = "warning")
    }
    if(abs(sum(rv$rv_input_comb_lib$frac) - 1) > 0.01) {
      showNotification("Fractional compositions don't add up to 1", type = "warning")
    }
    if(!0 %in% rv$rv_input_comb_lib$lib_id & abs(sum(rv$rv_input_comb_lib$frac) - 1) < 0.01) {
      temp_current <- file.path(config$project_dir, config$combined_library_metadata) %>% fread() %>% as_tibble()
      temp <- rv$rv_input_comb_lib %>%
        mutate(comb_lib_id = min(setdiff(1:99999, temp_current$comb_lib_id)))
      temp_current %>%
        bind_rows(temp) %>%
        write.csv(file.path(config$project_dir, config$combined_library_metadata), row.names = F)
      showNotification("Entry added")
    }
  })
  
  observeEvent(input$reload_comb_lib_meta, {
    rv$rv_comb_lib_meta <- file.path(config$project_dir, config$combined_library_metadata) %>% fread() %>% as_tibble()
    showNotification("Reloaded")
  })
  
  output$output_comb_lib_meta_download <- downloadHandler("comb_lib_meta.csv", 
                                                     content = function(file) {
                                                       write.csv(file.path(config$project_dir, config$combined_library_metadata) %>% 
                                                                   fread() %>% 
                                                                   as_tibble(), 
                                                                 file, row.names = F)
                                                     },
                                                     contentType = "text/csv")


  
  output$targ_meta_table <- DT::renderDataTable({
    DT::datatable(rv$rv_targ_meta %>%
                    mutate(targ_id = as.factor(targ_id)),
                  rownames = F, filter = "top", selection = "none",
                  class = "nowrap display")
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
      showNotification("Entry added (preview)")
    }
  })
  
  observeEvent(input$save_targ_meta, {
    if(input$user_targ_name == "") {
      showNotification("Specify the name of the target", type = "warning")
    }
    if(input$user_targ_name != "") {
      temp_current <- file.path(config$project_dir, config$target_metadata) %>% fread() %>% as_tibble()
      temp <- tibble(targ_id = min(setdiff(1:99999, temp_current$targ_id)),
                     targ_name = input$user_targ_name)
      temp_current %>%
        bind_rows(temp) %>%
        write.csv(file.path(config$project_dir, config$target_metadata), row.names = F)
      showNotification("Entry added")
    }
  })
  
  observeEvent(input$reload_targ_meta, {
    rv$rv_targ_meta <- file.path(config$project_dir, config$target_metadata) %>% fread() %>% as_tibble()
    showNotification("Reloaded")
  })
  
  output$output_targ_meta_download <- downloadHandler("targ_meta.csv", 
                                                      content = function(file) {
                                                        write.csv(file.path(config$project_dir, config$target_metadata) %>% 
                                                                    fread() %>% 
                                                                    as_tibble(), 
                                                                  file, row.names = F)
                                                      },
                                                      contentType = "text/csv")
  
  output$resin_meta_table <- DT::renderDataTable({
    DT::datatable(rv$rv_resin_meta %>%
                    mutate(resin_id = as.factor(resin_id)), 
                  rownames = F, filter = "top", selection = "none",
                  class = "nowrap display")
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
      showNotification("Entry added (preview)")
    }
  })
  
  observeEvent(input$save_resin_meta, {
    if(input$user_resin_name == "") {
      showNotification("Specify the name of the resin", type = "warning")
    }
    if(input$user_resin_name != "") {
      temp_current <- file.path(config$project_dir, config$resin_metadata) %>% fread() %>% as_tibble()
      temp <- tibble(resin_id = min(setdiff(1:99999, temp_current$resin_id)),
                     resin_name = input$user_resin_name)
      temp_current %>%
        bind_rows(temp) %>%
        write.csv(file.path(config$project_dir, config$resin_metadata), row.names = F)
      showNotification("Entry added")
    }
  })
  
  observeEvent(input$reload_resin_meta, {
    rv$rv_resin_meta <- file.path(config$project_dir, config$resin_metadata) %>% fread() %>% as_tibble()
    showNotification("Reloaded")
  })
  
  output$output_resin_meta_download <- downloadHandler("resin_meta.csv", 
                                                       content = function(file) {
                                                         write.csv(file.path(config$project_dir, config$resin_metadata) %>% 
                                                                     fread() %>% 
                                                                     as_tibble(), 
                                                                   file, row.names = F)
                                                       },
                                                       contentType = "text/csv")
  
  output$scr_meta_table <- DT::renderDataTable({
    DT::datatable(rv$rv_scr_meta %>%
                    mutate(scr_id = as.factor(scr_id)), 
                  rownames = F, filter = "top", selection = "none",
                  class = "nowrap display")
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
      showNotification("Entry added (preview)")
    }
  })
  
  observeEvent(input$save_scr_meta, {
    if(input$user_scr_name == "") {
      showNotification("Specify the name of the DEL screen", type = "warning")
    }
    if(input$user_scr_completion_date == "") {
      showNotification("Specify the completion date for the DEL screen", type = "warning")
    }
    if(input$user_scr_name != "" & input$user_scr_completion_date != "") {
      temp_current <- file.path(config$project_dir, config$screen_metadata) %>% fread() %>% as_tibble()
      temp <- tibble(scr_id = min(setdiff(1:99999, temp_current$scr_id)),
                     scr_name = input$user_scr_name,
                     scr_completion_date = input$user_scr_completion_date)
      temp_current %>%
        bind_rows(temp) %>%
        write.csv(file.path(config$project_dir, config$screen_metadata), row.names = F)
      showNotification("Entry added")
    }
  })
  
  observeEvent(input$reload_scr_meta, {
    rv$rv_scr_meta <- file.path(config$project_dir, config$screen_metadata) %>% fread() %>% as_tibble()
    showNotification("Reloaded")
  })
  
  output$output_scr_meta_download <- downloadHandler("scr_meta.csv", 
                                                     content = function(file) {
                                                       write.csv(file.path(config$project_dir, config$screen_metadata) %>% 
                                                                   fread() %>% 
                                                                   as_tibble(), 
                                                                 file, row.names = F)
                                                     },
                                                     contentType = "text/csv")
  
  output$samp_meta_table <- DT::renderDataTable({
    DT::datatable(rv$rv_samp_meta %>%
                    mutate(across(matches("_id"), as.factor)),
                  rownames = F, filter = "top", selection = "none")
  })
  
  observeEvent(input$save_samp_meta, {
    temp_current <- file.path(config$project_dir, config$sample_metadata) %>% fread() %>% as_tibble()
    unused_samp_ids <- setdiff(1:999999, temp_current$samp_id)
    temp <- as_tibble(fread(input$new_samples$datapath)) %>%
      mutate(samp_id = unused_samp_ids[1:nrow(.)])
    temp_current %>%
      bind_rows(temp) %>%
      write.csv(file.path(config$project_dir, config$sample_metadata), row.names = F)
    showNotification("Entries added")
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
    showNotification("Entries added (preview)")
  })
  
  output$run_meta_table <- DT::renderDataTable({
    DT::datatable(rv$rv_run_meta %>%
                    mutate(run_id = as.factor(run_id)), 
                  rownames = F, filter = "top", selection = "none",
                  class = "nowrap display")
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
      showNotification("Entry added (preview)")
    }
  })
  
  observeEvent(input$save_run_meta, {
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
      temp_current <- file.path(config$project_dir, config$sequencing_run_metadata) %>% fread() %>% as_tibble()
      temp <- tibble(run_id = min(setdiff(1:10000, temp_current$run_id)),
                     seq_type = input$user_seq_type,
                     seq_location = input$user_seq_location,
                     seq_completion_date = input$user_seq_completion_date,
                     fastq_folder_name = input$user_fastq_folder_name)
      temp_current %>%
        bind_rows(temp) %>%
        write.csv(file.path(config$project_dir, config$sequencing_run_metadata), row.names = F)
      showNotification("Entry added")
    }
  })
  
  observeEvent(input$reload_run_meta, {
    rv$rv_run_meta <- file.path(config$project_dir, config$sequencing_run_metadata) %>% fread() %>% as_tibble()
    showNotification("Reloaded")
  })
  
  output$output_run_meta_download <- downloadHandler("run_meta.csv", 
                                                     content = function(file) {
                                                       write.csv(file.path(config$project_dir, config$sequencing_run_metadata) %>% 
                                                                   fread() %>% 
                                                                   as_tibble(), 
                                                                 file, row.names = F)
                                                     },
                                                     contentType = "text/csv")
  
  output$new_samples_table <- DT::renderDataTable({
    if(is.null(input$new_samples)) {
      return(NULL)
    }
    DT::datatable(fread(input$new_samples$datapath),
                  rownames = F, selection = "none",
                  class = "nowrap display")
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

  output$get_structure_table <- DT::renderDataTable({
    DT::datatable(rv$rv_cpd_structure, rownames = F, selection = "none",
                  options = list(dom = "tp"))
  })
  
  observeEvent(input$cpd_struct_upload, {
    temp <- as_tibble(fread(input$cpd_struct_upload$datapath))
    rv$rv_cpd_structure <- rv$rv_cpd_structure %>%
      bind_rows(temp)
  })
  
  observeEvent(input$get_cpd_struct, {
    if(!as.character(input$input_cpd_struct_single_lib_id) %in% names(rv$rv_lib_enum_struct_map)) {
      if(file.exists(file.path(config$project_dir, config$enumerated_structure_location, 
                               sprintf("lib%03d.csv", input$input_cpd_struct_single_lib_id)))) {
        load_cpd_notif <- showNotification("Loading metadata...", duration = NULL)
        temp_meta <- file.path(config$project_dir, config$library_metadata_column_location, 
                               sprintf("lib%03d_agg07_meta_cols.csv", input$input_cpd_struct_single_lib_id)) %>% 
          fread() %>% 
          as_tibble()
        temp_struct <- file.path(config$project_dir, config$enumerated_structure_location, 
                                 sprintf("lib%03d.csv", input$input_cpd_struct_single_lib_id)) %>% 
          fread() %>% 
          as_tibble()
        
        temp_meta_filled <- temp_meta %>%
          expand(cycle1, cycle2, cycle3)
        temp <- temp_meta %>%
          select(-lib_id) %>%
          mutate(structure = temp_struct$structure) %>%
          right_join(temp_meta_filled) %>%
          arrange(cycle1, cycle2, cycle3)
        rv$rv_lib_enum_struct_map[[as.character(input$input_cpd_struct_single_lib_id)]] <- array(data = temp$structure, 
                                                                                                 dim = c(max(temp$cycle3), 
                                                                                                         max(temp$cycle2), 
                                                                                                         max(temp$cycle1))) %>% aperm()
        removeNotification(load_cpd_notif)
      }
      else {
        showNotification("Library doesn't exist", type = "warning")
      }
    }
    tryCatch(
      rv$rv_single_cpd_struct <- rv$rv_lib_enum_struct_map[[as.character(input$input_cpd_struct_single_lib_id)]][input$input_cpd_struct_single_cy1,
                                                                                                          input$input_cpd_struct_single_cy2,
                                                                                                          input$input_cpd_struct_single_cy3],
      error = function(e) {rv$rv_single_cpd_struct <- NA}
    )
  })
  
  output$output_get_structure_single_cpd_struct <- renderPlot({
    req(rv$rv_single_cpd_struct)
    dep <- get.depictor(width=400, height=400)
    img <- parse.smiles(rv$rv_single_cpd_struct)[[1]] %>%
      view.image.2d(dep)
    par(mar = c(0, 0, 0, 0))
    plot("", xlim = c(0, 1), ylim = c(0,1), xlab = "", ylab = "", xaxt = "n", yaxt = "n")
    rasterImage(img, 0, 0, 1, 1)
  },
  height = reactive(ifelse(!is.null(input$get_structure_single_cpd_panel_width),input$get_structure_single_cpd_panel_width*0.3,0)),
  width = reactive(ifelse(!is.null(input$get_structure_single_cpd_panel_width),input$get_structure_single_cpd_panel_width*0.3,0))
  )
  
  observeEvent(input$input_get_structure_single_cpd_smi_to_clip, {
    req(rv$rv_single_cpd_struct)
    clipr::write_clip(rv$rv_single_cpd_struct)
    showNotification("SMILES copied to clipboard")
  })
  
  observeEvent(input$get_cpd_structs, {
    temp_libs <- unique(rv$rv_cpd_structure$lib_id)
    for(lib in temp_libs) {
      if(!lib %in% rv$rv_enum_struct_loaded) {
        load_cpd_notif <- showNotification("Loading metadata...", duration = NULL)
        temp_meta <- file.path(config$project_dir, config$library_metadata_column_location, 
                               sprintf("lib%03d_agg07_meta_cols.csv", lib)) %>% 
          fread() %>% 
          as_tibble()
        temp_struct <- file.path(config$project_dir, config$enumerated_structure_location, 
                                 sprintf("lib%03d.csv", lib)) %>% 
          fread() %>% 
          as_tibble()
        
        temp_meta_filled <- temp_meta %>%
          expand(cycle1, cycle2, cycle3)
        temp <- temp_meta %>%
          select(-lib_id) %>%
          mutate(structure = temp_struct$structure) %>%
          right_join(temp_meta_filled) %>%
          arrange(cycle1, cycle2, cycle3)
        rv$rv_lib_enum_struct_map[[as.character(lib)]] <- array(data = temp$structure, 
                                                         dim = c(max(temp$cycle3), 
                                                                 max(temp$cycle2), 
                                                                 max(temp$cycle1))) %>% aperm()
        removeNotification(load_cpd_notif)
      }
    }
    for(i in 1:nrow(rv$rv_cpd_structure)) {
      tryCatch(
        rv$rv_cpd_structure$structure[i] <- rv$rv_lib_enum_struct_map[[as.character(rv$rv_cpd_structure$lib_id[i])]][rv$rv_cpd_structure$cycle1[i], 
                                                                                                              rv$rv_cpd_structure$cycle2[i], 
                                                                                                              rv$rv_cpd_structure$cycle3[i]],
        error = function(e) {rv$rv_cpd_structure$structure[i] <- NA}
      )
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
  
  
  
  observeEvent(input$input_create_an_select_all_cols, {
    updateCheckboxGroupInput(session, "input_create_an_show_vars", choices = names(rv$rv_merged_meta), selected = names(rv$rv_merged_meta))
  })
  
  observeEvent(input$input_create_an_clear_all_cols, {
    updateCheckboxGroupInput(session, "input_create_an_show_vars", choices = names(rv$rv_merged_meta), selected = NULL)
  })
  
  observeEvent(input$reload_merged_meta_create_an, {
    rv$rv_merged_meta <- rv$rv_ct_file_meta %>%
      left_join(rv$rv_run_meta, by = "run_id") %>%
      left_join(rv$rv_samp_meta, by = "samp_id") %>%
      left_join(rv$rv_scr_meta, by = "scr_id") %>%
      left_join(rv$rv_targ_meta, by = "targ_id") %>%
      left_join(rv$rv_resin_meta, by = "resin_id")
    
    showNotification("Reloaded")
  })
  
  output$input_create_an_tbl <- DT::renderDataTable({
    DT::datatable(rv$rv_merged_meta[,input$input_create_an_show_vars, drop = FALSE] %>%
                    mutate(across(everything(), as.character),
                           across(matches("_id"), as.factor)),
                  rownames = F, filter = "top",
                  options = list(dom = 'tp',
                                 pageLength = 4),
                  class = "nowrap display")
  })
  
  merged_create_an_tbl_proxy <- dataTableProxy("input_create_an_tbl")
  
  observe({
    rv$rv_create_an_tbl_search_columns <- structure(names = colnames(rv$rv_merged_meta),
                                                    rep("", ncol(rv$rv_merged_meta)))
  })
  
  observeEvent(input$input_create_an_tbl_search_columns, {
    temp <- structure(names = input$input_create_an_show_vars,
                      input$input_create_an_tbl_search_columns)
    for(i in names(temp)) {
      rv$rv_create_an_tbl_search_columns[i] <- temp[i]
    }
  })
  
  observeEvent(input$input_create_an_show_vars, {
    temp <- rv$rv_create_an_tbl_search_columns[names(rv$rv_create_an_tbl_search_columns) %in% input$input_create_an_show_vars]
    merged_create_an_tbl_proxy %>% updateSearch(keywords = list(columns = as.character(temp)))
  })
  
  output$create_an_input_n_samp_selected <- renderText({
    paste0(length(input$input_create_an_tbl_rows_selected), " samples selected.")
  })
  
  observeEvent(input$input_create_an_tbl_select_all, {
    merged_create_an_tbl_proxy %>% selectRows(input$input_create_an_tbl_rows_all)
    showNotification("All samples selected")
  })
  
  observeEvent(input$input_create_an_tbl_clear_all, {
    merged_create_an_tbl_proxy %>% selectRows(NULL)
    showNotification("Selection cleared")
  })
  
  observeEvent(input$input_create_an_tbl_assign_samp, {
    temp <- rv$rv_merged_meta[input$input_create_an_tbl_rows_selected,] %>%
      pull(ct_filename) %>%
      paste(collapse = ", ")
    rv$rv_create_an_tbl_samp_lst[[input$input_create_an_tbl_assign_samp_name]] <- temp
  })
  
  observe({
    updateSelectInput(session, "input_create_an_tbl_targ_choices", choices = names(rv$rv_create_an_tbl_samp_lst))
    updateSelectInput(session, "input_create_an_tbl_ctrl_choices", choices = names(rv$rv_create_an_tbl_samp_lst))
  })
  
  observeEvent(input$input_create_an_tbl_add_row, {
    req(input$input_create_an_tbl_targ_choices)
    req(input$input_create_an_tbl_ctrl_choices)
    new_row <- tibble(targ_ct_lst = rv$rv_create_an_tbl_samp_lst[[input$input_create_an_tbl_targ_choices]],
                      ctrl_ct_lst = rv$rv_create_an_tbl_samp_lst[[input$input_create_an_tbl_ctrl_choices]],
                      targ_description = input$input_create_an_tbl_targ_choices,
                      ctrl_description = input$input_create_an_tbl_ctrl_choices)
    rv$rv_create_an_tbl <- rv$rv_create_an_tbl %>%
      bind_rows(new_row)
  })
  
  output$create_an_input_tbl <- DT::renderDataTable({
    DT::datatable(rv$rv_create_an_tbl,
                  rownames = F,
                  options = list(dom = 'tp',
                                 pageLength = 5),
                  class = "nowrap display",
                  selection = "none")
  })
  
  observeEvent(input$input_create_an_tbl_remove_row, {
    if(nrow(rv$rv_create_an_tbl) == 0) {}
    else if(nrow(rv$rv_create_an_tbl) == 1) {
      rv$rv_create_an_tbl <- tibble()
    }
    else {
      rv$rv_create_an_tbl <- rv$rv_create_an_tbl[1:(nrow(rv$rv_create_an_tbl)-1),]
    }
  })
  
  output$create_an_input_tbl_save <- downloadHandler("run_an_adv_input.csv", 
                                                     content = function(file) {
                                                       write.csv(rv$rv_create_an_tbl, file, row.names = F)
                                                     },
                                                     contentType = "text/csv")
  
  
  
  
  
  
  
  
  
  
  output$run_an_simp_result_name <- renderValueBox({
    valueBox(
      value = rv$rv_run_an_simp_res_name, subtitle = "Assigned result ID", color = "teal",
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
    temp <- list.files(file.path(config$project_dir, config$result_metadata_location), recursive = F, full.names = F) %>% 
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
    if(rv$rv_run_an_simp_res_name %in% list.files(file.path(config$project_dir, config$result_metadata_location), recursive = F, full.names = F) %>% 
       str_remove_all(".csv")) {
      showNotification("Result ID already assigned. Assign a new result ID", type = "warning")
    }
    else if(!is.null(rv$rv_run_an_simp_targ_lst) &
            !is.null(rv$rv_run_an_simp_ctrl_lst) &
            !is.null(rv$rv_run_an_simp_res_name)) {
      run_an_simp_start_notif <- showNotification("Running analysis...", duration = NULL)
      
      tibble() %>%
        write.csv(file.path(config$project_dir, config$result_metadata_location, paste0(rv$rv_run_an_simp_res_name, ".csv")), row.names = F)
      
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
      
      dir.create(file.path(config$project_dir, config$result_location, res_name), showWarnings = F)
      
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
  
  output$run_an_adv_result_name <- renderValueBox({
    valueBox(
      value = rv$rv_run_an_adv_res_name, subtitle = "Assigned result ID", color = "teal",
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
    if(rv$rv_run_an_adv_res_name %in% list.files(file.path(config$project_dir, config$result_metadata_location), recursive = F, full.names = F) %>% 
       str_remove_all(".csv")) {
      showNotification("Result ID already assigned. Assign a new result ID", type = "warning")
    }
    else if("ctrl_ct_lst" %in% names(rv$rv_run_an_adv_tbl) &
            "targ_ct_lst" %in% names(rv$rv_run_an_adv_tbl) &
            !is.null(rv$rv_run_an_adv_res_name)) {
      run_an_adv_start_notif <- showNotification("Running analysis...", duration = NULL)
      
      tibble() %>%
        write.csv(file.path(config$project_dir, config$result_metadata_location, paste0(rv$rv_run_an_adv_res_name, ".csv")), row.names = F)
      
      ptm <- proc.time()
      report_name <- format(Sys.time(), "%Y%m%d_%H%M%S")
      report_lines <- c()
      
      # Import metadata
      source(file.path(config$project_dir, config$metadata_script))
      
      res_name <- rv$rv_run_an_adv_res_name
      
      dir.create(file.path(config$project_dir, config$result_location, res_name), showWarnings = F)
      
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
    if(input$input_completed_analyses == "All") {
      completed_res <- list.files(file.path(config$project_dir, config$result_metadata_location))
      temp <- tibble()
      for(i in completed_res) {
                temp <- temp %>%
          bind_rows(as_tibble(fread(file.path(config$project_dir, config$result_metadata_location, i))) %>% 
                      mutate(res_name = str_remove(i, ".csv")))
      }
      rv$rv_res_meta_viewer <- temp
    }
    else {
      rv$rv_res_meta_viewer <- file.path(config$project_dir, config$result_metadata_location, paste0(input$input_completed_analyses, ".csv")) %>%
        fread() %>%
        as_tibble()
    }
  })
  
  observe({
    updateSelectInput(session, "input_completed_analyses", choices = c("All", rv$rv_res_lst))
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
        left_join(rv$rv_lib_meta, by = "lib_id") %>%
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
                    mutate(across(matches("cycle"), as.factor)), 
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
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  observeEvent(input$input_vis_jitterplot_completed_analyses, {
    req(input$input_vis_jitterplot_completed_analyses)
    if(input$input_vis_jitterplot_completed_analyses == "All") {
      completed_res <- list.files(file.path(config$project_dir, config$result_metadata_location))
      temp <- tibble()
      for(i in completed_res) {
        temp <- temp %>%
          bind_rows(as_tibble(fread(file.path(config$project_dir, config$result_metadata_location, i))) %>% 
                      mutate(res_name = str_remove(i, ".csv")))
      }
      rv$rv_vis_jitterplot_completed_an_tbl <- temp
    }
    else {
      rv$rv_vis_jitterplot_completed_an_tbl <- file.path(config$project_dir, config$result_metadata_location, paste0(input$input_vis_jitterplot_completed_analyses, ".csv")) %>%
        fread() %>%
        as_tibble() %>%
        mutate(res_name = input$input_vis_jitterplot_completed_analyses)
    }
  })
  
  observe({
    updateSelectInput(session, "input_vis_jitterplot_completed_analyses", choices = c("All", rv$rv_res_lst))
  })
  
  output$vis_jitterplot_completed_an_tbl <- DT::renderDataTable({
    req(rv$rv_vis_jitterplot_completed_an_tbl)
    temp_lib_id <- sapply(strsplit(rv$rv_vis_jitterplot_completed_an_tbl$targ_ct_lst, "_lib|.csv"), `[`, 2) %>%
      as.numeric()
    DT::datatable(rv$rv_vis_jitterplot_completed_an_tbl %>%
                    mutate(lib_id = sprintf("%03d", temp_lib_id)) %>%
                    relocate(lib_id) %>%
                    select(-agg_id, -an_type_id, -an_id, -targ_ct_lst, -ctrl_ct_lst),
                  rownames = F, filter = "top",
                  options = list(dom = 'tp',
                                 pageLength = 4),
                  class = "nowrap display")
  })
  
  vis_jitterplot_an_tbl_proxy <- dataTableProxy("vis_jitterplot_completed_an_tbl")
  
  observeEvent(input$input_vis_jitterplot_completed_an_tbl_select_all, {
    vis_jitterplot_an_tbl_proxy %>% selectRows(input$vis_jitterplot_completed_an_tbl_rows_all)
    showNotification("All analyses selected")
  })
  
  observeEvent(input$input_vis_jitterplot_completed_an_tbl_clear_all, {
    vis_jitterplot_an_tbl_proxy %>% selectRows(NULL)
    showNotification("Selection cleared")
  })
  
  output$vis_jitterplot_n_samp_selected <- renderText({
    paste0(length(input$vis_jitterplot_completed_an_tbl_rows_selected), " analyses selected.")
  })
  
  observeEvent(input$vis_jitterplot_add_sample_to_list, {
    temp <- rv$rv_vis_jitterplot_completed_an_tbl[input$vis_jitterplot_completed_an_tbl_rows_selected,] %>%
      mutate(label = paste(targ_description, "vs", ctrl_description, "-", res_name))
    rv$rv_vis_jitterplot_selected_samples_tbl <- rv$rv_vis_jitterplot_selected_samples_tbl %>%
      bind_rows(temp) %>%
      unique()
    showNotification("Selection(s) added")
  })
  
  observe({
    req("label" %in% names(rv$rv_vis_jitterplot_selected_samples_tbl))
    updateSelectInput(session, "input_vis_jitterplot_sample_choices", 
                      choices = rv$rv_vis_jitterplot_selected_samples_tbl$label)
    updateSelectInput(session, "input_vis_jitterplot_color_choices", 
                      choices = rv$rv_vis_jitterplot_selected_samples_tbl$label)
  })
  
  observeEvent(input$input_vis_jitterplot_add_entry, {
    req(input$input_vis_jitterplot_sample_choices)
    req(input$input_vis_jitterplot_color_choices)
    if(nrow(rv$rv_vis_jitterplot_vis_tbl[[1]] > 0)){
      if(input$input_vis_jitterplot_nickname %in% rv$rv_vis_jitterplot_vis_tbl[[1]]$x_lab) {
        showNotification("Nickname already used. Choose another.", type = "warning")
        req(!input$input_vis_jitterplot_nickname %in% rv$rv_vis_jitterplot_vis_tbl[[1]]$x_lab)
      }
    }
    temp_tbl_y <- rv$rv_vis_jitterplot_selected_samples_tbl %>%
      filter(label == input$input_vis_jitterplot_sample_choices)
    temp_tbl_col <- rv$rv_vis_jitterplot_selected_samples_tbl %>%
      filter(label == input$input_vis_jitterplot_color_choices)
    temp_lib_id <- sapply(strsplit(temp_tbl_y$targ_ct_lst, "_lib|.csv"), `[`, 2) %>%
      as.numeric()
    temp_label <- temp_tbl_y %>%
      pull(label)
    temp <- temp_tbl_y %>%
      pull(targ_ct_lst) %>%
      str_split(", ") %>%
      unlist()
    targ_ct <- as_tibble(fread(file.path(config$project_dir, config$data_count_location, temp[1])))
    if(length(temp) > 1) {
      for(i in 2:length(temp)) {
        targ_ct <- as_tibble(targ_ct + as_tibble(fread(file.path(config$project_dir, config$data_count_location, temp[i]))))
      }
    }
    temp <- temp_tbl_y %>%
      pull(ctrl_ct_lst) %>%
      str_split(", ") %>%
      unlist()
    ctrl_ct <- as_tibble(fread(file.path(config$project_dir, config$data_count_location, temp[1])))
    if(length(temp) > 1) {
      for(i in 2:length(temp)) {
        ctrl_ct <- as_tibble(ctrl_ct + as_tibble(fread(file.path(config$project_dir, config$data_count_location, temp[i]))))
      }
    }
    
    temp_agg_meta <- as_tibble(fread(file.path(config$project_dir, config$library_metadata_column_location,
                                               sprintf("lib%03d_agg%02d_meta_cols.csv", temp_lib_id, 7))))
    targ_ct_agg <- as_tibble(fread(file.path(config$project_dir, config$library_metadata_column_location,
                                             sprintf("lib%03d_agg%02d_meta_cols.csv", temp_lib_id, 7))))
    y_vals <- as_tibble(fread(file.path(config$project_dir, config$result_analysis_location,
                                        sprintf("an%06d_lib%03d_agg%02d_type02.csv", temp_tbl_y %>% pull(an_id), temp_lib_id, 7))))
    col_vals <- as_tibble(fread(file.path(config$project_dir, config$result_analysis_location,
                                          sprintf("an%06d_lib%03d_agg%02d_type02.csv", temp_tbl_col %>% pull(an_id), temp_lib_id, 7))))
    
    agg_ids <- as_tibble(fread(file.path(config$project_dir, config$aggregation_metadata))) %>% pull(agg_id)
    for(temp_agg_id in agg_ids) {
      if(temp_agg_id == tail(agg_ids,1)) {
        temp_out <- temp_agg_meta %>%
          bind_cols(tibble(targ_ct = targ_ct %>% pull(value),
                           ctrl_ct = ctrl_ct %>% pull(value),
                           y_val = y_vals %>% pull(value),
                           col_val = col_vals %>% pull(value),
                           x_lab = input$input_vis_jitterplot_nickname))
        if(nrow(rv$rv_vis_jitterplot_vis_tbl[[temp_agg_id]]) == 0) {
          temp_out <- temp_out %>%
            mutate(x_val = 1)
        } else {
          temp_out <- temp_out %>%
            mutate(x_val = max(rv$rv_vis_jitterplot_vis_tbl[[temp_agg_id]]$x_val) + 1)
        }
      }
      else {
        cycs <- paste0("cycle", unlist(str_split(agg_meta %>% filter(agg_id == temp_agg_id) %>%
                                                   pull(cyc), ", ")))
        y_val_agg <- as_tibble(fread(file.path(config$project_dir, config$result_analysis_location,
                                               sprintf("an%06d_lib%03d_agg%02d_type02.csv", temp_tbl_y %>% pull(an_id), temp_lib_id, temp_agg_id)))) %>%
          pull(value)
        col_val_agg <- as_tibble(fread(file.path(config$project_dir, config$result_analysis_location,
                                                 sprintf("an%06d_lib%03d_agg%02d_type02.csv", temp_tbl_col %>% pull(an_id), temp_lib_id, temp_agg_id)))) %>%
          pull(value)
        temp_out <- temp_agg_meta %>%
          bind_cols(tibble(targ_ct = targ_ct %>% pull(value),
                           ctrl_ct = ctrl_ct %>% pull(value))) %>%
          group_by_at(cycs) %>%
          summarise(targ_ct = sum(targ_ct),
                    ctrl_ct = sum(ctrl_ct),
                    .groups = "keep") %>%
          ungroup() %>%
          mutate(y_val = y_val_agg,
                 col_val = col_val_agg,
                 x_lab = input$input_vis_jitterplot_nickname)
        if(nrow(rv$rv_vis_jitterplot_vis_tbl[[temp_agg_id]]) == 0) {
          temp_out <- temp_out %>%
            mutate(x_val = 1)
        } else {
          temp_out <- temp_out %>%
            mutate(x_val = max(rv$rv_vis_jitterplot_vis_tbl[[temp_agg_id]]$x_val) + 1)
        }
      }
        
      rv$rv_vis_jitterplot_vis_tbl[[temp_agg_id]] <- rv$rv_vis_jitterplot_vis_tbl[[temp_agg_id]] %>%
        bind_rows(temp_out)
    }
    showNotification("Plot updated")
  })
  
  observe({
    req("x_lab" %in% names(rv$rv_vis_jitterplot_vis_tbl[[1]]))
    updateSelectInput(session, "input_vis_jitterplot_remove_entry_choice", 
                      choices = unique(rv$rv_vis_jitterplot_vis_tbl[[1]]$x_lab))
  })
  
  observeEvent(input$input_vis_jitterplot_remove_entry, {
    agg_ids <- as_tibble(fread(file.path(config$project_dir, config$aggregation_metadata))) %>% pull(agg_id)
    for(temp_agg_id in agg_ids) {
      rv$rv_vis_jitterplot_vis_tbl[[temp_agg_id]] <- rv$rv_vis_jitterplot_vis_tbl[[temp_agg_id]] %>%
        filter(x_lab != input$input_vis_jitterplot_remove_entry_choice) %>%
        mutate(x_val <- as.numeric(as.factor(x_val)))
    }
    showNotification("Plot updated")
  })
  
  observeEvent(input$input_vis_jitterplot_remove_all, {
    rv$rv_vis_jitterplot_vis_tbl <- rep(list(tibble()), length(agg_lst))
    updateSelectInput(session, "input_vis_jitterplot_remove_entry_choice", 
                      choices = unique(rv$rv_vis_jitterplot_vis_tbl[[1]]$x_lab))
    showNotification("Plot updated")
  })
  
  observe({
    req("col_val" %in% names(rv$rv_vis_jitterplot_vis_tbl[[1]]))
    temp_agg_id <- as.numeric(unlist(str_split(input$vis_jitterplot_agg_id, ": "))[1])
    updateSliderInput(session, "vis_jitterplot_selectivity_cutoff", 
                      max = max(rv$rv_vis_jitterplot_vis_tbl[[temp_agg_id]]$col_val, na.rm = T) %>% ceiling_dec(1))
  })
  
  output$output_vis_jitterplot <- renderPlotly({
    req(rv$rv_vis_jitterplot_vis_tbl)
    # plotly_empty()
    temp_agg_id <- as.numeric(unlist(str_split(input$vis_jitterplot_agg_id, ": "))[1])
    temp <- rv$rv_vis_jitterplot_vis_tbl[[temp_agg_id]]
    if(nrow(temp) > 0) {
      temp <- rv$rv_vis_jitterplot_vis_tbl[[temp_agg_id]] %>% 
        filter(col_val > input$vis_jitterplot_selectivity_cutoff) %>%
        top_n(input$vis_jitterplot_num_points, y_val)
      missing_cycs <- paste0("cycle", unlist(str_split(agg_meta %>% filter(agg_id == max(agg_id)) %>% pull(cyc), ", "))) %>%
        setdiff(colnames(temp))
      for(cyc in missing_cycs) {
        temp[cyc] <- "?"
      }
      temp_labels <- temp[,order(colnames(temp))] %>%
        select(matches("cycle")) %>%
        apply(1, function(x) {paste0(x, collapse = ", ")})
      temp_labels <- paste0("(", temp_labels, ")")
      temp_samp_lst <- c()
      temp_labels <- paste0(temp_labels, "\n",
                            "Control counts: ", temp$ctrl_ct, "\n",
                            "Target counts: ", temp$targ_ct, "\n",
                            "Selectivity (lb): ", round(temp$col_val, 2))
      temp <- temp %>%
        mutate(hover_text = temp_labels)
      plot_ly(temp, 
              source = "plotly_jitterplot", 
              x = ~jitter(x_val,0.5), 
              y = ~y_val, 
              type = "scattergl", 
              mode = "markers",
              # colors = c("#132B43", "#56B1F7"),
              hoverinfo = "text", 
              text = ~hover_text,
              color = ~col_val,
              marker = list(
                size = ~col_val,
                sizeref = max(temp$col_val, na.rm = T) / 16
              )) %>%
        hide_colorbar() %>%
        layout(xaxis = list(title = "Targets",
                            nticks = length(unique(temp$x_val)),
                            tickvals = unique(temp$x_val),
                            ticktext = unique(temp$x_lab)),
               yaxis = list(title = "Enrichment ratio (lb)"),
               hoverlabel = list(bgcolor = "#deebf7",
                                 align = "left"),
               margin = list(l = 50, 
                             r = 0, 
                             b = 50, 
                             t = 0))
    }
  })

  output$output_vis_jitterplot_save_plot <- downloadHandler("vis_jitterplot.html",
                                                            content = function(file) {
                                                              temp_agg_id <- as.numeric(unlist(str_split(input$vis_jitterplot_agg_id, ": "))[1])
                                                              req(rv$rv_vis_jitterplot_vis_tbl[[temp_agg_id]])
                                                              
                                                              temp <- rv$rv_vis_jitterplot_vis_tbl[[temp_agg_id]]
                                                              if(nrow(temp) > 0) {
                                                                temp <- rv$rv_vis_jitterplot_vis_tbl[[temp_agg_id]] %>% 
                                                                  top_n(input$vis_jitterplot_num_points, y_val)
                                                                missing_cycs <- paste0("cycle", unlist(str_split(agg_meta %>% filter(agg_id == max(agg_id)) %>% pull(cyc), ", "))) %>%
                                                                  setdiff(colnames(temp))
                                                                for(cyc in missing_cycs) {
                                                                  temp[cyc] <- "?"
                                                                }
                                                                temp_labels <- temp[,order(colnames(temp))] %>%
                                                                  select(matches("cycle")) %>%
                                                                  apply(1, function(x) {paste0(x, collapse = ", ")})
                                                                temp_labels <- paste0("(", temp_labels, ")")
                                                                temp_samp_lst <- c()
                                                                temp_labels <- paste0(temp_labels, "\n",
                                                                                      "Controls: ", temp$ctrl_ct, "\n",
                                                                                      "Target counts: ", temp$targ_ct, "\n",
                                                                                      "Selectivity (lb): ", round(temp$col_val, 2))
                                                                temp <- temp %>%
                                                                  mutate(hover_text = temp_labels)
                                                                temp_plot <- plot_ly(temp, 
                                                                                     x = ~jitter(x_val,0.5), 
                                                                                     y = ~y_val, 
                                                                                     type = "scattergl", 
                                                                                     mode = "markers",
                                                                                     hoverinfo = "text", 
                                                                                     text = ~hover_text,
                                                                                     size = ~col_val, 
                                                                                     color = ~col_val) %>%
                                                                  hide_colorbar() %>%
                                                                  layout(xaxis = list(title = "Targets",
                                                                                      nticks = length(unique(temp$x_val)),
                                                                                      tickvals = unique(temp$x_val),
                                                                                      ticktext = unique(temp$x_lab)),
                                                                         yaxis = list(title = "Enrichment ratio (lb)"),
                                                                         hoverlabel = list(bgcolor = "#deebf7",
                                                                                           align = "left"),
                                                                         margin = list(l = 50, 
                                                                                       r = 0, 
                                                                                       b = 50, 
                                                                                       t = 0))
                                                              }
                                                              htmlwidgets::saveWidget(partial_bundle(temp_plot), file)
                                                            },
                                                            contentType = "html")
    
  output$output_vis_jitterplot_save_tbl <- downloadHandler("vis_jitterplot_tbl_save.csv",
                                                     content = function(file) {
                                                       temp_agg_id <- as.numeric(unlist(str_split(input$vis_jitterplot_agg_id, ": "))[1])
                                                       write.csv(rv$rv_vis_jitterplot_vis_tbl[[temp_agg_id]], file, row.names = F)
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
        left_join(rv$rv_lib_meta, by = "lib_id") %>%
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
    
    load_cpd_notif <- showNotification("Loading metadata...", duration = NULL)
    if(!as.character(temp_lib_id) %in% names(rv$rv_lib_enum_struct_map)) {
      temp_meta <- file.path(config$project_dir, config$library_metadata_column_location, 
                             sprintf("lib%03d_agg07_meta_cols.csv", temp_lib_id)) %>% 
        fread() %>% 
        as_tibble()
      temp_struct <- file.path(config$project_dir, config$enumerated_structure_location, 
                               sprintf("lib%03d.csv", temp_lib_id)) %>% 
        fread() %>% 
        as_tibble()
      
      temp_meta_filled <- temp_meta %>%
        expand(cycle1, cycle2, cycle3)
      temp <- temp_meta %>%
        select(-lib_id) %>%
        mutate(structure = temp_struct$structure) %>%
        right_join(temp_meta_filled) %>%
        arrange(cycle1, cycle2, cycle3)
      rv$rv_lib_enum_struct_map[[as.character(temp_lib_id)]] <- array(data = temp$structure, 
                                                                      dim = c(max(temp$cycle3), 
                                                                              max(temp$cycle2), 
                                                                              max(temp$cycle1))) %>% aperm()
    }
    
    if(!as.character(temp_lib_id) %in% names(rv$rv_lib_enum_prop_map)) {
      rv$rv_lib_enum_prop_map[[as.character(temp_lib_id)]] <- tibble(
        emw = file.path(config$project_dir, config$enumerated_properties_location, 
                        sprintf("lib%03d_emw.csv", temp_lib_id)) %>% fread() %>% as_tibble() %>% pull(value),
        fsp3 = file.path(config$project_dir, config$enumerated_properties_location, 
                         sprintf("lib%03d_fsp3.csv", temp_lib_id)) %>% fread() %>% as_tibble() %>% pull(value),
        nrb = file.path(config$project_dir, config$enumerated_properties_location, 
                        sprintf("lib%03d_nrb.csv", temp_lib_id)) %>% fread() %>% as_tibble() %>% pull(value),
        slogp = file.path(config$project_dir, config$enumerated_properties_location, 
                          sprintf("lib%03d_slogp.csv", temp_lib_id)) %>% fread() %>% as_tibble() %>% pull(value),
        filter_emw = T,
        filter_fsp3 = T,
        filter_nrb = T,
        filter_slogp = T
      )
    }
    
    temp_min <- floor(min(rv$rv_lib_enum_prop_map[[as.character(temp_lib_id)]]$emw))
    temp_max <- ceiling(max(rv$rv_lib_enum_prop_map[[as.character(temp_lib_id)]]$emw))
    updateSliderInput(session, "input_3d_er_sp_interactive_filter_mw_slider", 
                      min = temp_min, max = temp_max,
                      value = c(temp_min, temp_max))
    temp_min <- min(rv$rv_lib_enum_prop_map[[as.character(temp_lib_id)]]$fsp3) %>% floor_dec(2)
    temp_max <- max(rv$rv_lib_enum_prop_map[[as.character(temp_lib_id)]]$fsp3) %>% ceiling_dec(2)
    updateSliderInput(session, "input_3d_er_sp_interactive_filter_fsp3_slider", 
                      min = temp_min, max = temp_max,
                      value = c(temp_min, temp_max))
    temp_min <- min(rv$rv_lib_enum_prop_map[[as.character(temp_lib_id)]]$nrb)
    temp_max <- max(rv$rv_lib_enum_prop_map[[as.character(temp_lib_id)]]$nrb)
    updateSliderInput(session, "input_3d_er_sp_interactive_filter_nrb_slider", 
                      min = temp_min, max = temp_max,
                      value = c(temp_min, temp_max))
    temp_min <- floor(min(rv$rv_lib_enum_prop_map[[as.character(temp_lib_id)]]$slogp)) %>% floor_dec(2)
    temp_max <- ceiling(max(rv$rv_lib_enum_prop_map[[as.character(temp_lib_id)]]$slogp)) %>% ceiling_dec(2)
    updateSliderInput(session, "input_3d_er_sp_interactive_filter_slogp_slider", 
                      min = temp_min, max = temp_max,
                      value = c(temp_min, temp_max))
    
    removeNotification(load_cpd_notif)
    
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
      mutate(label = paste0("(", cycle1, ", ", cycle2, ", ", cycle3 ,")\n", "er_lb = ", round(value, 2))) %>%
      bind_cols(rv$rv_lib_enum_prop_map[[as.character(temp_lib_id)]]) %>%
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
  
  output$output_3d_er_sp_interactive_filter_ui <- renderUI({
    req(input$input_3d_er_sp_interactive_filter_cycle_id)
    req(input$input_3d_er_interactive_lib_id)
    temp_lib_id <- as.numeric(unlist(str_split(input$input_3d_er_interactive_lib_id, ": "))[1])
    temp_cycle_id <- as.numeric(input$input_3d_er_sp_interactive_filter_cycle_id)
    cycle_choices <- bb_meta %>%
      select(lib_id, cyc, tag_id) %>%
      filter(lib_id == temp_lib_id & cyc == temp_cycle_id)
    cycle_selections <- rv$rv_3d_er_sp_interactive_cycle_selections[[as.character(temp_lib_id)]][[temp_cycle_id]]
    checkboxGroupInput("interacitve_3d_sp_cycle_checkbox", NULL, choices = cycle_choices$tag_id, selected = cycle_selections)
  })
  
  observeEvent(input$interacitve_3d_sp_cycle_checkbox, {
    req(input$input_3d_er_sp_interactive_filter_cycle_id)
    req(input$input_3d_er_interactive_lib_id)
    temp_lib_id <- as.numeric(unlist(str_split(input$input_3d_er_interactive_lib_id, ": "))[1])
    temp_cycle_id <- as.numeric(input$input_3d_er_sp_interactive_filter_cycle_id)
    rv$rv_3d_er_sp_interactive_cycle_selections[[as.character(temp_lib_id)]][[temp_cycle_id]] <- input$interacitve_3d_sp_cycle_checkbox
  })
  
  observeEvent(input$input_3d_er_sp_interactive_filter_reset, {
    rv$rv_3d_er_sp_interactive_cycle_selections <- input_3d_er_sp_interactive_cycle_selections
  })
  
  observeEvent(input$input_3d_er_sp_interactive_filter_mw_slider, {
    req(rv$rv_3d_er_sp_interactive_tbl)
    rv$rv_3d_er_sp_interactive_tbl <- rv$rv_3d_er_sp_interactive_tbl %>%
      mutate(filter_emw = emw >= input$input_3d_er_sp_interactive_filter_mw_slider[1] & emw <= input$input_3d_er_sp_interactive_filter_mw_slider[2])
  })
  
  output$input_3d_er_sp_interactive_filter_mw_plot <- renderPlot({
    req(rv$rv_3d_er_sp_interactive_tbl)
    req(input$input_3d_er_sp_interactive_filter_mw_slider)
    temp <- rv$rv_3d_er_sp_interactive_tbl
    ggplot(temp, aes(x = emw, fill = filter_emw)) +
      geom_histogram(aes(y = ..count..),
                     binwidth = 10) +
      xlab("molecular weight") +
      ylab("") +
      theme_bw() +
      theme(axis.text = element_text(color = "black", size = 10),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.title = element_text(face = "plain", size = 10),
            axis.ticks.length = unit(1, "mm"),
            panel.grid.major.y = element_blank(),
            # panel.grid.minor = element_blank(),
            plot.margin=unit(c(0,0,0,0),"mm")) +
      scale_fill_manual(values = c("light gray", "dark gray")) +
      guides(fill = "none")
  }, height = 200)
  
  observeEvent(input$input_3d_er_sp_interactive_filter_fsp3_slider, {
    req(rv$rv_3d_er_sp_interactive_tbl)
    rv$rv_3d_er_sp_interactive_tbl <- rv$rv_3d_er_sp_interactive_tbl %>%
      mutate(filter_fsp3 = fsp3 >= input$input_3d_er_sp_interactive_filter_fsp3_slider[1] & fsp3 <= input$input_3d_er_sp_interactive_filter_fsp3_slider[2])
  })
  
  output$input_3d_er_sp_interactive_filter_fsp3_plot <- renderPlot({
    req(rv$rv_3d_er_sp_interactive_tbl)
    req(input$input_3d_er_sp_interactive_filter_fsp3_slider)
    temp <- rv$rv_3d_er_sp_interactive_tbl
    ggplot(temp, aes(x = fsp3, fill = filter_fsp3)) +
      geom_histogram(aes(y = ..count..),
                     binwidth = 0.01) +
      xlab("fraction sp3") +
      ylab("") +
      theme_bw() +
      theme(axis.text = element_text(color = "black", size = 10),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.title = element_text(face = "plain", size = 10),
            axis.ticks.length = unit(1, "mm"),
            panel.grid.major.y = element_blank(),
            # panel.grid.minor = element_blank(),
            plot.margin=unit(c(0,0,0,0),"mm")) +
      scale_fill_manual(values = c("light gray", "dark gray")) +
      guides(fill = "none")
  }, height = 200)
  
  observeEvent(input$input_3d_er_sp_interactive_filter_nrb_slider, {
    req(rv$rv_3d_er_sp_interactive_tbl)
    rv$rv_3d_er_sp_interactive_tbl <- rv$rv_3d_er_sp_interactive_tbl %>%
      mutate(filter_nrb = nrb >= input$input_3d_er_sp_interactive_filter_nrb_slider[1] & nrb <= input$input_3d_er_sp_interactive_filter_nrb_slider[2])
  })
  
  output$input_3d_er_sp_interactive_filter_nrb_plot <- renderPlot({
    req(rv$rv_3d_er_sp_interactive_tbl)
    req(input$input_3d_er_sp_interactive_filter_nrb_slider)
    temp <- rv$rv_3d_er_sp_interactive_tbl
    ggplot(temp, aes(x = nrb, fill = filter_nrb)) +
      geom_histogram(aes(y = ..count..),
                     binwidth = 1) +
      xlab("num rotatable bonds") +
      ylab("") +
      theme_bw() +
      theme(axis.text = element_text(color = "black", size = 10),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.title = element_text(face = "plain", size = 10),
            axis.ticks.length = unit(1, "mm"),
            panel.grid.major.y = element_blank(),
            # panel.grid.minor = element_blank(),
            plot.margin=unit(c(0,0,0,0),"mm")) +
      scale_fill_manual(values = c("light gray", "dark gray")) +
      guides(fill = "none")
  }, height = 200)
  
  observeEvent(input$input_3d_er_sp_interactive_filter_slogp_slider, {
    req(rv$rv_3d_er_sp_interactive_tbl)
    rv$rv_3d_er_sp_interactive_tbl <- rv$rv_3d_er_sp_interactive_tbl %>%
      mutate(filter_slogp = slogp >= input$input_3d_er_sp_interactive_filter_slogp_slider[1] & slogp <= input$input_3d_er_sp_interactive_filter_slogp_slider[2])
  })
  
  output$input_3d_er_sp_interactive_filter_slogp_plot <- renderPlot({
    req(rv$rv_3d_er_sp_interactive_tbl)
    req(input$input_3d_er_sp_interactive_filter_slogp_slider)
    temp <- rv$rv_3d_er_sp_interactive_tbl
    ggplot(temp, aes(x = slogp, fill = filter_slogp)) +
      geom_histogram(aes(y = ..count..),
                     binwidth = 0.01) +
      xlab("calculated logP") +
      ylab("") +
      theme_bw() +
      theme(axis.text = element_text(color = "black", size = 10),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.title = element_text(face = "plain", size = 10),
            axis.ticks.length = unit(1, "mm"),
            panel.grid.major.y = element_blank(),
            # panel.grid.minor = element_blank(),
            plot.margin=unit(c(0,0,0,0),"mm")) +
      scale_fill_manual(values = c("light gray", "dark gray")) +
      guides(fill = "none")
  }, height = 200)
  
  output$output_3d_er_sp_interactive <- renderPlotly({
    req(rv$rv_3d_er_sp_interactive_tbl)
    req(rv$rv_3d_er_sp_interactive_tbl$scaled_size)
    temp_lib_id <- unlist(str_split(input$input_3d_er_interactive_lib_id, ": "))[1]
    plot_ly(rv$rv_3d_er_sp_interactive_tbl %>%
              filter(cycle1 %in% rv$rv_3d_er_sp_interactive_cycle_selections[[temp_lib_id]][[1]] &
                       cycle2 %in% rv$rv_3d_er_sp_interactive_cycle_selections[[temp_lib_id]][[2]] &
                       cycle3 %in% rv$rv_3d_er_sp_interactive_cycle_selections[[temp_lib_id]][[3]]) %>%
              filter(filter_emw & filter_fsp3 & filter_nrb & filter_slogp) %>%
              top_n(input$visualize_3d_er_sp_interactive_num_points) %>%
              mutate(x = rv$rv_3d_er_sp_interactive_x_map[as.character(cycle1)],
                     y = rv$rv_3d_er_sp_interactive_y_map[as.character(cycle2)],
                     z = rv$rv_3d_er_sp_interactive_z_map[as.character(cycle3)]), 
            x = ~x, y = ~y, z = ~z, text = ~label,
            hoverinfo = "text",
            # height = input$visualize_3d_er_sp_interactive_size,
            height = input$main_panel_innerWidth_3d_er_sp_interactive*0.35,
            marker = list(
              color = ~value, 
              size = ~scaled_size, 
              colorscale = "Reds",
              showscale = F,
              line = list(width = 0.5, color = "#666666")
            ),
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
                          aspectmode = "cube"),
             margin = list(l = 0, r = 10, b = 0, t = 10)) %>%
      # event_register("plotly_hover") %>%
      # event_register("plotly_unhover") %>%
      event_register("plotly_click")
  })
  
  # plotly_3d_er_interactive_hover_event <- reactive({
  #   req(rv$rv_3d_er_sp_interactive_tbl)
  #   event_data(event = "plotly_hover", source = "plotly_3d_er_interactive")
  # })
  # 
  # plotly_3d_er_interactive_unhover_event <- reactive({
  #   req(rv$rv_3d_er_sp_interactive_tbl)
  #   event_data(event = "plotly_unhover", source = "plotly_3d_er_interactive")
  # })
  
  plotly_3d_er_interactive_click_event <- reactive({
    req(rv$rv_3d_er_sp_interactive_tbl)
    event_data(event = "plotly_click", source = "plotly_3d_er_interactive")
  })
  
  observeEvent(plotly_3d_er_interactive_click_event(), {
    req(rv$rv_lib_enum_struct_map)
    temp_lib_id <- as.numeric(unlist(str_split(input$input_3d_er_interactive_lib_id, ": "))[1])
    cy1 <- which(rv$rv_3d_er_sp_interactive_x_map == plotly_3d_er_interactive_click_event()$x) %>% names() %>% as.numeric()
    cy2 <- which(rv$rv_3d_er_sp_interactive_y_map == plotly_3d_er_interactive_click_event()$y) %>% names() %>% as.numeric()
    cy3 <- which(rv$rv_3d_er_sp_interactive_z_map == plotly_3d_er_interactive_click_event()$z) %>% names() %>% as.numeric()
    rv$rv_3d_er_sp_interactive_selected_struct <- rv$rv_lib_enum_struct_map[[as.character(temp_lib_id)]][cy1, cy2, cy3]
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
  height = reactive(ifelse(!is.null(input$side_panel_innerWidth_3d_er_sp_interactive),input$side_panel_innerWidth_3d_er_sp_interactive*0.25-120,0)),
  width = reactive(ifelse(!is.null(input$side_panel_innerWidth_3d_er_sp_interactive),input$side_panel_innerWidth_3d_er_sp_interactive*0.25-120,0))
  )
  
  output$output_3d_er_sp_interactive_selected_struct_props <- renderText({
    cy1 <- which(rv$rv_3d_er_sp_interactive_x_map == plotly_3d_er_interactive_click_event()$x) %>% names() %>% as.numeric()
    cy2 <- which(rv$rv_3d_er_sp_interactive_y_map == plotly_3d_er_interactive_click_event()$y) %>% names() %>% as.numeric()
    cy3 <- which(rv$rv_3d_er_sp_interactive_z_map == plotly_3d_er_interactive_click_event()$z) %>% names() %>% as.numeric()
    if(length(cy1) == 1 & length(cy2) == 1 & length(cy3) == 1) {
      temp <- rv$rv_3d_er_sp_interactive_tbl %>%
        filter(cycle1 == cy1 & cycle2 == cy2 & cycle3 == cy3)
      paste0("Molecular weight: ", temp %>% pull(emw) %>% round(2), "\n",
             "Fraction sp3: ", temp %>% pull(fsp3) %>% round(2), "\n", 
             "Num rotatable bonds: ", temp %>% pull(nrb), "\n", 
             "LogP: ", temp %>% pull(slogp) %>% round(2))
    }
  })
  
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
  
  observe({
    req(input$input_interactive_color_dropdown)
    req(rv$rv_interactive_meta)
    temp_cols <- rep("Other", nrow(rv$rv_interactive_meta))
    if(input$input_interactive_color_dropdown == "Manually") {
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
    else if(input$input_interactive_color_dropdown == "By table") {
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
                height = as.numeric(input$interactive_2d_sp_box_height)-440) %>%
          layout(xaxis = list(title = temp_x_lab),
                 yaxis = list(title = temp_y_lab),
                 hoverlabel = list(align = "left"))
      }
      else if(input$input_interactive_x_log & !input$input_interactive_y_log) {
        if(input$input_interactive_x_type == "enrichment") {
          plot_ly(rv$rv_interactive_tbl, source = "interactive", x = ~x_value, y = ~y_value, color = ~color, colors = c("#132B43", "#56B1F7"),
                  type = "scattergl", mode = "markers", text = ~label,
                  hoverinfo = "text",
                  height = as.numeric(input$interactive_2d_sp_box_height)-440) %>%
            layout(xaxis = list(title = temp_x_lab, type = "log"),
                   yaxis = list(title = temp_y_lab),
                   hoverlabel = list(align = "left"))
        }
        else {
          plot_ly(rv$rv_interactive_tbl, source = "interactive", x = ~x_value, y = ~y_value, color = ~color, colors = c("#132B43", "#56B1F7"),
                  type = "scattergl", mode = "markers", text = ~label,
                  hoverinfo = "text",
                  height = as.numeric(input$interactive_2d_sp_box_height)-440) %>%
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
                  height = as.numeric(input$interactive_2d_sp_box_height)-440) %>%
            layout(xaxis = list(title = temp_x_lab),
                   yaxis = list(title = temp_y_lab, type = "log"),
                   hoverlabel = list(align = "left"))
        }
        else {
          plot_ly(rv$rv_interactive_tbl, source = "interactive", x = ~x_value, y = ~y_value, color = ~color, colors = c("#132B43", "#56B1F7"),
                  type = "scattergl", mode = "markers", text = ~label,
                  hoverinfo = "text",
                  height = as.numeric(input$interactive_2d_sp_box_height)-440) %>%
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
                  height = as.numeric(input$interactive_2d_sp_box_height)-440) %>%
            layout(xaxis = list(title = temp_x_lab, type = "log"),
                   yaxis = list(title = temp_y_lab),
                   hoverlabel = list(align = "left"))
        }
        else if(input$input_interactive_x_type != "enrichment") {
          plot_ly(rv$rv_interactive_tbl, source = "interactive", x = ~x_value, y = ~y_value, color = ~color, colors = c("#132B43", "#56B1F7"),
                  type = "scattergl", mode = "markers", text = ~label,
                  hoverinfo = "text",
                  height = as.numeric(input$interactive_2d_sp_box_height)-440) %>%
            layout(xaxis = list(title = temp_x_lab),
                   yaxis = list(title = temp_y_lab, type = "log"),
                   hoverlabel = list(align = "left"))
        }
        else {
          plot_ly(rv$rv_interactive_tbl, source = "interactive", x = ~x_value, y = ~y_value, color = ~color, colors = c("#132B43", "#56B1F7"),
                  type = "scattergl", mode = "markers", text = ~label,
                  hoverinfo = "text",
                  height = as.numeric(input$interactive_2d_sp_box_height)-440) %>%
            layout(xaxis = list(title = temp_x_lab, type = "log"),
                   yaxis = list(title = temp_y_lab, type = "log"),
                   hoverlabel = list(align = "left"))
        }
      }
    }
  })

  output$output_interactive_sp_save_plot_html <- downloadHandler("plot_name.html", 
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
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  observeEvent(input$input_vis_mono_heatmap_completed_analyses, {
    req(input$input_vis_mono_heatmap_completed_analyses)
    if(input$input_vis_mono_heatmap_completed_analyses == "All") {
      completed_res <- list.files(file.path(config$project_dir, config$result_metadata_location))
      temp <- tibble()
      for(i in completed_res) {
        temp <- temp %>%
          bind_rows(as_tibble(fread(file.path(config$project_dir, config$result_metadata_location, i))) %>% 
                      mutate(res_name = str_remove(i, ".csv")))
      }
      rv$rv_vis_mono_heatmap_completed_an_tbl <- temp
    }
    else {
      rv$rv_vis_mono_heatmap_completed_an_tbl <- file.path(config$project_dir, config$result_metadata_location, paste0(input$input_vis_mono_heatmap_completed_analyses, ".csv")) %>%
        fread() %>%
        as_tibble() %>%
        mutate(res_name = input$input_vis_mono_heatmap_completed_analyses)
    }
  })
  
  observe({
    updateSelectInput(session, "input_vis_mono_heatmap_completed_analyses", choices = c("All", rv$rv_res_lst))
  })
  
  output$vis_mono_heatmap_completed_an_tbl <- DT::renderDataTable({
    req(rv$rv_vis_mono_heatmap_completed_an_tbl)
    temp_lib_id <- sapply(strsplit(rv$rv_vis_mono_heatmap_completed_an_tbl$targ_ct_lst, "_lib|.csv"), `[`, 2) %>%
      as.numeric()
    DT::datatable(rv$rv_vis_mono_heatmap_completed_an_tbl %>%
                    mutate(lib_id = sprintf("%03d", temp_lib_id)) %>%
                    relocate(lib_id) %>%
                    select(-agg_id, -an_type_id, -an_id, -targ_ct_lst, -ctrl_ct_lst),
                  rownames = F, filter = "top", selection = "single",
                  options = list(dom = 'tp',
                                 pageLength = 4),
                  class = "nowrap display")
  })
  
  vis_mono_heatmap_an_tbl_proxy <- dataTableProxy("vis_mono_heatmap_completed_an_tbl")

  observeEvent(input$input_vis_mono_heatmap_add_entry, {
    if(nrow(rv$rv_vis_mono_heatmap_vis_tbl[[1]] > 0)){
      if(input$input_vis_mono_heatmap_nickname %in% rv$rv_vis_mono_heatmap_vis_tbl[[1]]$y_lab) {
        showNotification("Nickname already used. Choose another.", type = "warning")
        req(!input$input_vis_mono_heatmap_nickname %in% rv$rv_vis_mono_heatmap_vis_tbl[[1]]$y_lab)
      }
    }
    temp_tbl <- rv$rv_vis_mono_heatmap_completed_an_tbl[input$vis_mono_heatmap_completed_an_tbl_rows_selected,]
    temp_lib_id <- sapply(strsplit(temp_tbl$targ_ct_lst, "_lib|.csv"), `[`, 2) %>%
      as.numeric()
    temp_an_id <- temp_tbl %>%
      pull(an_id)
    
    for(temp_agg_id in 1:3) {
      temp_agg_meta <- as_tibble(fread(file.path(config$project_dir, config$library_metadata_column_location,
                                                 sprintf("lib%03d_agg%02d_meta_cols.csv", temp_lib_id, temp_agg_id)))) %>%
        rename(x_val = 2) %>%
        select(x_val)
      temp_out <- as_tibble(fread(file.path(config$project_dir, config$result_analysis_location,
                                            sprintf("an%06d_lib%03d_agg%02d_type02.csv", temp_an_id, temp_lib_id, temp_agg_id)))) %>%
        bind_cols(temp_agg_meta) %>%
        mutate(x_lab = paste("Cycle", temp_agg_id),
               y_lab = input$input_vis_mono_heatmap_nickname)
      if(nrow(rv$rv_vis_mono_heatmap_vis_tbl[[temp_agg_id]]) == 0) {
        temp_out <- temp_out %>%
          mutate(y_val = 1)
      } else {
        temp_out <- temp_out %>%
          mutate(y_val = max(rv$rv_vis_mono_heatmap_vis_tbl[[temp_agg_id]]$y_val) + 1)
      }
      
      rv$rv_vis_mono_heatmap_vis_tbl[[temp_agg_id]] <- rv$rv_vis_mono_heatmap_vis_tbl[[temp_agg_id]] %>%
        bind_rows(temp_out)
    }
    showNotification("Plot updated")
  })
  
  observe({
    req("y_lab" %in% names(rv$rv_vis_mono_heatmap_vis_tbl[[1]]))
    updateSelectInput(session, "input_vis_mono_heatmap_remove_entry_choice", 
                      choices = unique(rv$rv_vis_mono_heatmap_vis_tbl[[1]]$y_lab))
  })
  
  observeEvent(input$input_vis_mono_heatmap_remove_entry, {
    for(temp_agg_id in 1:3) {
      rv$rv_vis_mono_heatmap_vis_tbl[[temp_agg_id]] <- rv$rv_vis_mono_heatmap_vis_tbl[[temp_agg_id]] %>%
        filter(y_lab != input$input_vis_mono_heatmap_remove_entry_choice) %>%
        mutate(y_val <- as.numeric(as.factor(y_val)))
    }
    showNotification("Plot updated")
  })
  
  observeEvent(input$input_vis_mono_heatmap_remove_all, {
    rv$rv_vis_mono_heatmap_vis_tbl <- rep(list(tibble()), length(agg_lst))
    updateSelectInput(session, "input_vis_mono_heatmap_remove_entry_choice", 
                      choices = unique(rv$rv_vis_mono_heatmap_vis_tbl[[1]]$y_lab))
    showNotification("Plot updated")
  })
  
  output$output_vis_mono_heatmap <- renderPlotly({
    req(rv$rv_vis_mono_heatmap_vis_tbl)
    temp_agg_id <- as.numeric(unlist(str_split(input$vis_mono_heatmap_agg_id, ": "))[1])
    temp <- rv$rv_vis_mono_heatmap_vis_tbl[[temp_agg_id]]
    if(nrow(temp) > 0) {
      temp <- temp %>%
        mutate(log2_value = log2(value))
      temp$log2_value[temp$log2_value < input$input_vis_mono_heatmap_col_scale * -1] <- input$input_vis_mono_heatmap_col_scale * -1
      temp$log2_value[temp$log2_value > input$input_vis_mono_heatmap_col_scale] <- input$input_vis_mono_heatmap_col_scale
      temp_labels <- paste0(temp$y_lab, "\n",
                            "Cycle ID: ", temp$x_val, "\n",
                            "Enrichment: ", round(temp$value, 2), "\n")
      temp <- temp %>%
        mutate(hover_text = temp_labels)
      
      col3 <- colorRamp(c("red", "white", "blue"))
      colorlength <- 100
      
      null_value <- (0 - min(temp$log2_value, na.rm = T)) / (max(temp$log2_value, na.rm = T) - min(temp$log2_value, na.rm = T))
      border <- as.integer(null_value * colorlength)
      
      colorscale <- as.list(1:colorlength)
      
      # colorscale below zero
      s <- scales::seq_gradient_pal("#2c7bb6", "#FFFFFF", "Lab")(seq(0, 1, length.out = border))
      for(i in 1:border) {
        colorscale[[i]] <- c((i - 1) / colorlength, s[i])
      }
      
      # colorscale above zero
      s <- scales::seq_gradient_pal("#FFFFFF", "#d7191c", "Lab")(seq(0, 1, length.out = colorlength - border))
      for(i in 1:(colorlength - border)) {
        colorscale[[i + border]] <- c((i + border) / colorlength, s[i])
      }

      plot_ly(temp, 
              source = "plotly_mono_heatmap", 
              x = ~x_val, 
              y = ~y_val, 
              z = ~log2_value,
              type = "heatmap",
              colorscale = colorscale,
              hoverinfo = "text",
              text = ~hover_text) %>%
        colorbar(title = "log2(enrichment)") %>%
        layout(xaxis = list(title = "Cycle ID"),
               yaxis = list(title = "Targets",
                            nticks = length(unique(temp$y_val)),
                            tickvals = unique(temp$y_val),
                            ticktext = unique(temp$y_lab)),
               hoverlabel = list(bgcolor = "#deebf7",
                                 align = "left"),
               margin = list(l = 50, 
                             r = 0, 
                             b = 50, 
                             t = 0))
    }
  })
  
  output$output_vis_mono_heatmap_save_plot <- downloadHandler("vis_mono_heatmap.html",
                                                            content = function(file) {
                                                              temp_agg_id <- as.numeric(unlist(str_split(input$vis_mono_heatmap_agg_id, ": "))[1])
                                                              req(rv$rv_vis_mono_heatmap_vis_tbl[[temp_agg_id]])
                                                              
                                                              temp <- rv$rv_vis_mono_heatmap_vis_tbl[[temp_agg_id]]
                                                              if(nrow(temp) > 0) {
                                                                temp <- temp %>%
                                                                  mutate(log2_value = log2(value))
                                                                temp$log2_value[temp$log2_value < input$input_vis_mono_heatmap_col_scale * -1] <- input$input_vis_mono_heatmap_col_scale * -1
                                                                temp$log2_value[temp$log2_value > input$input_vis_mono_heatmap_col_scale] <- input$input_vis_mono_heatmap_col_scale
                                                                temp_labels <- paste0(temp$y_lab, "\n",
                                                                                      "Cycle ID: ", temp$x_val, "\n",
                                                                                      "Enrichment: ", round(temp$value, 2), "\n")
                                                                temp <- temp %>%
                                                                  mutate(hover_text = temp_labels)
                                                                
                                                                col3 <- colorRamp(c("red", "white", "blue"))
                                                                colorlength <- 100
                                                                
                                                                null_value <- (0 - min(temp$log2_value, na.rm = T)) / (max(temp$log2_value, na.rm = T) - min(temp$log2_value, na.rm = T))
                                                                border <- as.integer(null_value * colorlength)
                                                                
                                                                colorscale <- as.list(1:colorlength)
                                                                
                                                                # colorscale below zero
                                                                s <- scales::seq_gradient_pal("#2c7bb6", "#FFFFFF", "Lab")(seq(0, 1, length.out = border))
                                                                for(i in 1:border) {
                                                                  colorscale[[i]] <- c((i - 1) / colorlength, s[i])
                                                                }
                                                                
                                                                # colorscale above zero
                                                                s <- scales::seq_gradient_pal("#FFFFFF", "#d7191c", "Lab")(seq(0, 1, length.out = colorlength - border))
                                                                for(i in 1:(colorlength - border)) {
                                                                  colorscale[[i + border]] <- c((i + border) / colorlength, s[i])
                                                                }
                                                                
                                                                temp_plot <- plot_ly(temp, 
                                                                        source = "plotly_mono_heatmap", 
                                                                        x = ~x_val, 
                                                                        y = ~y_val, 
                                                                        z = ~log2_value,
                                                                        type = "heatmap",
                                                                        colorscale = colorscale,
                                                                        hoverinfo = "text",
                                                                        text = ~hover_text) %>%
                                                                  colorbar(title = "log2(enrichment)") %>%
                                                                  layout(xaxis = list(title = "Cycle ID"),
                                                                         yaxis = list(title = "Targets",
                                                                                      nticks = length(unique(temp$y_val)),
                                                                                      tickvals = unique(temp$y_val),
                                                                                      ticktext = unique(temp$y_lab)),
                                                                         hoverlabel = list(bgcolor = "#deebf7",
                                                                                           align = "left"),
                                                                         margin = list(l = 50, 
                                                                                       r = 0, 
                                                                                       b = 50, 
                                                                                       t = 0))
                                                              }
                                                              htmlwidgets::saveWidget(partial_bundle(temp_plot), file)
                                                            },
                                                            contentType = "html")
  
  output$output_vis_mono_heatmap_save_tbl <- downloadHandler("vis_mono_heatmap_tbl_save.csv",
                                                           content = function(file) {
                                                             temp_agg_id <- as.numeric(unlist(str_split(input$vis_mono_heatmap_agg_id, ": "))[1])
                                                             write.csv(rv$rv_vis_mono_heatmap_vis_tbl[[temp_agg_id]], file, row.names = F)
                                                           },
                                                           contentType = "text/csv")
  
}

shinyApp(ui, server)
  