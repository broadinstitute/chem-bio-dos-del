# Set working directory ####
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
# setwd(utils::getSrcDirectory()[1])
setwd("../..")

# Import packages ####
packages <- c("data.table", "tidyverse", "plotly")
package_check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = T)) {
      install.packages(x, dependencies = T)
      library(x, character.only = T)
    }
  }
)
# suppressPackageStartupMessages(lapply(packages, library, character.only = T))

# Import functions ####
source("./code/fun/functions.R")

# Import metadata ####
source("./code/fun/metadata.R")

# save meta_cols for all libraries (7/8/2021) ####
for(x in lib_meta$lib_id) {
  temp_enum <- lib_enum %>%
    filter(lib_id == x) %>%
    select(lib_id, cycle1, cycle2, cycle3)
  for(i in 1:nrow(agg_meta)) {
    agg_id <- agg_meta$agg_id[i]
    cycs <- paste0("cycle", str_split(agg_meta$cyc[i], ", ")[[1]])
    temp_enum %>%
      select(c(lib_id, cycs)) %>%
      distinct() %>%
      write.csv(paste0("./meta/lib/meta_cols/lib", sprintf("%03d", x), "_agg", sprintf("%02d", agg_id), "_meta_cols.csv"), row.names = F)
  }
}

# outdated below
# for regular libraries
for(x in lib_meta$lib_id) {
  if(!x %in% c(200,201)) {
    temp_enum <- crossing(
      x,
      bb_meta %>% filter(lib_id == x & cyc == 1) %>% pull(tag_id),
      bb_meta %>% filter(lib_id == x & cyc == 2) %>% pull(tag_id),
      bb_meta %>% filter(lib_id == x & cyc == 3) %>% pull(tag_id)
    ) %>%
      rename(lib_id = 1, cycle1 = 2, cycle2 = 3, cycle3 = 4)
    for(i in 1:nrow(agg_meta)) {
      agg_id <- agg_meta$agg_id[i]
      cycs <- paste0("cycle", str_split(agg_meta$cyc[i], ", ")[[1]])
      temp_enum %>%
        select(c(lib_id, cycs)) %>%
        distinct() %>%
        write.csv(paste0("./meta/lib/meta_cols/lib", sprintf("%03d", x), "_agg", sprintf("%02d", agg_id), "_meta_cols.csv"), row.names = F)
    }
  }
}
# for libraries with spike-ins
for(x in lib_meta$lib_id) {
  if(x %in% c(200,201)) {
    cy1 <- bb_meta %>% filter(lib_id == x & cyc == 1) %>% pull(tag_id) %>% head(-1)
    cy2 <- bb_meta %>% filter(lib_id == x & cyc == 2) %>% pull(tag_id) %>% head(-1)
    cy3 <- bb_meta %>% filter(lib_id == x & cyc == 3) %>% pull(tag_id) %>% head(-1)
    temp_enum <- crossing(
      x,
      cy1,
      cy2,
      cy3
    ) %>%
      rename(lib_id = 1, cycle1 = 2, cycle2 = 3, cycle3 = 4)
    for(i in 1:nrow(agg_meta)) {
      agg_id <- agg_meta$agg_id[i]
      cycs <- paste0("cycle", str_split(agg_meta$cyc[i], ", ")[[1]])
      temp_enum %>%
        select(c(lib_id, cycs)) %>%
        distinct() %>%
        write.csv(paste0("./meta/lib/meta_cols/lib", sprintf("%03d", x), "_agg", sprintf("%02d", agg_id), "_meta_cols.csv"), row.names = F)
    }
  }
}

# generate separate counts files for prototype folder structure - protein/bead titration (7/9/2021) ####
ex_run_id <- 8
data_file_name <- "titration_v1v2_hiseq.csv"

for(i in 1) {
  temp <- as_tibble(fread(paste0("./data/count/", data_file_name))) %>%
    select(-matches("^V1|unnamed|cpd|cycle|_tot"))
  for(x in unique(temp$library_id)) {
    for(j in 1:24) {
      temp %>%
        filter(library_id == x) %>%
        select(-library_id) %>%
        select(all_of(j)) %>%
        rename(value = 1) %>%
        write.csv(paste0("./data/count/run", sprintf("%03d", ex_run_id), 
                         "_samp", sprintf("%06d", j), 
                         "_lib", sprintf("%03d", x),
                         ".csv"), row.names = F)
    }
  }
}

# protein/bead titration - first pass analysis (7/16/2021) ####
ptm <- proc.time()
report_name <- format(Sys.time(), "%Y%m%d_%H%M%S")
report_lines <- c()

meta <- run_meta %>%
  left_join(ct_file_meta) %>%
  left_join(samp_meta) %>%
  left_join(targ_meta)

# input from user - lists for tgts and ctrls
targ_ct_lst <- paste0("run008_samp", sprintf("%06d", 7:24), "_lib100.csv")
ctrl_ct_lst <- paste0("run008_samp", sprintf("%06d", 3:4), "_lib100.csv")
res_name <- "res0002"

submeta <- meta %>%
  filter(ct_filename %in% c(targ_ct_lst, ctrl_ct_lst)) %>%
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
        id1 <- ids$ct_filename[i]
        id2 <- ids$ct_filename[j]
        generate_replicate_scatterplots(
          id1,
          id2,
          plot_title = "replicate scatterplot", 
          width = 6, # in inches
          height = 6, # in inches
          out_loc = paste0("./result/", res_name, "/plot")
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

curr_targ_ct_lst <- c()
for(gid in gids) {
  tids <- submeta %>%
    filter(group_id == gid) %>%
    pull(ct_filename)
  curr_targ_ct_lst <- c(curr_targ_ct_lst, paste0(tids, collapse = ", "))
}

in_tbl <- crossing(
  curr_targ_ct_lst,
  paste0(ctrl_ct_lst, collapse = ", "),
  paste0(agg_meta$agg_id, collapse = ", "),
  paste0(an_type_meta$an_type_id[an_type_meta$an_type_name == "enrichment"], collapse = ", "),
  0.95
) %>%
  setNames(c("targ_ct_lst", "ctrl_ct_lst", "agg_id", "an_type_id", "ci"))

in_tbl <- tibble(
  targ_ct_lst = curr_targ_ct_lst,
  ctrl_ct_lst = paste0(ctrl_ct_lst, collapse = ", "),
  agg_id = paste0(agg_meta$agg_id, collapse = ", "),
  an_type_id = paste0(an_type_meta$an_type_id[an_type_meta$an_type_name == "enrichment"], collapse = ", "),
  ci = 0.95
)

# generate/assign analysis ids and update an_meta
in_tbl <- assign_an_ids(in_tbl)

# perform analyses to get lists of values for plotting
temp_ptm <- proc.time()
prfm_an_expts(in_tbl)
temp_dur <- round((proc.time() - temp_ptm)[[3]],2)
report_lines <- c(report_lines, paste0("Time to perform analysis experiments: ", temp_dur, " seconds"))

# generate enrichment scatterplots, y:lb, x:rank of ctrl counts
temp_ptm <- proc.time()
dir.create(paste0("./result/", res_name, "/plot/enrichment_scatterplot_labeled"), recursive = T)
generate_enrichment_scatterplots_labeled(in_tbl, res_name)
temp_dur <- round((proc.time() - temp_ptm)[[3]],2)
report_lines <- c(report_lines, paste0("Time to generate enrichment scatter plots: ", temp_dur, " seconds"))

# generate interactive 3D enrichment scatterplots
temp_ptm <- proc.time()
dir.create(paste0("./result/", res_name, "/plot/enrichment_scatterplot_3d"), recursive = T)
generate_enrichment_scatterplots_3d(in_tbl, res_name)
temp_dur <- round((proc.time() - temp_ptm)[[3]],2)
report_lines <- c(report_lines, paste0("Time to generate interactive 3D enrichment scatterplots: ", temp_dur, " seconds"))

# generate count scatterplots for top trisynthons and save to table
temp_ptm <- proc.time()
dir.create(paste0("./result/", res_name, "/plot/count_scatterplot_top"), recursive = T)
dir.create(paste0("./result/", res_name, "/table"), recursive = T)
generate_ts_ct_scatterplots(in_tbl, res_name)
temp_dur <- round((proc.time() - temp_ptm)[[3]],2)
report_lines <- c(report_lines, paste0("Time to generate count scatter plots for top trisynthons with associated tables: ", temp_dur, " seconds"))

temp_dur <- round((proc.time() - ptm)[[3]],2)
report_lines <- c(report_lines, paste0("Time to complete script: ", temp_dur, " seconds"))

write_lines(report_lines, paste0("./result/", res_name, "/report_", report_name, ".txt"))

# MetaViewer (7/21/2021) ####
library(shiny)
library(shinythemes)
library(DT)
library(jsonlite)

merged_meta <- run_meta %>%
  left_join(ct_file_meta) %>%
  left_join(samp_meta) %>%
  left_join(targ_meta) %>%
  left_join(resin_meta)

for(i in 1) {
  ui <- fluidPage(
    theme = shinytheme("simplex"),
    # themeSelector(),
    titlePanel("MetaViewer 1.0"),
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
          actionButton("reload_lib_meta", "Reload"),
        ),
        conditionalPanel(
          'input.dataset === "comb_lib_meta"',
          div("Metadata associated with the pooled DNA-encoded libraries. Each unique combination of pools, including the ratios, is assigned a unique ID."),
          br(),
          div("Submit this form to add new entries:"),
          numericInput("user_comb_lib_id", "comb_lib_id", ""),
          numericInput("user_lib_id", "lib_id", ""),
          numericInput("user_frac", "frac", ""),
          actionButton("append_comb_lib_meta", "Add entry"),
          actionButton("save_comb_lib_meta", "Save changes"),
          actionButton("reload_comb_lib_meta", "Reload"),
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
          'input.dataset === "targ_meta"',
          div("Metadata associated with the targets."),
          br(),
          div("Submit this form to add new entries:"),
          textInput("user_targ_name", "targ_name", ""),
          actionButton("append_targ_meta", "Add entry"),
          actionButton("save_targ_meta", "Save changes"),
          actionButton("reload_targ_meta", "Reload"),
        ),
        conditionalPanel(
          'input.dataset === "resin_meta"',
          div("Metadata associated with the resins (i.e., beads)."),
          br(),
          div("Submit this form to add new entries:"),
          textInput("user_resin_name", "resin_name", ""),
          actionButton("append_resin_meta", "Add entry"),
          actionButton("save_resin_meta", "Save changes"),
          actionButton("reload_resin_meta", "Reload"),
        ),
        conditionalPanel(
          'input.dataset === "samp_meta"',
          div("Metadata associated with the samples."),
          actionButton("reload_samp_meta", "Reload"),
        ),
        conditionalPanel(
          'input.dataset === "merged_meta"',
          div("Merged metadata to aid in defining analysis experiments."),
          checkboxGroupInput("show_vars", "Columns in merged_meta to show:",
                             names(merged_meta), selected = names(merged_meta)),
          actionButton("reload_merged_meta", "Reload")
        ),
        conditionalPanel(
          'input.dataset === "create_json"',
          div("Create config file (json) from csv table."),
          fileInput("config_upload", "Choose CSV File",
                    multiple = FALSE,
                    accept = c("text/csv",
                               "text/comma-separated-values,text/plain",
                               ".csv")),
          downloadButton("download_json", "Create json")
        )
      ),
      mainPanel(
        tabsetPanel(
          id = 'dataset',
          tabPanel("lib_meta", DT::dataTableOutput("mytable1")),
          tabPanel("comb_lib_meta", DT::dataTableOutput("mytable2")),
          tabPanel("run_meta", DT::dataTableOutput("mytable3")),
          tabPanel("targ_meta", DT::dataTableOutput("mytable4")),
          tabPanel("resin_meta", DT::dataTableOutput("mytable5")),
          tabPanel("samp_meta", DT::dataTableOutput("mytable6")),
          tabPanel("merged_meta", DT::dataTableOutput("mytable7")),
          tabPanel("create_json", DT::dataTableOutput("mytable8"))
        )
      )
    )
  )
}

server <- function(input, output) {
  
  rv <- reactiveValues(rv_lib_meta = lib_meta,
                       rv_comb_lib_meta = comb_lib_meta,
                       rv_run_meta = run_meta,
                       rv_targ_meta = targ_meta,
                       rv_resin_meta = resin_meta,
                       rv_samp_meta = samp_meta,
                       rv_merged_meta = run_meta %>%
                         left_join(ct_file_meta) %>%
                         left_join(samp_meta) %>%
                         left_join(targ_meta) %>%
                         left_join(resin_meta))

  output$mytable1 <- DT::renderDataTable({
    DT::datatable(rv$rv_lib_meta, rownames = F, filter = "top")
  })

  output$mytable2 <- DT::renderDataTable({
    DT::datatable(rv$rv_comb_lib_meta, rownames = F, filter = "top")
  })

  output$mytable3 <- DT::renderDataTable({
    DT::datatable(rv$rv_run_meta, rownames = F, filter = "top")
  })
  
  output$mytable4 <- DT::renderDataTable({
    DT::datatable(rv$rv_targ_meta, rownames = F, filter = "top")
  })
  
  output$mytable5 <- DT::renderDataTable({
    DT::datatable(rv$rv_resin_meta, rownames = F, filter = "top")
  })
  
  output$mytable6 <- DT::renderDataTable({
    DT::datatable(rv$rv_samp_meta, rownames = F, filter = "top")
  })
  
  output$mytable7 <- DT::renderDataTable({
    DT::datatable(rv$rv_merged_meta[, input$show_vars, drop = FALSE],
                  rownames = F, filter = "top")
  })
  
  output$mytable8 <- DT::renderDataTable({
    if(is.null(input$config_upload)) {
      return(NULL)
    }
    DT::datatable(fread(input$config_upload$datapath),
                   rownames = F)
  })

  observeEvent(input$append_lib_meta, {
    temp <- tibble(lib_id = min(setdiff(1:10000, rv$rv_lib_meta$lib_id)),
                   lib_seq = input$user_lib_seq,
                   lib_name = input$user_lib_name,
                   chemist = input$user_chemist,
                   lib_completion_date = input$user_lib_completion_date,
                   lib_loc = input$user_lib_loc)
    rv$rv_lib_meta <- rv$rv_lib_meta %>%
      bind_rows(temp)
  })
  
  observeEvent(input$save_lib_meta, {
    write.csv(rv$rv_lib_meta, "./meta/lib/lib_meta.csv", row.names = F)
  })
  
  observeEvent(input$reload_lib_meta, {
    rv$rv_lib_meta <- as_tibble(fread("./meta/lib/lib_meta.csv"))
  })
  
  observeEvent(input$append_comb_lib_meta, {
    temp <- tibble(comb_lib_id = input$user_comb_lib_id,
                   lib_id = input$user_lib_id,
                   frac = input$user_frac)
    rv$rv_comb_lib_meta <- rv$rv_comb_lib_meta %>%
      bind_rows(temp)
  })
  
  observeEvent(input$save_comb_lib_meta, {
    write.csv(rv$rv_comb_lib_meta, "./meta/lib/comb_lib_meta.csv", row.names = F)
  })
  
  observeEvent(input$reload_comb_lib_meta, {
    rv$rv_comb_lib_meta <- as_tibble(fread("./meta/lib/comb_lib_meta.csv"))
  })
  
  observeEvent(input$append_run_meta, {
    temp <- tibble(run_id = min(setdiff(1:10000, rv$rv_run_meta$run_id)),
                   seq_type = input$user_seq_type,
                   seq_location = input$user_seq_location,
                   seq_completion_date = input$user_seq_completion_date,
                   fastq_folder_name = input$user_fastq_folder_name)
    rv$rv_run_meta <- rv$rv_run_meta %>%
      bind_rows(temp)
  })
  
  observeEvent(input$save_run_meta, {
    write.csv(rv$rv_run_meta, "./meta/run/run_meta.csv", row.names = F)
  })
  
  observeEvent(input$reload_run_meta, {
    rv$rv_run_meta <- as_tibble(fread("./meta/run/run_meta.csv"))
  })
  
  observeEvent(input$append_targ_meta, {
    temp <- tibble(targ_id = min(setdiff(1:10000, rv$rv_targ_meta$targ_id)),
                   targ_name = input$user_targ_name)
    rv$rv_targ_meta <- rv$rv_targ_meta %>%
      bind_rows(temp)
  })
  
  observeEvent(input$save_targ_meta, {
    write.csv(rv$rv_targ_meta, "./meta/samp/targ_meta.csv", row.names = F)
  })
  
  observeEvent(input$reload_targ_meta, {
    rv$rv_targ_meta <- as_tibble(fread("./meta/samp/targ_meta.csv"))
  })
  
  observeEvent(input$append_resin_meta, {
    temp <- tibble(resin_id = min(setdiff(1:10000, rv$rv_resin_meta$resin_id)),
                   resin_name = input$user_resin_name)
    rv$rv_resin_meta <- rv$rv_resin_meta %>%
      bind_rows(temp)
  })
  
  observeEvent(input$save_resin_meta, {
    write.csv(rv$rv_resin_meta, "./meta/samp/resin_meta.csv", row.names = F)
  })
  
  observeEvent(input$reload_resin_meta, {
    rv$rv_resin_meta <- as_tibble(fread("./meta/samp/resin_meta.csv"))
  })
  
  observeEvent(input$reload_samp_meta, {
    rv$rv_samp_meta <- as_tibble(fread("./meta/samp/samp_meta.csv"))
  })
  
  observeEvent(input$reload_merged_meta, {
    rv$rv_merged_meta <- rv$rv_run_meta %>%
      left_join(ct_file_meta) %>%
      left_join(rv$rv_samp_meta) %>%
      left_join(rv$rv_targ_meta) %>%
      left_join(rv$rv_resin_meta)
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
                                                fname = pull(temp_out, seq_filename)[i] %>% unbox(),
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
                                              save_often = 2 %>% unbox(),
                                              outfile = str_remove_all(input$config_upload$name, ".csv") %>% unbox(),
                                              reprocess_reads = F %>% unbox(),
                                              reprocess_counts = F %>% unbox(),
                                              reprocess_enrichments = T %>% unbox(),
                                              njobs = 4 %>% unbox()
                                            )
                                            
                                            write(toJSON(cfg), file)
                                          },
                                          contentType = "json")
}

shinyApp(ui, server)


# jitter plot (3/19/2021) ####
meta <- seq_meta %>%
  left_join(seq_smp_meta) %>%
  left_join(smp_meta) %>%
  left_join(tgt_meta)

# manually generate table of an_expts
bl_ids <- meta %>%
  filter(seq_type == "HiSeq" & scrn_id != 3) %>%
  pull(seq_smp_id)

in_tbl <- tibble(tgt_seq_smp_lst = c(list(c(603,604)), list(c(605,606)), list(c(607,608)), list(c(603,604)), list(c(605,606)), list(c(607,608))),
                 ctrl_seq_smp_lst = c(list(c(601,602)), list(c(601,602)), list(c(601,602)), list(bl_ids), list(bl_ids), list(bl_ids)),
                 lvl = 3,
                 an_type_id = c(2,2,2,5,5,5),
                 ci = 0.95) %>%
  mutate(tgt_seq_smp_lst = as.character(tgt_seq_smp_lst),
         ctrl_seq_smp_lst = as.character(ctrl_seq_smp_lst)) %>%
  rowwise() %>%
  mutate(tgt_seq_smp_lst = list(eval(parse(text = tgt_seq_smp_lst))),
         ctrl_seq_smp_lst = list(eval(parse(text = ctrl_seq_smp_lst)))) %>%
  ungroup()

# generate/assign analysis ids and update an_meta
in_tbl <- assign_an_ids(in_tbl)

# perform analyses to get lists of values for plotting
an_meta <- as_tibble(fread("./meta/an/an_meta.csv")) %>%
  rowwise() %>%
  mutate(tgt_seq_smp_lst = list(eval(parse(text = tgt_seq_smp_lst))),
         ctrl_seq_smp_lst = list(eval(parse(text = ctrl_seq_smp_lst)))) %>%
  ungroup()

prfm_an_expts(in_tbl)

# generate jitter plot
in_tbl <- tibble(y_an_id = c(2,3,4),
                 col_an_id = c(5,6,10))

vs_df <- c()
for(i in 1:nrow(in_tbl)) {
  temp <- as_tibble(fread(paste0("./result/an/", in_tbl$y_an_id[i], ".csv"))) %>%
    rename(y = value) %>%
    left_join(as_tibble(fread(paste0("./result/an/", in_tbl$col_an_id[i], ".csv")))) %>%
    rename(col = value) %>%
    mutate(an_id = in_tbl$y_an_id[i])
  vs_df <- bind_rows(vs_df, temp)
}

df <- vs_df %>%
  left_join(an_meta)

genetare_jitterplot(
  df, 
  lvl = 3,
  plot_title = "", 
  x_lab = "", 
  y_lab = "", 
  col_lab = "",
  top_n = 100000,
  width = 12, # in inches
  height = 6, # in inches
  out_loc = "./result/vs/", 
  out_name = "plt001_002"
)

# monosynthon heatmap (3/19/2021) ####
meta <- seq_meta %>%
  left_join(seq_smp_meta) %>%
  left_join(smp_meta) %>%
  left_join(tgt_meta)

# manually generate table of an_expts
in_tbl <- tibble(tgt_seq_smp_lst = c(list(c(603,604)), list(c(605,606)), list(c(607,608))),
                 ctrl_seq_smp_lst = c(list(c(601,602)), list(c(601,602)), list(c(601,602))),
                 lvl = 1,
                 an_type_id = 2,
                 ci = 0.95) %>%
  mutate(tgt_seq_smp_lst = as.character(tgt_seq_smp_lst),
         ctrl_seq_smp_lst = as.character(ctrl_seq_smp_lst)) %>%
  rowwise() %>%
  mutate(tgt_seq_smp_lst = list(eval(parse(text = tgt_seq_smp_lst))),
         ctrl_seq_smp_lst = list(eval(parse(text = ctrl_seq_smp_lst)))) %>%
  ungroup()

# generate/assign analysis ids and update an_meta
in_tbl <- assign_an_ids(in_tbl)

# perform analyses to get lists of values for plotting
an_meta <- as_tibble(fread("./meta/an/an_meta.csv")) %>%
  rowwise() %>%
  mutate(tgt_seq_smp_lst = list(eval(parse(text = tgt_seq_smp_lst))),
         ctrl_seq_smp_lst = list(eval(parse(text = ctrl_seq_smp_lst)))) %>%
  ungroup()

prfm_an_expts(in_tbl)

# generate monosynthon heatmap
vs_cy_no <- 1
an_ids <- 7:9

vs_df <- c()
for(an_id in an_ids) {
  temp <- as_tibble(fread(paste0("./result/an/", an_id, ".csv"))) %>%
    mutate(an_id = an_id)
  vs_df <- bind_rows(vs_df, temp)
}
vs_df <- vs_df %>%
  left_join(in_tbl) %>%
  left_join(an_meta)

df <- vs_df

generate_heatmap_ms(
  df, 
  vs_cy_no = 1,
  plot_title = "", 
  x_lab = "", 
  y_lab = "", 
  col_lab = "",
  width = 12, # in inches
  height = 6, # in inches
  out_loc = "./result/vs/", 
  out_name = "plt_001_001"
)



# first pass analysis (3/22/2021) ####
meta <- seq_meta %>%
  left_join(seq_smp_meta) %>%
  left_join(smp_meta) %>%
  left_join(tgt_meta)

# input from user - lists for tgts and ctrls
tgt_seq_smp_lst <- 603:608
ctrl_seq_smp_lst <- 601:602
res_name <- "res0001"

# generate replicate scatter plots
submeta <- meta %>%
  filter(seq_smp_id %in% c(tgt_seq_smp_lst, ctrl_seq_smp_lst)) %>%
  group_by(across(-c("seq_smp_id", "smp_id"))) %>%
  mutate(group_id = cur_group_id()) %>%
  ungroup()

gids <- unique(submeta$group_id)
for(gid in gids) {
  ids <- submeta %>%
    filter(group_id == gid)
  for(i in 1:(nrow(ids)-1)) {
    for(j in (i+1):nrow(ids)) {
      id1 <- ids$seq_smp_id[i]
      id2 <- ids$seq_smp_id[j]
      generate_replicate_scatterplots(
        id1,
        id2,
        plot_title = "replicate scatterplot", 
        width = 6, # in inches
        height = 6, # in inches
        out_loc = paste0("./result/", res_name, "/replacate_scatterplot/")
      )
    }
  }
}

# parse through input from user and automatically generate relevant analysis experiments
submeta <- meta %>%
  filter(seq_smp_id %in% tgt_seq_smp_lst) %>%
  group_by(across(-c("seq_smp_id", "smp_id"))) %>%
  mutate(group_id = cur_group_id()) %>%
  ungroup()

gids <- unique(submeta$group_id)
curr_tgt_seq_smp_lst <- c()
for(gid in gids) {
  tids <- submeta %>%
    filter(group_id == gid) %>%
    pull(seq_smp_id) %>%
    as.double() %>%
    list()
  curr_tgt_seq_smp_lst <- c(curr_tgt_seq_smp_lst, tids)
}

in_tbl <- tibble(tgt_seq_smp_lst = rep(curr_tgt_seq_smp_lst, 3),
                 ctrl_seq_smp_lst = list(as.double(ctrl_seq_smp_lst)),
                 lvl = rep(1:3, each = length(curr_tgt_seq_smp_lst)),
                 an_type_id = 2,
                 ci = 0.95) %>%
  mutate(tgt_seq_smp_lst = as.character(tgt_seq_smp_lst),
         ctrl_seq_smp_lst = as.character(ctrl_seq_smp_lst)) %>%
  rowwise() %>%
  mutate(tgt_seq_smp_lst = list(eval(parse(text = tgt_seq_smp_lst))),
         ctrl_seq_smp_lst = list(eval(parse(text = ctrl_seq_smp_lst)))) %>%
  ungroup()

# generate/assign analysis ids and update an_meta
in_tbl <- assign_an_ids(in_tbl)

# perform analyses to get lists of values for plotting
an_meta <- as_tibble(fread("./meta/an/an_meta.csv")) %>%
  rowwise() %>%
  mutate(tgt_seq_smp_lst = list(eval(parse(text = tgt_seq_smp_lst))),
         ctrl_seq_smp_lst = list(eval(parse(text = ctrl_seq_smp_lst)))) %>%
  ungroup()

prfm_an_expts(in_tbl)

# generate monosynthon heatmaps
an_ids <- in_tbl %>%
  filter(lvl == 1) %>%
  pull(an_id)

for(curr_cy_no in 1:3) {
  vs_df <- c()
  for(an_id in an_ids) {
    temp <- as_tibble(fread(paste0("./result/an/", sprintf("%06d", an_id), ".csv"))) %>%
      mutate(an_id = an_id)
    vs_df <- bind_rows(vs_df, temp)
  }
  vs_df <- vs_df %>%
    left_join(in_tbl) %>%
    left_join(an_meta)
  
  generate_heatmap_ms(
    df = vs_df, 
    vs_cy_no = curr_cy_no,
    plot_title = paste0("cycle ", curr_cy_no, " monosynthon heatmap"), 
    x_lab = "building block", 
    y_lab = "sample(s)", 
    col_lab = "log2(enrichment, lb)",
    width = 12, # in inches
    height = 6, # in inches
    out_loc = paste0("./result/", res_name, "/heatmap/"), 
    out_name = paste0("cycle", curr_cy_no)
  )
}

# generate jitter plots
for(syn_lvl in 2:3) {
  an_ids <- in_tbl %>%
    filter(lvl == syn_lvl) %>%
    pull(an_id)
  
  vs_df <- c()
  for(i in 1:length(an_ids)) {
    temp <- as_tibble(fread(paste0("./result/an/", sprintf("%06d", an_ids[i]), ".csv"))) %>%
      rename(y = value) %>%
      mutate(an_id = an_ids[i])
    vs_df <- bind_rows(vs_df, temp)
  }
  
  df <- vs_df %>%
    left_join(an_meta)
  
  generate_jitterplot(
    df, 
    lvl = syn_lvl,
    coloring = F,
    plot_title = paste0(syn2name(syn_lvl), " jitter plot"), 
    x_lab = "sample(s)", 
    y_lab = "enrichment, lb", 
    top_n = 100000,
    width = 12, # in inches
    height = 6, # in inches
    out_loc = paste0("result/", res_name, "/jitter_plot/"), 
    out_name = paste0(syn2name(syn_lvl), "_jitter_plot")
  )
}

# get folder structure and file column headers (3/29/2021) ####
folders <- list.dirs(getwd()) %>%
  str_replace_all(getwd(),
                  "/chembio/datasets/proj/del")
tibble(folder_name = folders) %>%
  write.csv("del_project_folder_structure_v2.csv", row.names = F)

csv <- list.files(getwd(), pattern = "csv", recursive = T)
img <- list.files(getwd(), pattern = "png", recursive = T)

csv_tbl <- c()
for(i in 1:length(csv)) {
  tcsv <- csv[i]
  splitloc <- gregexpr("\\/", tcsv)[[1]]
  fname <- substring(tcsv, splitloc[length(splitloc)]+1)
  fdir <- substring(tcsv, 1, splitloc[length(splitloc)])
  temp <- tibble(data_dir_name = fdir,
                 data_file_name = fname,
                 column_header = names(read.csv(tcsv, nrows = 1)))
  csv_tbl <- bind_rows(csv_tbl, temp)
}

img_tbl <- c()
for(i in 1:length(img)) {
  timg <- img[i]
  splitloc <- gregexpr("\\/", timg)[[1]]
  fname <- substring(timg, splitloc[length(splitloc)]+1)
  fdir <- substring(timg, 1, splitloc[length(splitloc)])
  temp <- tibble(data_dir_name = fdir,
                 data_file_name = fname,
                 column_header = "")
  img_tbl <- bind_rows(img_tbl, temp)
}

bind_rows(csv_tbl, img_tbl) %>%
  mutate(data_dir_name = paste0("/chembio/datasets/proj/del/",
                                data_dir_name)) %>%
  write.csv("del_project_columns_v2.csv", row.names = F)


