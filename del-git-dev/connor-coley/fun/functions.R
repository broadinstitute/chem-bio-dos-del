config <- read_json("config_file.json")
date_to_quarter <- function(x) {
  temp <- str_split(x, "/") %>% unlist()
  yr <- substr(temp[3], nchar(temp[3])-1, nchar(temp[3]))
  temp_map <- structure(rep(1:4, each = 3), names = 1:12)
  quart <- temp_map[as.character(as.numeric(temp[1]))]
  return(paste0("'", yr, "/Q", quart))
} # converts a date in form MM/DD/YYYY to 'YYQ[1-4]
floor_dec <- function(x, level=1) {
  round(x - 5*10^(-level-1), level)
} # round down at the decimal level
ceiling_dec <- function(x, level=1) {
  round(x + 5*10^(-level-1), level)
} # round up at the decimal level
get_range <- function(x) {
  if(grepl("[[:digit:]]-[[:digit:]]",x)) {
    temp <- x %>% str_split("-") %>% unlist()
    return(seq(temp[1], temp[2]))
  }
  else {
    return(str_remove_all(x, "-"))
  }
} # converts a text range (e.g. 2-5) to a list
z_func_prob <- function(k1, k2, n1, n2, r) {
  if(length(r) > 1) {
    res <- sapply(r, function(r) { dnorm(2*(sqrt(k1+3/8)-sqrt(n1/n2*r)*sqrt(k2+3/8))/sqrt(1+n1/n2*r)) })
    return(res)
  }
} # calculates relative probability of a given enrichment ratio
get_er <- function(k1, k2, n1, n2, z = 0) {
  k1 <- as.numeric(k1)
  k2 <- as.numeric(k2)
  n1 <- as.numeric(n1)
  n2 <- as.numeric(n2)
  if(z == 0) {
    res <- n2/n1*(k1 + 3/8)/(k2 + 3/8)
  }
  else if(z > 0) {
    res <- (-sqrt(4*n1^2*n2^2*(-64*k1^2+32*k1*z^2-48*k1-4*z^4+12*z^2-9)*(64*k2^2-32*k2*z^2+48*k2+4*z^4-12*z^2+9)+
                    n1^2*n2^2*(-128*k1*k2-32*k1*z^2-48*k1-32*k2*z^2-48*k2+8*z^4-24*z^2-18)^2)-
              n1*n2*(-128*k1*k2-32*k1*z^2-48*k1-32*k2*z^2-48*k2+8*z^4-24*z^2-18))/(2*n1^2*(64*k2^2-32*k2*z^2+48*k2+4*z^4-12*z^2+9))
  }
  else {
    res <- (sqrt(4*n1^2*n2^2*(-64*k1^2+32*k1*z^2-48*k1-4*z^4+12*z^2-9)*(64*k2^2-32*k2*z^2+48*k2+4*z^4-12*z^2+9)+
                   n1^2*n2^2*(-128*k1*k2-32*k1*z^2-48*k1-32*k2*z^2-48*k2+8*z^4-24*z^2-18)^2)-
              n1*n2*(-128*k1*k2-32*k1*z^2-48*k1-32*k2*z^2-48*k2+8*z^4-24*z^2-18))/(2*n1^2*(64*k2^2-32*k2*z^2+48*k2+4*z^4-12*z^2+9))
  }
  res[2*(sqrt(k1+3/8)) <= z] <- 0
  res[k1 == 0 & k2 == 0] <- NA
  return(res)
} # calculates enrichment of (1) over (2)
assign_an_ids <- function(in_tbl) {
  an_meta <- as_tibble(fread(file.path(config$project_dir, config$analysis_metadata)))
  
  if(nrow(an_meta) == 0) {
    in_tbl <- in_tbl %>%
      mutate(an_id = NA)
  } else {
    in_tbl <- in_tbl %>%
      left_join(an_meta)
  }
  
  # Checks to see if any of the analysis experiments haven't been defined before
  new_ids <- which(is.na(in_tbl$an_id))
  
  if(!is.numeric(in_tbl$an_id)) {
    in_tbl$an_id <- as.numeric(in_tbl$an_id)
  }
  if(length(new_ids) > 0) {
    in_tbl[is.na(in_tbl$an_id),"an_id"] <- seq(max(an_meta$an_id, 0) + 1, length(new_ids) + max(an_meta$an_id, 0))
  }
  
  if(nrow(an_meta) == 0) {
    in_tbl %>%
      select(-targ_description, -ctrl_description) %>%
      write.csv(file.path(config$project_dir, config$analysis_metadata), row.names = F)
  } else {
    an_meta %>%
      bind_rows(in_tbl[new_ids,] %>% select(-targ_description, -ctrl_description)) %>%
      write.csv(file.path(config$project_dir, config$analysis_metadata), row.names = F)
  }
  return(in_tbl)
} # assigns an_ids and updates an_meta; expects rows for targ_ct_lst, ctrl_ct_lst, agg_id, an_type_id, ci
prfm_an_expts <- function(in_tbl) {
  an_complete <- list.files(file.path(config$project_dir, config$result_analysis_location))
  if(nrow(in_tbl > 0)) {
    for(i in 1:nrow(in_tbl)) {
      targ_ids <- unlist(str_split(in_tbl$targ_ct_lst[i], ", "))
      ctrl_ids <- unlist(str_split(in_tbl$ctrl_ct_lst[i], ", "))
      agg_ids <- as.numeric(unlist(str_split(in_tbl$agg_id[i], ", ")))
      an_type_ids <- as.numeric(unlist(str_split(in_tbl$an_type_id[i], ", ")))
      an_id <- in_tbl$an_id[i]
      filename_parts <- str_split(ctrl_ids[1], "_") %>% unlist() %>% str_remove(".csv")
      lib_id <- filename_parts[grep("lib", filename_parts)] %>% str_remove("lib") %>% as.numeric()
      ci <- in_tbl$ci[i]
      an_request <- paste0("an", sprintf("%06d", an_id), "_lib", sprintf("%03d", lib_id), 
                           "_agg", sprintf("%02d", agg_ids), "_type",  sprintf("%02d", an_type_ids), ".csv")
      # Only performs if any of the requests analyses are new
      if(length(setdiff(an_request, an_complete)) > 0) {
        # Pull counts for the traget. If more than one counts file is assigned to the sample, then the counts are summed.
        if(length(targ_ids) == 1) {
          targ_cts <- as_tibble(fread(file.path(config$project_dir, config$data_count_location, targ_ids))) %>%
            pull(value)
        } else {
          targ_cts <- 0
          for(targ_id in targ_ids) {
            targ_cts <- as_tibble(fread(file.path(config$project_dir, config$data_count_location, targ_id))) %>%
              pull(value) + targ_cts
          }
        }
        
        # Pull counts for the control. If more than one counts file is assigned to the sample, then the counts are summed.
        if(length(ctrl_ids) == 1) {
          ctrl_cts <- as_tibble(fread(file.path(config$project_dir, config$data_count_location, ctrl_ids))) %>%
            pull(value)
        } else {
          ctrl_cts <- 0
          for(ctrl_id in ctrl_ids) {
            ctrl_cts <- as_tibble(fread(file.path(config$project_dir, config$data_count_location, ctrl_id))) %>%
              pull(value) + ctrl_cts
          }
        }
        
        # Combined counts with the associated metadata for the counts (i.e., the cycle IDs)
        temp_cts <- as_tibble(fread(file.path(config$project_dir, config$library_metadata_column_location,
                                              sprintf("lib%03d_agg%02d_meta_cols.csv", lib_id, max(agg_meta$agg_id))))) %>%
          mutate(targ_ct = targ_cts,
                 ctrl_ct = ctrl_cts)
        
        for(j in agg_ids) {
          # Get the cycles defined by the aggregation level
          cycs <- paste0("cycle", unlist(str_split(agg_meta %>% filter(agg_id == j) %>% pull(cyc), ", ")))
          
          # Aggregate the counts based on the relevant cycles
          temp_agg <- temp_cts %>%
            group_by_at(c("lib_id", cycs)) %>%
            summarise(targ_ct = sum(targ_ct),
                      ctrl_ct = sum(ctrl_ct),
                      .groups = "keep")
          
          # Calculate the enrichment ratios
          for(k in an_type_ids) {
            if(an_type_meta$an_type_modifier[an_type_meta$an_type_id == k] == "e") {
              an_out <- temp_agg %>%
                ungroup() %>%
                mutate(value = get_er(targ_ct, ctrl_ct, sum(targ_ct), sum(ctrl_ct))) %>%
                select(value)
            }
            if(an_type_meta$an_type_modifier[an_type_meta$an_type_id == k] == "lb") {
              an_out <- temp_agg %>%
                ungroup() %>%
                mutate(value = get_er(targ_ct, ctrl_ct, sum(targ_ct), sum(ctrl_ct), qnorm((1+ci)/2))) %>%
                select(value)
            }
            if(an_type_meta$an_type_modifier[an_type_meta$an_type_id == k] == "ub") {
              an_out <- temp_agg %>%
                ungroup() %>%
                mutate(value = get_er(targ_ct, ctrl_ct, sum(targ_ct), sum(ctrl_ct), qnorm((1-ci)/2))) %>%
                select(value)
              
            }
            an_out %>%
              write.csv(file.path(config$project_dir, config$result_analysis_location,
                                  sprintf("an%06d_lib%03d_agg%02d_type%02d.csv", an_id, lib_id, j, k)), row.names = F)
          }
        }
      }
    }
  }
} # outputs a list of dataframes with requests enrichment calculations; expects rows for targ_ct_lst, ctrl_ct_lst, agg_id, an_type_id, ci

generate_replicate_scatterplots <- function(
  id1,
  id2,
  xlab = NA,
  ylab = NA,
  plot_title = "replicate scatterplot", 
  width = 6, # in inches
  height = 6, # in inches
  out_loc = NA
) {
  stopifnot(!(nrow(df) == 0 | is.na(out_loc)))
  dir.create(file.path(out_loc, "count_replicate_scatterplot"), recursive = T, showWarnings = F)
  ct1 <- as_tibble(fread(file.path(config$project_dir, config$data_count_location, id1))) %>%
    pull(value)
  ct2 <- as_tibble(fread(file.path(config$project_dir, config$data_count_location, id2))) %>%
    pull(value)
  if(is.na(xlab)) {
    xlab <- id1
  }
  if(is.na(ylab)) {
    xlab <- id2
  }
  tibble(x = ct1,
         y = ct2) %>%
    distinct() %>%
    ggplot(aes(x = x,
               y = y)) +
    geom_point(shape = 16,
               stroke = 0) +
    theme_bw() +
    ggtitle(plot_title) +
    xlab(paste0(str_remove_all(xlab, ".csv"))) +
    ylab(paste0(str_remove_all(ylab, ".csv"))) +
    theme(axis.text = element_text(color = "black", size = 10),
          axis.text.x = element_text(color = "black", size = 10),
          strip.text = element_text(size = 10),
          axis.title = element_text(face = "plain", size = 10),
          plot.title = element_text(face = "plain", size = 12, hjust = 0.5),
          legend.title = element_text(face = "bold", size = 12),
          legend.text = element_text(size = 10),
          legend.key = element_blank(),
          panel.border = element_rect(color = "black", size = .5),
          axis.line = element_blank(),
          axis.ticks = element_line(color = "black", size = .5),
          axis.ticks.length = unit(1.5, "mm"),
          panel.grid.minor = element_blank(),
          plot.margin=unit(c(10,5,5,5),"mm"))
  ggsave(paste0(out_loc, "/count_replicate_scatterplot/", str_remove_all(id1, ".csv"), "_", 
                str_remove_all(id2, ".csv"), ".png"), width=width, height=height, dpi = 300)
}