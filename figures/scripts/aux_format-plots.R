#------------------------------------------------------------------------------
# Packages
#------------------------------------------------------------------------------

load_packages <- function(package_list){
  for (p in package_list){
    suppressMessages(suppressWarnings(library(p, character.only = TRUE)))
  }
}

packages_default <- c("tidyverse", "cowplot", "metR", "reshape2", "scales",
                      "viridis", "grid")

load_packages(packages_default)

#------------------------------------------------------------------------------
# Fonts
#------------------------------------------------------------------------------

font <- "sans" # Main figure font
titlefont <- font # Font for axis titles etc. (if different)
fontsize_base <- 12 # Basic figure font size
fontscale_title <- 1.5 # Default axis-title scale relative to regular font
fontscale_main <- 1.5 # Default plot-title scale
fontscale_mid <- 1.25
fontscale_label <- 2 # Default subfigure-label scale (A, B, etc)
fontscale_legend <- 1

#------------------------------------------------------------------------------
# Plotting theme information
#------------------------------------------------------------------------------

theme_base <-   theme_bw() + theme(
  legend.position = "bottom",
  axis.text = element_text(size = fontsize_base, family = font, 
                           colour = "black"),
  axis.title.y = element_text(size = fontsize_base * fontscale_title,
                              family = titlefont, colour = "black",
                              margin = margin(r=3,unit="mm")),
  axis.title.x = element_text(size = fontsize_base * fontscale_title,
                              family = titlefont, colour = "black",
                              margin = margin(t=5,unit="mm")),
  legend.text = element_text(size = fontsize_base * fontscale_legend, 
                             family = font, colour = "black"),
  legend.title = element_text(size = fontsize_base * fontscale_legend,
                              family = font, colour = "black",
                              face = "bold", vjust=0.5, 
                              margin=margin(r=3, unit="mm")),
  plot.title = element_text(size = fontsize_base * fontscale_main,
                            family = titlefont, colour="black",
                            hjust = 0.5, vjust = 0.5,
                            margin = margin(b=4, unit="mm"),
                            face = "bold"),
  plot.margin = margin(t=1.5, l=0.5, r=0.5, b = 0.5, unit="cm"),
  strip.background = element_blank(),
  strip.text.x = element_text(size = fontsize_base * fontscale_title,
                              family = titlefont, colour = "black",
                              margin = margin(t = 3, b = 3, unit="mm")),
  strip.text.y = element_text(size = fontsize_base * fontscale_title,
                              family = titlefont, colour = "black",
                              margin = margin(l=3, r=3, unit="mm")),
  legend.justification = "center",
  panel.border = element_blank(),
  axis.line = element_line(size=0.5, colour="black")
)

#------------------------------------------------------------------------------
# Labelling functions
#------------------------------------------------------------------------------

label_backtrace <- function(x){
  ifelse(x==0, "Forward tracing only",
         "Bidirectional tracing")
  # "Forward and reverse\ntracing")
}

label_backtrace_brief <- function(x){
  ifelse(x==0, "Forward-only",
         "Bidirectional")
}

label_backtrace_gen <- function(x){
  ifelse(x==0, "Forward-only",
         ifelse(x == Inf, "Bidirectional (unlimited generations)",
                ifelse(x == 1, "Bidirectional (1 generation)",
                       paste0("Bidirectional (", x, "generations"))))
}

label_limits <- function(x){
  ifelse(as.numeric(x) == Inf, "No trace limit", paste0(x, "-day trace limit"))
}

label_limits_manual <- function(x){
  ifelse(as.numeric(x) == Inf, "No manual limit", paste0(x, "-day manual limit"))
}

label_windows <- function(x){
  ifelse(as.numeric(x) == Inf, "No trace limit", paste0(x, "-day tracing window"))
}

label_windows_manual <- function(x){
  ifelse(as.numeric(x) == Inf, "No manual limit", paste0(x, "-day manual window"))
}



label_vars <- function(x){
  ifelse(x=="p_data_sharing_auto", "% of cases sharing data",
         "% of cases with\ntrace-enabled smartphones")
}

label_pc <- function(x) round(as.numeric(x)*100)

label_test <- function(x){
  ifelse(!x, "Test not required", "Test required")
}

label_ident <- function(x) x

#------------------------------------------------------------------------------
# Plotting functions
#------------------------------------------------------------------------------

# Plot types

format_ctrl <- function(g, err_width=0.05, pt_size = 2, coord_ratio = 1,
                        ylab = "% of outbreaks controlled",
                        ylimits = c(0,1), ybreaks = seq(0,1,0.2)){
  #' Add general formatting to control plot
  g <- g + geom_errorbar(width=err_width) + geom_line() + geom_point(size=pt_size) +
    scale_y_continuous(name = ylab, limits = ylimits,
                       breaks = ybreaks, labels = label_pc)
  if (!is.null(coord_ratio)) g <- g + coord_fixed(ratio = coord_ratio)
  return(g)
}

format_reff <- function(g, r0 = NULL, pt_size = 2, coord_ratio = 1,
                        r0_col = "red", r0_lab_x = 1,
                        r0_lab_yscale = 0.97, r0_lab_size = 4.6,
                        ylab = "Mean "){
  #' Add general formatting to R_eff plot
  g <- g + geom_line() + geom_point(size=pt_size) +
    geom_hline(yintercept=1, linetype="dotted", size=1)
  if (!is.null(r0)){
    g <- g + geom_hline(yintercept = r0, linetype = "dotted", size=1,
               colour = r0_col) +
      scale_y_continuous(name = expression(paste("Mean ",italic(R)[eff])),
                         breaks = seq(0,100,0.5), limits = c(0, r0)) +
      annotate("text", x=r0_lab_x, y=r0*r0_lab_yscale, hjust=1, vjust = 1,
               colour = r0_col, size=r0_lab_size, label=expression(italic(R)[0]))
    if (!is.null(coord_ratio)) g <- g + coord_fixed(ratio = coord_ratio/r0)
  } else {
    g <- g + 
      scale_y_continuous(name = expression(paste("Mean ",italic(R)[eff])),
                         breaks = seq(0,100,0.5), limits = c(0, NA))
    if (!is.null(coord_ratio)) g <- g + coord_fixed(ratio = coord_ratio)
  }
  return(g)
}

# X-axis schemes

x_ptrace <- function(g){
  g + scale_x_continuous(name="Probability of trace\nsuccess (%)",
                         limits=c(0,1.03), breaks=seq(0,1,0.2), labels=label_pc)
}

x_ascertainment <- function(g){
  g + scale_x_continuous(name="% symptomatic cases\nidentified",
                         limits=c(0,1.03), breaks=seq(0,1,0.2), labels=label_pc)
}

x_sensitivity <- function(g){
  g + scale_x_continuous(name="Test sensitivity (%)",
                         limits=c(0,1.03), breaks=seq(0,1,0.2), labels=label_pc)
}

x_r0 <- function(g){
  g + scale_x_continuous(name=expression(italic(R)[0]))
}

# Colour schemes

colour_backtrace <- function(g, name=NULL, ncol = 1, order = 1){
  #' Colour plot geoms according to backtrace distance
  g + scale_colour_brewer(palette = "Dark2", name = name,
                          labels = label_backtrace) +
    guides(colour=guide_legend(ncol = ncol, order = order))
}

colour_ttype <- function(g, name="Trace type", ncol = 2, order = 1,
                         label = label_ident){
  #' Colour plot geoms according to backtrace distance
  g + scale_colour_brewer(palette = "Set1", name = name,
                          labels = label) +
    guides(colour=guide_legend(ncol = ncol, order = order))
}

fill_backtrace <- function(g, name=NULL, ncol = 1, order = 1){
  #' Fill plot geoms according to backtrace distance
  g + scale_fill_brewer(palette = "Dark2", name = name,
                          labels = label_backtrace) +
    guides(fill=guide_legend(ncol = ncol, order = order))
}


# Line schemes

linetype_window <- function(g, name=NULL, ncol = 1, order = 2,
                            label=label_limits){
  #' Colour plot geoms according to backtrace distance
  g + scale_linetype_discrete(name = name, labels = label) +
    scale_shape_discrete(name = name, labels = label) +
    guides(linetype = guide_legend(ncol = ncol, order = order),
           shape = guide_legend(ncol = ncol, order = order))
}

linetype_test <- function(g, name=NULL, ncol = 1, order = 2,
                            label=label_test){
  #' Colour plot geoms according to testing strategy
  g + scale_linetype_discrete(name = name, labels = label) +
    scale_shape_discrete(name = name, labels = label) +
    guides(linetype = guide_legend(ncol = ncol, order = order),
           shape = guide_legend(ncol = ncol, order = order))
}

linetype_backtrace <- function(g, name=NULL, ncol = 1, order = 2,
                            label=label_backtrace_brief){
  #' Colour plot geoms according to backtrace distance
  g + scale_linetype_discrete(name = name, labels = label) +
    scale_shape_discrete(name = name, labels = label) +
    guides(linetype = guide_legend(ncol = ncol, order = order),
           shape = guide_legend(ncol = ncol, order = order))
}


# Themes

theme_internal <- function(g, pos = c(0.02, 1), vspace=0,
                           yjust = "top", xjust="left"){
  g + theme_base + theme(
      legend.justification = c(xjust,yjust),
      legend.background = element_rect(fill=alpha("white", 0.5)),
      legend.text = element_text(margin = margin(t=1.1, b=1.1, unit="mm")),
      legend.position = pos,
      legend.spacing.y = unit(vspace, "mm"),
    )
}

theme_external <- function(g, pos="bottom"){
  g + theme_base + theme(legend.position=pos)
}

strip_upper <- function(g, join_scale = 20){
  g + theme(axis.text.x = element_blank(),
            axis.title.x = element_blank(),
            plot.margin = margin(b = -join_scale),
            legend.position = "none")
}

strip_lower <- function(g, join_scale = 20){
  g + theme(strip.text.x = element_blank(),
            strip.background.x = element_blank(),
            plot.margin = margin(t = -join_scale))
}

#------------------------------------------------------------------------------
# I/O
#------------------------------------------------------------------------------

save_plot <- function(path, plot, plot_scale,
                      plot_ratio, device=c("png")){
  for (d in device){
    width <- plot_scale * plot_ratio
    ggsave(filename=path, plot=plot, device=d, units="cm",
           dpi=320, limitsize=FALSE, width=width, height=plot_scale)
  }
}

#------------------------------------------------------------------------------
# Palettes for contour plotting
#------------------------------------------------------------------------------

# Basic palette settings
col_bidir <- "#d95f02"
col_fwd   <- "#1b9e77"
step_dense <- 0.05
step_sparse <- 0.1

make_contour_palette <- function(colour_high, step_size, data = NULL,
                                 contour_var = NULL, bt_dist = NULL, 
                                 palette_max = NULL, palette_min = NULL,
                                 reverse_scale = FALSE){
  # Calibrate range of contour palette
  if (is.null(palette_max) | is.null(palette_min)){
    data_range <- data %>% filter(backtrace_distance == bt_dist) %>%
      pull(contour_var) %>% unique
  }
  if (is.null(palette_max)){
    palette_max <- data_range %>% max %>%
      (function(x) ceiling(x/step_size)*step_size)
  }
  if (is.null(palette_min)){
    palette_min <- data_range %>% min %>%
      (function(x) floor(x/step_size)*step_size)
  }
  palette_val <- seq(palette_min, palette_max, step_size)
  # Calculate hue values
  palette_fn  <- scales::seq_gradient_pal(low="white", high=colour_high, "Lab")
  palette_hue <- palette_fn((seq(0,1,length.out=length(palette_val))))
  if (reverse_scale) palette_hue <- rev(palette_hue)
  return(tibble(level = palette_val, colour = palette_hue))
}

make_contour_palette_ctrl <- function(colour_high, step_size){
  make_contour_palette(colour_high, step_size, palette_min = 0, palette_max = 1) %>%
    mutate(level = round(level * 100))
}

make_contour_palette_reff <- function(colour_high, step_size, data,
                                      bt_dist){
  make_contour_palette(colour_high, step_size, data = data, bt_dist = bt_dist,
                       reverse_scale = TRUE, contour_var = "effective_r0_mean")
}

#------------------------------------------------------------------------------
# Data processing for contour plotting
#------------------------------------------------------------------------------

ttype_levels <- c("Manual only", "Manual + digital",
                  "Digital only", "No tracing")

label_data <- function(data, p_traced = 0.9){
  data %>% filter(p_traced_auto %in% c(-1, 0, p_traced),
                  p_traced_manual %in% c(-1, 0, p_traced)) %>%
    mutate(bt_type = ifelse(backtrace_distance == 0, "Forward-only",
                            "Bidirectional"),
           trace_window_manual = paste0(contact_limit_manual, "-day manual window"),
           trace_window_auto = paste0(contact_limit_manual, "-day digital window"),
           trace_type = ifelse(p_traced_auto == 0,
                               ifelse(p_traced_manual == 0, "No tracing", "Manual only"),
                               ifelse(p_traced_manual == 0, "Digital only", "Manual + digital"))) %>%
    mutate(trace_type = factor(trace_type, levels = ttype_levels))
}

index_data <- function(data, key_x = "p_asymptomatic", key_y = "rel_r0_asymptomatic"){
  # Compute index levels
  levels_x <- sort(unique(data[[key_x]]))
  levels_y <- sort(unique(data[[key_y]]))
  # Assign indices
  data_indexed <- data %>%
    mutate(index_x = match(!!sym(key_x), levels_x),
           index_y = match(!!sym(key_y), levels_y))
  return(data_indexed)
}

smooth_data_single <- function(data){
  indices_x <- sort(unique(data$index_x))
  indices_y <- sort(unique(data$index_y))
  index_tab <- expand_grid(index_x = indices_x, index_y = indices_y)
  smoothed_control <- sapply(1:nrow(index_tab), function(n)
    data %>% filter(abs(index_x - index_tab$index_x[n]) <= 1,
                    abs(index_y - index_tab$index_y[n]) <= 1) %>%
      pull(p_controlled) %>% mean)
  smoothed_r_eff <- sapply(1:nrow(index_tab), function(n)
    data %>% filter(abs(index_x - index_tab$index_x[n]) <= 1,
                    abs(index_y - index_tab$index_y[n]) <= 1) %>%
      pull(effective_r0_mean) %>% mean)
  index_tab_smoothed <- mutate(index_tab,
                               p_controlled_smoothed = smoothed_control,
                               r_eff_smoothed = smoothed_r_eff)
  data_smoothed <- inner_join(data, index_tab_smoothed,
                              by=c("index_x", "index_y"))
  return(data_smoothed)
}

smooth_data <- function(data, group_vars = c("trace_type", "row_label")){
  data_split <- data %>% group_by_at(group_vars) %>% group_split()
  data_smoothed_split <- lapply(data_split, smooth_data_single)
  return(bind_rows(data_smoothed_split) %>% ungroup %>%
           mutate(pc_controlled_smoothed = p_controlled_smoothed*100))
}

prep_contour_data <- function(data, p_traced, key_x, key_y,
                              group_vars = c("trace_type")){
  data %>% label_data(p_traced) %>% index_data(key_x, key_y) %>%
    smooth_data(group_vars)
}

#------------------------------------------------------------------------------
# Contour plots of processed data
#------------------------------------------------------------------------------

contour_plot_base <- function(data, x_key, y_key, z_key, contour_breaks,
                              palette_df,
                              label_map = aes(label = ..level..),
                              coord_ratio = 1, lwd = 0.8){
  g_contour <- data %>%
    ggplot(aes_string(x=x_key, y=y_key, z=z_key)) +
    geom_contour_filled(colour="black", breaks = contour_breaks) +
    geom_text_contour(colour="black", stroke = 0.2, breaks = contour_breaks,
                      check_overlap = TRUE, mapping = label_map) +
    coord_fixed(ratio = coord_ratio) +
    theme_base + theme(
      legend.position = "none",
    )
  min_level <- data %>% pull(z_key) %>% min
  palette <- palette_df %>% mutate(level = level - min_level) %>%
    group_by(level > 0) %>% filter(level > 0 | row_number() == n()) %>% pull(colour)
  return(g_contour + scale_fill_manual(values = palette))
}

#------------------------------------------------------------------------------
# Absolute values
#------------------------------------------------------------------------------

contour_plot_ctrl_abs <- function(data, x_key, y_key, palette,
                                  breaks=seq(0,100,10),
                                  coord_ratio=1, lwd=0.8){
  contour_plot_base(data, x_key, y_key, "pc_controlled_smoothed",
                    breaks, palette,
                    coord_ratio = coord_ratio, lwd = lwd)
}

contour_plot_reff_abs <- function(data, x_key, y_key, palette,
                                  breaks=seq(0,10,0.1),
                                  coord_ratio=1, lwd=0.8){
  contour_plot_base(data, x_key, y_key, "r_eff_smoothed",
                    breaks, palette,
                    coord_ratio = coord_ratio, lwd = lwd)
}

#==============================================================================
# Make contour plots
#==============================================================================

#------------------------------------------------------------------------------
# Base functions
#------------------------------------------------------------------------------

contour_plot_reff <- function(data, inc_manual, bt_dist,
                              palette_tib = bidir_palette_dense,
                              contour_skip = 0, step_size = 0.05){
  # Prepare data
  data_filtered <- data %>% filter((p_traced_manual != 0) == inc_manual,
                                   backtrace_distance == bt_dist)
  data_indexed <- data_filtered %>%
    group_by(p_smartphone_overall) %>% arrange(p_smartphone_overall) %>% 
    mutate(index_smartphone = group_indices()) %>%
    group_by(p_data_sharing_auto) %>% arrange(p_data_sharing_auto) %>% 
    mutate(index_sharing = group_indices())
  # Smooth data
  smoothed_reff <- sapply(1:max(data_indexed$index_smartphone), function(p)
    sapply(1:max(data_indexed$index_sharing), function(s)
      data_indexed %>% filter(abs(index_smartphone-p) <= 1,
                              abs(index_sharing-s) <= 1) %>% 
        pull(effective_r0_mean) %>% mean)) %>%
    melt() %>% as_tibble %>%
    rename(index_smartphone = Var2, index_sharing = Var1, r_eff_smoothed = value)
  data_smoothed <- inner_join(data_indexed, smoothed_reff, 
                              by=c("index_sharing", "index_smartphone"))
  # Prepare palette
  min_level <- data_smoothed %>% pull(r_eff_smoothed) %>% 
    (function(x) floor(x*10)/10) %>% min
  palette <- palette_tib %>% filter(level >= min_level) %>% pull(colour)
  # Make plot
  g <- ggplot(data_smoothed, aes(x=p_smartphone_overall, y=p_data_sharing_auto,
                                 z = r_eff_smoothed)) +
    geom_contour_filled(colour="black", breaks=seq(0,10,step_size)) +
    geom_text_contour(colour="black", stroke=0.2, breaks=seq(0,10,step_size),
                      skip=contour_skip) +
    scale_x_continuous(name = "% of cases with\nchirping smartphones",
                       breaks = seq(0,1,0.2), labels = function(x) round(x*100)) +
    scale_y_continuous(name = "% of cases sharing data", breaks = seq(0,1,0.2),
                       labels = function(x) round(x*100)) +
    scale_fill_manual(values = palette) +
    coord_fixed() +
    theme_base + theme(legend.position = "none")
  return(g)
}