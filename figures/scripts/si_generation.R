#==============================================================================
# Preamble
#==============================================================================

logfile <- file(snakemake@log[[1]], open = "wt")
sink(logfile ,type = "output")
sink(logfile, type = "message")

cat("Preparing to run script...")

# Source formatting
source("scripts/aux_format-plots.R")

# Specify parameters
simul_nrep <- snakemake@params[["simul_nrep"]]
plot_scale_cm <- snakemake@params[["panel_scale"]]
# simul_nrep <- 1e6
# plot_scale_cm <- 9
plot_scale_in <- plot_scale_cm/2.54
label_vjust <- 2

# Specify output path
output_path <- snakemake@output[[1]]
# output_path <- "output_files/dev_si_generation.png"

#==============================================================================
# Auxiliary functions
#==============================================================================

# Incubation time distribution
inc <- function(n=simul_nrep){
  rlnorm(n=n, meanlog = 1.644, sdlog = 0.363)
}

# Generation time distribution
gen <- function(a, n=simul_nrep){
  sn::rsn(n=n, omega=2, alpha=a)
} 

# Combined distribution
incgen <- function(a, n = simul_nrep){
  i <- inc(n)
  g <- gen(a,n)
  h <- pmax(i + g, 0)
  return(h)
}

# Data table
incgen_tab <- function(a, scenario, n = simul_nrep){
  tibble(n=1:n, g=incgen(a,n), scenario=scenario)
}

cat("done.\n")

#==============================================================================
# Per-scenario histograms
#==============================================================================

cat("\nPreparing generation-time histograms...")

incgen_median <- incgen_tab(0.064, "median")
incgen_optim  <- incgen_tab(0.397, "optimistic")
incgen_pessim <- incgen_tab(-0.095, "pessimistic")
incgen_all <- bind_rows(incgen_median, incgen_optim, incgen_pessim)

g_hist <- ggplot(incgen_all) + 
  geom_histogram(aes(x=g, fill=scenario), binwidth=0.5, boundary=0) + 
  facet_grid(scenario~.) + 
  scale_x_continuous(name="Generation time (days)", limits=c(0,20),
                     expand = c(0,0)) +
  scale_y_continuous(name="% of cases", labels=function(y) y/simul_nrep * 100,
                     limits=c(0,10/100*simul_nrep),
                     breaks=seq(0,10/100*simul_nrep,2/100*simul_nrep),
                     expand = c(0,0)) +
  scale_fill_brewer(palette="Dark2", name="Scenario") +
  theme_base + theme(strip.text.y = element_blank(),
                     panel.spacing.y = unit(0.6, "cm"),
                     aspect.ratio = 1/4)

cat("done.\n")

#==============================================================================
# Rate as a function of alpha
#==============================================================================

cat("\nPreparing line plots...")

alphas <- c(-1.38, -0.73, -0.325, 0, 0.325, 0.73, 1.38, 0.064, 0.397, -0.095)
p_presym <- c(0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.48, 0.38, 0.53)
a_presym <- tibble(alpha = alphas, p_presym = p_presym)
  
incgen_alphas <- lapply(alphas, function(a) incgen_tab(a, a)) %>%
  bind_rows %>% rename(alpha = scenario) %>%
  inner_join(a_presym, by="alpha")
p_zero_alpha <- incgen_alphas %>% group_by(alpha, p_presym) %>%
  summarise(`= 0` = mean(g == 0), `≤ ½` = mean(g<=0.5)) %>%
  gather(condition, rate, -alpha, -p_presym)
g_zero_alpha <- ggplot(p_zero_alpha, aes(x=alpha, y=rate, colour=condition)) +
  geom_line() + geom_point(size=2) +
  geom_vline(xintercept=c(0.064, 0.397, -0.095), linetype="dashed",
             colour=c("#1B9E77", "#D95F02", "#7570B3")) +
  scale_y_continuous(name="% of cases", labels=function(y) y*100,
                     limits=c(0,0.06), breaks=seq(0,1,0.02),
                     expand = c(0,0)) +
  scale_x_continuous(name=expression(alpha)) +
  scale_colour_brewer(palette="Set1", name=NULL,
                      labels = c("Generation time = 0", "Generation time ≤ ½")) +
  theme_base + theme(aspect.ratio = 0.8)
g_zero_presym <- ggplot(p_zero_alpha, aes(x=p_presym, y=rate, colour=condition)) +
  geom_line() + geom_point(size=2) +
  geom_vline(xintercept=c(0.48, 0.38, 0.53), linetype="dashed",
             colour=c("#1B9E77", "#D95F02", "#7570B3")) +
  scale_y_continuous(name="% of cases", labels=function(y) y*100,
                     limits=c(0,0.06), breaks=seq(0,1,0.02),
                     expand = c(0,0)) +
  scale_x_continuous(name="% pre-symptomatic\ntransmission", 
                     labels=function(x) x*100) +
  scale_colour_brewer(palette="Set1", name=NULL,
                      labels = c("Generation time = 0", "Generation time ≤ ½")) +
  # guides(colour = guide_legend(nrow = 2)) +
  theme_base + theme(aspect.ratio = 0.8)

#==============================================================================
# Plot of alpha vs presymptomatic transmission
#==============================================================================

g_alpha <- p_zero_alpha %>%
  ggplot(aes(x=alpha, y=p_presym)) +
  geom_line() + geom_point(size=2) +
  scale_x_continuous(name=expression(alpha)) +
  scale_y_continuous(name="% pre-symptomatic\ntransmission", 
                     labels=function(x) x*100) +
  #coord_fixed(ratio = 3) +
  theme_base + theme(aspect.ratio = 0.95)

cat("done.\n")

#==============================================================================
# Grid plot
#==============================================================================

cat("\nAssembling output grid...")

# First row: alpha plot and histograms
legend_hist <- get_legend(g_hist)
grid_alpha_hist <- plot_grid(g_alpha, g_hist + theme(legend.position = "none"), 
                             ncol=2, labels = c("a", "b"), align = "hv", 
                             vjust=label_vjust,
                       axis = "lrtb", label_size = fontsize_base * fontscale_label,
                       label_fontfamily = titlefont, label_colour = "black")

# Second row: % of zero cases
legend_zero <- get_legend(g_zero_alpha)
grid_zero <- plot_grid(g_zero_alpha + theme(legend.position = "none"), 
                       g_zero_presym + theme(legend.position = "none"), ncol=2,
                       labels = c("c","d"), align = "hv", vjust=label_vjust,
                       axis = "lrtb", label_size = fontsize_base * fontscale_label,
                       label_fontfamily = titlefont, label_colour = "black")

# Combine and write
grid_out <- plot_grid(grid_alpha_hist, legend_hist, grid_zero, legend_zero, 
                      nrow = 4,
                      rel_heights = c(1, 0.08, 1, 0.08))

cat("done.\n")

#==============================================================================
# Save output
#==============================================================================

cat("\nSaving output...")

cowplot::save_plot(filename=output_path, plot=grid_out,
                   ncol = 3, nrow = 3, base_height = plot_scale_in,
                   base_asp = 1)
cat("done.\n")
