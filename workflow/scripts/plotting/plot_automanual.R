library(tidyverse)

# Read in data
scenario_path <- "saved_data/automanual_scenario.tsv.gz"
scenario_data <- suppressMessages(read_tsv(scenario_path)) %>%
  filter(data_limit_auto == 14, p_smartphone_link == 0.8)

# Make plot
# label_limit <- function(x) {
#   x1 <- ifelse(as.numeric(x) == Inf, "No",
#                paste0(x, "-day"))
#   return(paste(x1, "limit"))
# }
# label_gen <- function(x) paste0(x, "-generation delay")
g <- ggplot(scenario_data, aes(x=p_traced_manual, y=p_controlled,
                               colour = factor(p_traced_auto),
                               fill = factor(p_traced_auto))) +
  geom_ribbon(aes(ymin=p_controlled_lower, ymax=p_controlled_upper),
              alpha=0.3, colour = NA) +
  geom_line() + geom_point(size=2) +
  scale_y_continuous(name = "% of outbreaks controlled", limits = c(0,1),
                     breaks = seq(0,1,0.2), labels = function(x) round(x*100)) +
  scale_x_continuous(name = "% of unlinked contacts traced", breaks = seq(0,1,0.2),
                     labels = function(x) round(x*100)) +
  #facet_grid(p_smartphone_link~data_limit_auto) +
  scale_colour_brewer(type = "div", palette = "Dark2", labels=function(x)as.numeric(x)*100,
                      name="% of smartphone-linked\ncontacts traced") +
  scale_fill_brewer(type = "div", palette = "Dark2", labels=function(x)as.numeric(x)*100,
                      name="% of smartphone-linked\ncontacts traced") +
  theme_bw() + theme(
    legend.position = "right",
    panel.spacing.x = unit(0.3, "cm")
  )

ggsave(filename="saved_data/plot_automanual.png",
       plot = g, device="png",
       width=22, height=12, units="cm", dpi=320, limitsize=FALSE)
