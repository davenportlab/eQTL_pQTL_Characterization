#----------------------------------------------------------
# GGPlot2 Theme
# Created: 28 October 2021
#----------------------------------------------------------

ggplot_theme <- theme_bw(base_size=18) +
  theme(
    panel.border=element_blank(),
    panel.grid=element_blank(),
    axis.line.x.bottom=element_line(color="black", size=0.25),
    axis.line.y.left=element_line(color="black", size=0.25),
    legend.position="bottom",
    strip.background=element_rect(fill="#EEEEEE", color="white", size=0.25)
  )
