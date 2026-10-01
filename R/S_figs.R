library(ggplot2)
library(ggtext)

# Figure 1 ----------------------------------------------------------------
fig1 <- ggplot(
  lat_summary,
  aes(
    x = horizon,
    y = mean_lat,
    color = species,
    linetype = futureScenario,
    shape = shape,
    group = interaction(species, futureScenario, year)
  )
) +
  scale_shape_manual(
    values = c("Present" = 15, "2050" = 16, "2090" = 17),
    breaks = c("2050", "2090")   # só anos na legenda
  )+
  geom_line(linewidth = 0.5) +
  geom_point(size = 1.5) +
  geom_ribbon(
    aes(ymin = lwr, ymax = upr, fill = species),
    alpha = .15,
    color = NA
  ) +
  facet_wrap(~ spGroup) +
  theme_bw(base_size = 8) +
  labs(
    x = NULL,
    y = "Latitude (mean ± CI)",
    color = "Species",
    linetype = "Scenario",
    shape = "Year",
  ) +
  guides(fill = "none",
         linetype = guide_legend(order = 1),
         shape    = guide_legend(order = 2),
         color    = guide_legend(order = 3))+
  scale_color_discrete(
    labels = function(x) paste0("<i>", x, "</i>")
  ) +
  theme(
    legend.position='none',
    legend.text = element_markdown(size = 7),
    legend.title = element_text(size = 8),
    strip.text = element_text(size = 8),    # títulos dos facets
    axis.text = element_text(size = 7),
    axis.title = element_text(size = 8)
  )

fig1

ggsave("Fig 1.png",
       device = png,
       plot = fig1,
       path = here::here("outputs", "figures"),
       width = 168,
       height = 80,
       units = "mm",
       dpi = 300,
)

# Figure 2 - area ---------------------------------------------------------

fig2 <- ggplot(
  areas,
  aes(
    x = horizon,
    y = log(area_km2),
    linetype = futureScenario,
    color = species,
    shape = shape,
    group = interaction(species, futureScenario, year)
  )) +
  scale_shape_manual(
    values = c("Present" = 15, "2050" = 16, "2090" = 17),
    breaks = c("2050", "2090")   # só anos na legenda
  )+
  geom_line(linewidth=0.5) +
  geom_point(size = 1.5) +
  theme_bw() +
  labs(
    x = NULL,
    y = "Log Predicted Area (km²)",
    linetype = "Scenario",
    color = "Species",
    shape = "Year"
  ) +
  facet_wrap(~spGroup)+
  guides(fill = "none",
         linetype = guide_legend(order = 1),
         shape = guide_legend(order = 2),
         color    = guide_legend(order = 3))+
  scale_color_discrete(
    labels = function(x) paste0("<i>", x, "</i>")
  ) +
theme(
  legend.position='none',
  legend.text = element_markdown(size = 7),
  legend.title = element_text(size = 8),
  strip.text = element_text(size = 8),    # títulos dos facets
  axis.text = element_text(size = 7),
  axis.title = element_text(size = 8)
)

fig2

ggsave("Fig 2.png",
       device = png,
       plot = fig2,
       path = here::here("outputs", "figures"),
       width = 168,
       height = 80,
       units = "mm",
       dpi = 300,
)


# Figure 3 - elevation ----------------------------------------------------
fig3 <- ggplot(
  elev_summary,
  aes(
    x = horizon,
    y = mean_elevation,
    color = species,
    linetype = futureScenario,
    shape = shape,
    group = interaction(species, futureScenario, year)
  )
) +
  scale_shape_manual(
    values = c("Present" = 15, "2050" = 16, "2090" = 17),
    breaks = c("2050", "2090")   # só anos na legenda
  )+
  geom_line(linewidth = 0.5) +
  geom_point(size = 1.5) +
  geom_ribbon(
    aes(ymin = lwr, ymax = upr, fill = species),
    alpha = .15,
    color = NA
  ) +
  facet_wrap(~ spGroup) +
  theme_bw(base_size = 8) +
  labs(
    x = "Time horizon",
    y = "Elevation (mean ± CI)",
    color = "Species",
    linetype = "Scenario",
    shape = "Year",
  ) +
  guides(fill = "none",
         linetype = guide_legend(order = 1),
         shape    = guide_legend(order = 2),
         color    = guide_legend(order = 3))+
  scale_color_discrete(
    labels = function(x) paste0("<i>", x, "</i>")
  ) +
  theme(
    legend.position='none',
    legend.direction='vertical',
    legend.text = element_markdown(size = 7),
    legend.title = element_text(size = 8),
    strip.text = element_text(size = 8),    # títulos dos facets
    axis.text = element_text(size = 7),
    axis.title = element_text(size = 8)
  )

fig3
?theme
ggsave("Fig 3.png",
       device = png,
       plot = fig3,
       path = here::here("outputs", "figures"),
       width = 168,
       height = 80,
       units = "mm",
       dpi = 300,
)




# -------------------------------------------------------------------------
#legend <- cowplot::get_legend(fig1)
#plot(legend)

fig <- cowplot::plot_grid(fig1, fig2, fig3, legend,
                  nrow=4, ncol=1, #align = 'h',
                  labels = c("(a)","(b)","(c)"), label_size=10, 
                  vjust=2, hjust=-0.5)

fig 


legend <- cowplot::get_legend(fig3)
plot(legend)

ggsave("Fig_trends_bottom.png",
       plot = fig,
       path = here::here("outputs", "figures"),
       width = 80,
       height = 240,
       units = "mm",
       dpi = 300,
       bg="transparent"
)
