# Cell press guidelines ==================================================================
# Define variables for figures formatting

# Note Cell guidelines:
# Always include/embed fonts and only use Helvetica or Arial fonts.
#Line or stroke width should not be narrower than half a point. 
# Line weights should range from 0.35 to 1.5 pt.
# Gray fills should be kept at least 20% different from other fills and no lighter than 10% or darker than 80%
# Each figure should be able to fit on a single 8.5 x 11 inch page.
# Please do not send figure panels as individual files.

# We use three standard widths for figures:
# 1 column, 85 mm;
one_col = 85
col_unit = 'mm'
# 1.5 column, 114 mm;
one_half_col = 114
# and 2 column, 174 mm (the full width of the page).
two_col = 174
# Define a size converter
# When defining size (geom_text, annotate...), it's in 'mm' so multiply by this to convert in points ('pt')
size_converter = 0.352777778 

# Define a general theme function
# Basically, take the theme_bw and divide everything by two
# Remove things pertaining to color, etc
theme_cell <- theme(
    axis.line = element_blank(),
    axis.text = element_text(margin = unit(0.2, "lines")),
    axis.text.x = element_text(size = unit(7 * 0.85, 'pt'), lineheight = 0.7, vjust = .5),
    axis.text.y = element_text(size = 7 * 0.85, lineheight = 0.7, hjust = .5),
    axis.ticks = element_line(size = 0.2),
    axis.title.x = element_text(size = 7, vjust = .5),
    axis.title.y = element_text(size = 7, angle = 90, vjust = 0.25),
    axis.ticks.length = unit(0.15, "lines"),
    # axis.title = element_text(size = unit(6, 'pt'), margin = margin(0, 0, 0, 0, unit = 'pt')),
    # axis.text = element_text(size = unit(5, 'pt'), margin = margin(0, 0, 0, 0, unit = 'pt')),
    
    legend.key.size = unit(.5, "lines"),
    legend.text = element_text(size = 7 * 0.85, margin = margin(0, 0, 0, 2, unit = 'pt')),
    legend.title = element_text(size = 7 * 0.9, face = "bold", hjust = 0,
                                margin = margin(0, 0, 2, 0, unit = 'pt')),
    legend.position = "right",
    # legend.text = element_text(size = unit(6, 'pt'), margin = margin(0, 0, 0, 0, unit = 'pt')),
    # legend.title = element_text(size = unit(6, 'pt'), margin = margin(0, 0, 0, 0, unit = 'pt')),
    # legend.key.size = unit(2, 'mm'),
    # legend.key.height = unit(2, 'mm'),
    # legend.key.width = unit(2, 'mm'),
    # legend.box.margin = margin(0, 0, 0, 0, unit = 'pt'),
    # legend.margin = margin(0, 0, 0, 0, unit = 'pt'),
    # legend.box.spacing = unit(0, 'pt'),
    # legend.spacing = unit(0, 'pt')),
    
    panel.grid.major = element_line(size = 0.5*size_converter),
    panel.grid.minor = element_line(size = 0.3*size_converter),
    panel.spacing = unit(0.1, "lines"),
    
    strip.text.x = element_text(size = 7 * 0.85),
    strip.text.y = element_text(size = 7 * 0.85, angle = -90),
    # strip.text = element_text(size = unit(6, 'pt'), margin = margin(0, 0, 0, 0, unit = 'pt')),
    
    plot.title = element_text(size = 7),
    plot.margin = unit(c(.5, .5, 0.25, 0.25), "lines")
    # plot.title = element_text(size = unit(7, 'pt'), margin = margin(0, 0, 0, 0, unit = 'pt')),
    # plot.subtitle = element_text(size = unit(6, 'pt'), margin = margin(0, 0, 0, 0, unit = 'pt')),
    #plot.margin = margin(0, 0, 0, 0, unit = 'pt')
  )

