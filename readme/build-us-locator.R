# Rebuild the original U.S. locator map from archived vector data.
library(data.table)
library(sf)
library(ggplot2)
ink <- "#18262D"
meta <- fread("readme/data/us_forests/sites.csv", colClasses=c(shot="character"))
map_furniture <- function(b, km, dark = ink, north_left = FALSE) {
  width <- diff(b[1:2]); height <- diff(b[3:4])
  x <- b[2] - width * .045 - km * 1000
  y <- b[3] + height * .044
  nx <- if (north_left) b[1] + width * .079 else b[2] - width * .079
  ny <- b[4] - height * .072
  list(annotate("rect", xmin = x - width*.016, xmax = x + km*1000 + width*.015,
                ymin = y - height*.017, ymax = y + height*.051, fill = "white", alpha = .90),
       annotate("segment", x = x, xend = x + km*1000, y = y, yend = y, colour = dark, linewidth = .5),
       annotate("segment", x = x, xend = x + km*500, y = y, yend = y,
                colour = dark, linewidth = .92),
       annotate("text", x = x + km*500, y = y + height*.027, label = paste(km, "km"),
                colour = dark, size = 2.45),
       annotate("rect", xmin = nx - width*.043, xmax = nx + width*.043,
                ymin = ny - height*.056, ymax = ny + height*.048, fill = "white", alpha = .90),
       annotate("segment", x = nx, xend = nx, y = ny - height*.039, yend = ny + height*.004,
                colour = dark, linewidth = .4, arrow = grid::arrow(length = grid::unit(1.7, "mm"), type = "closed")),
       annotate("text", x = nx, y = ny + height*.029, label = "N", size = 2.55, colour = dark, fontface = "bold"))
}

points <- st_transform(st_as_sf(meta, coords = c("lon", "lat"), crs = 4326), 5070)
xy <- st_coordinates(points)
points$dx <- fcase(points$landscape == "Northern Maine", -220000,
                   points$landscape == "Olympic Peninsula", 450000,
                   points$landscape == "Everglades", -180000,
                   default = 155000)
points$dy <- fcase(points$landscape == "Olympic Peninsula", -180000,
                   points$landscape == "Everglades", 140000,
                   default = 45000)
points$x_label <- xy[, 1] + points$dx
points$y_label <- xy[, 2] + points$dy
points$label <- sprintf("%s  %s", points$id, points$landscape)


states <- st_read("readme/data/us_forests/states.geojson",quiet=TRUE)
regions <- st_read("readme/data/us_forests/ecoregions.geojson",quiet=TRUE)
country <- st_union(st_make_valid(states))
id <- unique(st_drop_geometry(st_join(points, regions[, "NA_L1NAME"], left = FALSE))$NA_L1NAME)
regions$region <- ifelse(regions$NA_L1NAME %in% id, regions$NA_L1NAME, "Other ecoregions")
colours <- c('EASTERN TEMPERATE FORESTS'='#99BFD3',
 'MARINE WEST COAST FOREST'='#8AB3A6','NORTH AMERICAN DESERTS'='#DBCA92',
 'NORTHERN FORESTS'='#D5A3A2','NORTHWESTERN FORESTED MOUNTAINS'='#ABA5CB',
 'MEDITERRANEAN CALIFORNIA'='#B3C88F','GREAT PLAINS'='#C9C2A5',
 'TROPICAL WET FORESTS'='#B3C88F','Other ecoregions'='#E9E7E1')
box <- unname(st_bbox(country)[c("xmin", "xmax", "ymin", "ymax")])
box <- box + c(-2e5, 1.5e5, -1e5, 2e5)
map <- ggplot() +
  geom_sf(data = regions, aes(fill = region), colour = NA) +
  geom_sf(data = states, fill = NA, colour = "#FFFFFF", linewidth = .27) +
  geom_sf(data = country, fill = NA, colour = "#708286", linewidth = .34) +
  geom_segment(data = cbind(st_drop_geometry(points), x = xy[, 1], y = xy[, 2]),
               aes(x, y, xend = x + dx * .92, yend = y_label), colour = "#4F5657", linewidth = .26) +
  geom_sf(data = points, colour = "#112D32", fill = "white", shape = 21, size = 2.45, stroke = .75) +
  geom_label(data = st_drop_geometry(points), aes(x_label, y_label, label = label), hjust = .5,
             size = 3.0, label.padding = grid::unit(.10, "lines"), linewidth = 0, fill = "#FFFEFE") +
  scale_fill_manual(values = colours, name = "Level I ecoregions",
                    labels = function(x) tools::toTitleCase(tolower(x))) +
  coord_sf(xlim = box[1:2], ylim = box[3:4], expand = FALSE, datum = NA) +
  map_furniture(box, 750, north_left = TRUE) +
  guides(fill = guide_legend(ncol = 1, override.aes = list(colour = "#9DA7A8"))) +
  labs(title = "A   GEDI and USGS 3DEP comparison sites") +
  theme_void(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 12),
        legend.position = "right", legend.title = element_text(face = "bold"),
        legend.text = element_text(size = 10))

ggsave("readme/fig-us-forest-sites.png", map, width=10.5, height=4.5, dpi=400, bg="white")