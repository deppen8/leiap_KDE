library(here)
library(rgdal)
library(rgeos)
library(maptools)
library(spatstat)
library(RColorBrewer)
library(sp)
library(raster)
library(rasterVis)
library(gridExtra)
library(ggplot2)
library(grid)
library(scales)
library(ggthemes)
library(ggrepel)
library(ggsn)
library(ggpolypath)
library(spatialEco)

kde_per_production <- function(prod_file, title){
  pts <- read.csv(here('data', 'prods', prod_file))
  coordinates(pts) <- ~ Easting + Northing
  
  ch <- shapefile(here('data', 'geo', 'convex_hull_2018.shp'))
  muni_raw <- readOGR(here('data', 'geo', 'SS_municipio.shp'))
  muni <- gSimplify(muni_raw, tol=0.01, topologyPreserve=TRUE)
  muni <- fortify(muni)
  sites <- readOGR(here('data', 'geo', 'known_sites.shp'))
  sites.clip <- sites[ch, ]
  sites.clip <- as.data.frame(sites.clip)
  colnames(sites.clip) <- c("Name", "Northing", "Easting")
  
  fields <- readOGR(here('data', 'geo', 'geo_fields_dissolve.shp'))
  fields <- fields[ch, ]
  
  ch.poly <- as(ch, "SpatialPolygons")
  bb <- bbox(ch.poly)
  w <- as(ch.poly, 'owin')
  ch <- fortify(ch)
  
  # pts.ppp <- ppp(pts[,'Easting'], pts[,'Northing'], window=w)
  
  # main density calculations
  sigmas = c(50, 100, 150, 300)
  densities <- vector("list", 4)
  for (i in 1:length(sigmas)){
    # sp.kde(x, y, bw, newdata, n, standardize = FALSE, scale.factor)
    d <- sp.kde(x=pts, y=pts@data[,3], bw=sigmas[i], n=1000, standardize=TRUE)
    densities[[i]] <- d
  }
  
  # for (i in 1:length(sigmas)){
  #   d <- density(pts.ppp, weights=(pts[,5]), sigma=sigmas[i], eps=5.0, positive=TRUE)
  #   densities[[i]] <- d
  # }
  
  # Prep some plotting elements
  ## color for pts
  pal = brewer.pal(n = 9, name = "YlOrRd")
  ## scalebar locations
  ### kdes
  x1 <- min(ch$long)
  y1 <- min(ch$lat) + (max(ch$lat)-min(ch$lat))*.05
  x2 <- x1+1000
  y2 <- y1
  ### inset
  x4 <- max(muni$long)
  y4 <- min(muni$lat) + (max(muni$lat)-min(muni$lat))*.05
  x3 <- x4 - 2000
  y3 <- y4
  
  ## Title
  title <- title
  
  ## Total weight
  wt <- sum(pts@data[,3])
  
  ## Total frags
  frags <- sum(pts@data[,2])
  
  # Create plots
  ## Text info plot
  data_text = paste("pes total: ",wt,"g\n",
                    "fragments: ",frags)
  text_plot <- ggplot() + geom_point() + xlim(0, 1) + ylim(0,1) +
    annotate("text", x = 0.5, y=0.6, size=5, label = title) +
    annotate("text", x = 0.5, y=0.50, size=4, label="2014-2018") +
    annotate("text", x = 0.5, y=0.30, size=4, label = data_text) +
    theme_map() +
    theme(axis.line=element_blank(),
          axis.title.x = element_blank(), 
          axis.ticks.x = element_blank(), 
          axis.text.x = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank(), 
          axis.text.y = element_blank())
  
  ## Son Servera area and convex hull
  inset <- ggplot(muni) + ggtitle(paste('àrea analitzada')) +
    geom_polygon(data=muni, aes(x=long, y=lat, group=group), alpha=0.1) +  # Son Servera outline
    geom_polygon(data=ch, aes(x=long, y=lat, group=group),  # convex hull used in analysis
                 color='darkgray', alpha=0.7, fill=NA, size=0.4) +
    geom_polypath(data=fields, aes(x=long, y=lat, group=group), alpha=0.5,  # buffered fields
                  color=NA, fill='gray75', size=0.1) +
    # geom_point(data=sites.clip, aes(x=Easting, y=Northing),  # site locations
    #            size=1, alpha=0.7, shape=3) +
    geom_segment(aes(x=x3, y=y3, xend=x4, yend=y4),   # scale bar
                 color='lightgray', size=0.3) +
    geom_text(aes(x=mean(c(x3,x4)), y=y3+200, label='2km'),  # scale bar label
              size=3, color='lightgray') +
    coord_equal() +
    theme_map() +
    theme(plot.title = element_text(hjust = 0.5),
          axis.line=element_blank(),
          axis.title.x = element_blank(), 
          axis.ticks.x = element_blank(), 
          axis.text.x = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank(), 
          axis.text.y = element_blank())
  
  ## KDE plots
  kde <- vector("list", 4)
  for (i in 1:length(densities)){
    d <- densities[[i]]
    # r <- raster(d)
    # dm <- mask(r, fields)
    # raster_spdf <- as(dm, "SpatialPixelsDataFrame")
    dm <- mask(d, fields)
    raster_spdf <- as(dm, "SpatialPixelsDataFrame")
    raster_df <- as.data.frame(raster_spdf)
    colnames(raster_df) <- c("kde", "x", "y")
    
    p <- ggplot() + ggtitle(paste('r = ',sigmas[i],'m')) +
      geom_tile(data=raster_df, aes(x=x, y=y, fill=kde), alpha=0.8) +  # KDE result
      # geom_polygon(data=fields, aes(x=long, y=lat, group=group),  # fields outline
      #              alpha=1.0, color='black', fill=NA, size=0.1) +
      geom_point(data=sites.clip, aes(x=Easting, y=Northing),  # site locations
                 size=1, alpha=0.7, shape=3) +
      geom_text_repel(data=sites.clip, aes(x=Easting, y=Northing, label=Name),  # site labels
                      segment.color='lightgray', size=3, alpha=0.5) +
      geom_segment(aes(x=x1, y=y1, xend=x2, yend=y2),  # scale bar
                   color='gray', size=0.3) +
      geom_text(aes(x=mean(c(x1, x2)), y=y1+150, label='1km'),  # scale bar text
                size=3, color='gray') +
      scale_fill_gradientn(colours= pal,                       # density legend
                           breaks=seq(0, 1, by=0.1),
                           labels=sprintf("%0.1f", seq(0,1,by=0.1)),
                           guide = guide_legend(title=NULL, direction="vertical",
                                                ncol=1, label.position="right",
                                                label.vjust=1.1, reverse=TRUE
                           )) +
      coord_equal() +
      theme_map() +
      theme(legend.justification=c(0,0), legend.position=c(1,0),
            legend.key.height = unit(0.9, "cm"),
            legend.background = element_rect(colour="#00000000", fill="#00000000"),
            plot.title = element_text(hjust = 0.5),
            axis.line=element_blank(),
            axis.title.x = element_blank(), 
            axis.ticks.x = element_blank(), 
            axis.text.x = element_blank(),
            axis.title.y = element_blank(),
            axis.ticks.y = element_blank(), 
            axis.text.y = element_blank())
    
    kde[[i]] <- p
  }
  
  # Save plot outputs in a grid
  png(filename = here('outputs', 'tests', paste(title, '.png', sep='')), 
      width = 14, height = 8, units = 'in', res = 100)
  grid.arrange(kde[[1]], kde[[2]], text_plot, kde[[3]], kde[[4]], inset, ncol=3, nrow=2)
  dev.off()
}


for (fn in list.files(here('data', 'prods'))[1]){
  title <- sub('.csv', '', fn, fixed=TRUE)  # drop .csv
  kde_per_production(fn, title)
}