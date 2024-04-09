### Mapping using ggOceanmaps
##https://mikkovihtakari.github.io/ggOceanMaps/articles/ggOceanMaps.html

library(ggOceanMapsData)
library(ggOceanMaps)
library(sf)
library(ggspatial)

##
raw.metadata <- read.csv("metadata/CoreMetaData.csv")
metadata <- raw.metadata[raw.metadata$CoreID!="GC06",]


basemap(limits = c(-25, -12, 63, 67), bathymetry = TRUE)

limits = c(-25, -12, 63, 67)


pdf("Maps/OverviewV1.pdf",height = 6,width = 9)
basemap(limits = c(-25, -12, 63, 67), bathy.style = "rcb",rotate = TRUE)+
  #basemap(limits = c(-25, -12, 63, 67),rotate = TRUE)+
    ggspatial::geom_spatial_point(aes(x = Lon, y = Lat),data = metadata, color = "red",crs = 4326)+
  #ggspatial::geom_spatial_text(aes(x = Lon, y = Lat+0.04,label=CoreID),data = metadata, color = "darkred",cex=2)+
  xlab("Longitude")+
  ylab("Latitude")+
  annotation_scale(location = "bl")+ 
  annotation_north_arrow(location = "tr", which_north = "true", height = unit(0.9, "cm"), width = unit(0.9, "cm"))
dev.off()


basemap(limits = c(-25, -12, 63, 67), bathy.style = "rcb",rotate = TRUE)+
  ggspatial::geom_spatial_point(aes(x = Lon, y = Lat),data = metadata, color = "red")+
  xlab("Longitude")+
  ylab("Latitude")




 ggspatial::geom_spatial_point(aes(x = Lon, y = Lat),data = metadata, color = "red")
 
 
 
basemap(limits = c(-21.5, -16.5, 65.5, 67), bathymetry = TRUE,rotate = TRUE)+
  geom_point(data = transform_coord(metadata), aes(x = Lon, y = Lat), color = "red")
?basemap

