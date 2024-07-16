### Mapping using ggOceanmaps
##https://mikkovihtakari.github.io/ggOceanMaps/articles/ggOceanMaps.html

library(ggOceanMapsData)
library(ggOceanMaps)
library(sf)
library(ggspatial)

##
raw.metadata <- read.csv("metadata/CoreMetaData.csv")
metadata <- raw.metadata[raw.metadata$CoreID!="GC06"&raw.metadata$CoreID!="PC022",]


basemap(limits = c(-25, -12, 63, 67), bathymetry = TRUE)

limits = c(-25, -12, 63, 67)


pdf("Maps/OverviewV1.pdf",height = 6,width = 9)
basemap(limits = c(-25, -12, 63, 67), bathy.style = "rcb",rotate = TRUE)+
  #basemap(limits = c(-25, -12, 63, 67),rotate = TRUE)+
    ggspatial::geom_spatial_point(aes(x = Lon, y = Lat),data = metadata, color = c("black"),crs = 4326,cex=2.5)+
    ggspatial::geom_spatial_point(aes(x = Lon, y = Lat),data = metadata, color = c("dodgerblue","darkred"),crs = 4326,cex=2)+
    ggspatial::geom_spatial_text(aes(x = Lon, y = Lat+0.1,label=CoreID),data = metadata,cex=2)+
  xlab("Longitude")+
  ylab("Latitude")+
  annotation_scale(location = "bl")+ 
  annotation_north_arrow(location = "tr", which_north = "true", height = unit(0.9, "cm"), width = unit(0.9, "cm"))
dev.off()

pdf("Maps/OverviewV1mini.pdf",height = 3.5,width = 4.2)
basemap(limits = c(-25, -13, 63, 67), bathy.style = "rcb",rotate = TRUE)+
  ggspatial::geom_spatial_point(aes(x = Lon, y = Lat),data = metadata, color = c("black"),crs = 4326,cex=4.5)+
  ggspatial::geom_spatial_point(aes(x = Lon, y = Lat),data = metadata, color = c("#52868B","#282E69"),crs = 4326,cex=4)+
  ggspatial::geom_spatial_text(aes(x = Lon, y = Lat+0.3,label=c("GC01","PC19")),data = metadata,cex=3)+
  xlab("")+
  ylab("")+
  theme(legend.position="right",legend.box.spacing=unit(0.5, "lines"),
        legend.key.size = unit(0.3, 'cm'), #change legend key size
        legend.key.height = unit(0.5, 'cm'), #change legend key height
        legend.key.width = unit(0.3, 'cm'), #change legend key width
        legend.title = element_text(size=8), #change legend title font size
        legend.text = element_text(size=6))+
  #theme(legend.position.inside = c(1,1))+
  annotation_scale(location = "bl")+ 
  annotation_north_arrow(location = "tl", which_north = "true", height = unit(0.9, "cm"), width = unit(0.9, "cm"))
dev.off()



basemap(limits = c(-25, -12, 63, 67), bathy.style = "rcb",rotate = TRUE)+
  ggspatial::geom_spatial_point(aes(x = Lon, y = Lat),data = metadata, color = "red")+
  xlab("Longitude")+
  ylab("Latitude")




 ggspatial::geom_spatial_point(aes(x = Lon, y = Lat),data = metadata, color = "red")
 
 
 
basemap(limits = c(-21.5, -16.5, 65.5, 67), bathymetry = TRUE,rotate = TRUE)+
  geom_point(data = transform_coord(metadata), aes(x = Lon, y = Lat), color = "red")
?basemap

