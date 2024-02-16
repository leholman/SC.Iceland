## Map with coordinates ##
#Creating a map of relevant sites at Iceland for WP4
#Colours #f1a340 (orange), #f7f7f7 (white), #998ec3 (purple)

#Load the packages
library("ggplot2")
theme_set(theme_bw())
library("sf")
library("rnaturalearth")
library("rnaturalearthdata")
library("rgeos")
library("maps")
library("ggmap")
library("ggsn")
#First we define the frame of the map by attitude and longitude.
location <- c(left=-19.7, bottom=65.5, right=-17.2, top=66.8)

#Now we'll source a stamen map, the options are "terrain-background", "toner-lite" and "watercolor". 
map <- (get_stamenmap(location, zoom = 9, maptype = "terrain-background"))
ggmap(map)

#Now we create point data or sites
site_archeological <- data.frame(longitude= c(-18.26193077, -18.45846275, -18.02457025, -17.99972715), latitude = c(65.74753883, 65.67201469, 66.54877891, 66.53454765))
site_core <- data.frame(longitude = c(-19.3034, -17.4199, -17.700265, -17.700281, -19.50553), latitude = c(66.3009, 66.331, 66.551726, 66.551721, 66.500766))
#Map the sites
ggmap(map)+
  ggtitle("North East Iceland")+
  xlab("Longitude") + ylab("Latitude")+
  geom_point(data = site_archeological, aes(x = longitude, y = latitude), size = 3,
             shape = 23, fill = "coral3")+
  geom_point(data = site_core, aes(x = longitude, y = latitude), size = 3,
             shape = 23, fill = "darkgoldenrod1")+
  theme(text = element_text(size = 15))+
  coord_equal()+
  ggsn::scalebar(dist = 10, model = "WGS84", dist_unit = "km", st.dist = 100, height = 0.3, transform =TRUE, x.min = -21.25, x.max = -19.25,  y.min = 65.51, y.max = 65.61)+
  ggsn::north(x.min = -17.75, x.max = -17.25, 
              y.min = 65.96, y.max = 66.86, scale = 0.35)+
  coord_fixed(xlim = c(-19.7, -17.2), ratio = 2/1)



#Map of Iceland
location <- c(left=-24.996982, bottom=63.256773, right=-12.531783, top=66.790210)
map <- (get_stamenmap(location, zoom = 9, maptype = "terrain-background"))
ggmap(map)

ggmap(map)+
  ggtitle("Iceland")+
  xlab("Longitude") + ylab("Latitude")+
  theme(text = element_text(size = 15))+
  coord_equal()+
  coord_fixed(xlim = c(-24.996982, -12.531783), ratio = 2/1)


scalebar(map, dist = 100, dist_unit = "km",
         transform = TRUE, model = "WGS84")
