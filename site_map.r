# load required libraries
library(maps)
library(sf)
library(sp)

# georeferencing blue marble data
# gdal_translate -a_ullr -180 90 180 -90 \
# -a_srs EPSG:4326 world.topo.bathy.200408.3x5400x2700.jpg tmp.tif

df <- read.table("sites.csv", sep = ",", header = TRUE)

# function to slice and dice a map and convert it to an sp() object
maps2sp = function(xlim, ylim, l.out = 100, clip = TRUE) {
  stopifnot(require(maps))
  m = map(xlim = xlim, ylim = ylim, plot = FALSE, fill = TRUE)
  p = rbind(cbind(xlim[1], seq(ylim[1],ylim[2],length.out = l.out)),
            cbind(seq(xlim[1],xlim[2],length.out = l.out),ylim[2]),
            cbind(xlim[2],seq(ylim[2],ylim[1],length.out = l.out)),
            cbind(seq(xlim[2],xlim[1],length.out = l.out),ylim[1]))
  LL = CRS("+init=epsg:4326")
  IDs <- sapply(strsplit(m$names, ":"), function(x) x[1])
  stopifnot(require(maptools))
  m = map2SpatialPolygons(m, IDs=IDs, proj4string = LL)
  bb = SpatialPolygons(list(Polygons(list(Polygon(list(p))),"bb")), proj4string = LL)

  if (!clip)
    m
  else {
    stopifnot(require(rgeos))
    gIntersection(m, bb)
  }
}

# set colours for map grid
grid.col.light = rgb(0.5,0.5,0.5,0.8)
grid.col.dark = rgb(0.5,0.5,0.5)

# coordinate systems
polar = CRS("+init=epsg:3995")
longlat = CRS("+init=epsg:4326")

if(!exists("r")){
  # read in the raster map and
  # set the extent, crop to extent and reproject to polar
  #r = raster::brick("blue_marble.tif")
  r = raster::brick("blue_marble.tif")
  e = raster::extent(c(-180,180,55,90))
  r_crop = raster::crop(r,e)

  # traps NA values and sets them to 1
  r_crop[is.na(r_crop)] = 1
  r_polar = raster::projectRaster(r_crop, crs = polar, method = "bilinear")

  # some values are not valid after transformation
  # (rgb range = 1 - 255) set these back to 1
  # as they seem to be the black areas
  r_polar[r_polar < 1 ] = 1
}

pts=SpatialPoints(rbind(c(-180,55),c(0,55),c(180,85),c(180,85)), CRS("+init=epsg:4326"))
gl = gridlines(pts, easts = seq(-180,180,30), norths = seq(50,85,10), ndiscr = 100)
gl.polar = spTransform(gl, polar)

pts=SpatialPoints(rbind(c(-180,55),c(0,55),c(180,80),c(180,80)), CRS("+init=epsg:4326"))
ll = SpatialLines(list(Lines(Line(cbind(seq(-180,180,0.5),rep(55,721))), ID="outer")), CRS("+init=epsg:4326"))

# make a table of the sites to evaluate (should be done in a separate file for later use, eventually)
sites.df = cbind(df$lon,df$lat)

# create a spatial points object and transform to polar coords
sites.latlong = SpatialPoints(sites.df,longlat)
sites.polar = spTransform(sites.latlong,polar)

# crop a map (make the x component a bit larger not to exclude)
# some of the eastern islands (the centroid defines the bounding box)
m = maps2sp(c(-180,200),c(55,90),clip = TRUE)

pdf("map.pdf")

# plotting routine
# set margins
par(mar=rep(1,4))

# plot the grid, to initiate the area
# plotRGB() overrides margin settings in default plotting mode
plot(spTransform(gl, polar), lwd=2, lty=2,col="white")

# plot the blue marble raster
raster::plotRGB(r_polar, add = TRUE)

# plot grid lines
lines(spTransform(gl, polar), add = TRUE, lwd=2, lty=2,col=grid.col.light)

# plot outer margin of the greater circle
lines(spTransform(ll, polar), lwd = 3, lty = 1, col=grid.col.dark)

# plot continent outlines
plot(spTransform(m, polar), lwd = 1,
     lty = 1, col = "transparent", border=grid.col.dark, add = TRUE)

# plot longitude labels
l = labels(gl.polar, longlat, side = 1)
l$pos = NULL
text(l, cex = 1, adj = c( 0.5, 2 ),  col = "black")

# plot latitude labels
l = labels(gl.polar, longlat, side = 2)
l$srt = 0
l$pos = NULL
text(l, cex = 1, adj = c(1.2, -1), col = "white")

# site locations and names
points(sites.polar,pch=17, cex=1.2, lwd = 2, col = "yellow")
text(sites.polar@coords,
     labels=df$sitename,
     pos = df$pos,
     col="white",
     cex = 1
     )
dev.off()
