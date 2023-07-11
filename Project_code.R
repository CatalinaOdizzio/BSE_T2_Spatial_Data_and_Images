library(SpatialEpi)
library(sp)
library(spdep)
library(INLA)
library(leaflet)
library(DT)
library(ggplot2)

#Loading data: Lip cancer - Scotland from SpatialEpi

data(scotland)
data(scotland_sf)

#Exploring dataset

names(scotland)
names(scotland_sf)

df <- scotland$data #Creating dataframe with the data

head(df)

df %>% datatable() #Creating table to include in the report

#Map
map <- scotland$spatial.polygon

proj4string(map) <- "+proj=tmerc +lat_0=49 +lon_0=-2
+k=0.9996012717 +x_0=400000 +y_0=-100000 +datum=OSGB36
+units=km +no_defs"

##Transform to latitude/longitude
map <- spTransform(map, CRS("+proj=longlat +datum=WGS84 +no_defs"))

##Variables
df$areaid <- df$county.names
df$Y <- df$cases 
df$E <- df$expected
df$covariate <- df$AFF
df$SIR <- df$Y/df$E #Creating SIR variable

##Adding data 
rownames(df) <- df$county
map <- SpatialPolygonsDataFrame(map, df, match.ID = TRUE)
map$areaid <- map$county.names
order <- match(map$areaid, df$areaid)
map@data <- df[order, c("areaid", "Y", "E", "covariate", "SIR")]
head(map@data)

plot(map)

#Mapping SIR
lm <- leaflet(map) %>% addTiles() #Creating leaflet map

palette <- colorNumeric(palette = "YlGn", domain = map$SIR) #Defining color palette

lm %>% addPolygons(color = "black", weight = 1, fillColor = ~palette(SIR), fillOpacity = 0.6) %>%
  addLegend(pal = palette, values = ~SIR, opacity = 0.5, title = "SIR", position = "bottomright") #Adding polygon shapes to the map and coloring them based on the SIR variable

#Modelling

##Neighborhood matrix
nm <- poly2nb(map) #Generating neighborhood matrix from map
head(nm)

nb2INLA("map.adj", nm) #Converting neighborhood matrix into a format that can be used with INLA
mg <- inla.read.graph(filename = "map.adj") 

##Creating u and v
map$u <- 1:nrow(map@data)
map$v <- 1:nrow(map@data)

##Inference
formula <- Y ~ covariate + f(u, model = "besag", graph = mg, scale.model = TRUE) + f(v, model = "iid") #Defining regression. Response variable is Y. Predictors: AFF, u and v

results <- inla(formula, family = "poisson", data = map@data, E = E, control.predictor = list(compute = TRUE)) #INLA on the previous regression, where the response variable is modeled as a Poisson process

##Results
summary(results)

##Distribution of AFF coefficient

marginal <- inla.smarginal(results$marginals.fixed$covariate)
marginal <- data.frame(marginal)
ggplot(marginal, aes(x = x, y = y)) + geom_line() + labs(x = expression(beta[1]), y = "Density") +
  geom_vline(xintercept = 0, col = "black") + theme_bw() + ggtitle("AFF Coefficient Distribution") +
  theme(plot.title = element_text(hjust = 0.5))

##Fitted values
results$summary.fitted.values %>% 
  round(2) %>% 
  setNames(c("Mean", "Standard Deviation", "Lower Quantile" , "Median", "Upper Quantile", "Mode")) %>% 
  datatable(extensions = c('Scroller', 'FixedColumns'))

#Mapping disease risk

map$RR <- results$summary.fitted.values[, "mean"]
map$LL <- results$summary.fitted.values[, "0.025quant"]
map$UL <- results$summary.fitted.values[, "0.975quant"]

palette <- colorNumeric(palette = "YlOrRd", domain = map$RR)

labels <- sprintf("<strong> %s </strong> <br/> Observed: %s <br/> Expected: %s <br/>
                  AFF proportion: %s <br/>SIR: %s <br/>RR: %s (%s, %s)",
                  map$areaid, map$Y,  round(map$E, 2),  map$covariate, round(map$SIR, 2),
                  round(map$RR, 2), round(map$LL, 2), round(map$UL, 2)) %>%
  lapply(htmltools::HTML)


leaflet(map) %>% 
  addTiles() %>%
  addPolygons(color = "black", weight = 1, fillColor = ~palette(RR),  fillOpacity = 0.7,
              highlightOptions = highlightOptions(weight = 4),
              label = as.character(labels),
              labelOptions = labelOptions(style = list("font-weight" = "normal", padding = "3px 8px"),
                                          textsize = "15px", direction = "auto")) %>%
  addLegend(pal = palette, values = ~RR, opacity = 0.5, title = "RR", position = "bottomright")

#Summary table

risk <- data.frame(County = map$areaid,
                              Observed = map$Y, 
                              Expected = round(map$E, 2), 
                              AFF = map$covariate, 
                              SIR = round(map$SIR, 2),
                              RR = round(map$RR, 2), 
                              LL = round(map$LL, 2), 
                              UL = round(map$UL, 2))

datatable(risk, extensions = c('Scroller', 'FixedColumns'))
