library(maptools)
library(mapproj)
library(rgeos)
library(rgdal)
library(RColorBrewer)
library(ggplot2)
library(grid)
library(gridExtra)

# Load SpatialPolygonsDataFrame of US State borders
# Move Alaska and Hawaii to be in plot area
LoadAndCreateMap<-function(dir,dsn,layer){
  
  # Change to location of file
  setwd(dir)
  
  # Load file
  states50 <-readOGR(dsn=dsn,layer=layer)
  
  # Change projection
  states50 <- spTransform(states50, CRS("+proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m +no_defs"))
  states50@data$id <- rownames(states50@data)
  
  # Rotate, Scale and Move Alaska
  # Alaska State ID Number=02
  # https://www.census.gov/geo/reference/docs/state.txt
  alaska <- states50[states50$STATEFP=="02",]
  alaska <- elide(alaska, rotate=-50)
  alaska <- elide(alaska, scale=max(apply(bbox(alaska), 1, diff)) / 2.2)
  alaska <- elide(alaska, shift=c(-2100000, -2500000))
  # Force projection to be that of original plot
  proj4string(alaska) <- proj4string(states50)
  
  # Rotate and Move Hawaii, ID=15
  hawaii <- states50[states50$STATEFP=="15",]
  hawaii <- elide(hawaii, rotate=-35)
  hawaii <- elide(hawaii, shift=c(5400000, -1400000))
  # Force projection to be that of original plot
  proj4string(hawaii) <- proj4string(states50)
  
  # remove old Alaska/Hawaii and also DC, virgin islands, mariana, puerto rico, etc...
  states48 <- states50[!states50$STATEFP %in% c("02", "15", "72","66","60","69","74","78",'11'),]
  states.final <- rbind(states48, alaska, hawaii)

  return(states.final)
}

# Load Census data
# Combine data with shapefile
CombineDataWithMap<-function(dir,shapefile){

  # Change to location of data
  setwd(dir)
  
  # Load Census Data
  x<-read.table('ACS_13_1YR_DP02_with_ann.csv',sep=',',header=T,na.strings = 'N')[-435,]
  x$state = tolower(x$GEO.display.label)
  
  # Only keep a small number of columns that we're intersted in.
  # and rename our columns
  trim.data=with(x, data.frame(state=state,
                               pctArab=HC03_VC187,
                               pctEnglish=HC03_VC191,
                               pctGerman=HC03_VC194,
                               pctIrish=HC03_VC197,
                               pctItalian=HC03_VC198,
                               pctPolish=HC03_VC201,
                               pctSwedish=HC03_VC208,
                               pctRussian=HC03_VC203,
                               pctAmerican=HC03_VC186,
                               pctCzech=HC03_VC188,
                               pctDanish=HC03_VC189,
                               pctDutch=HC03_VC190,
                               pctFrench=HC03_VC192,
                               pctGreek=HC03_VC195,
                               pctScotchIrish=HC03_VC204,
                               pctScottish=HC03_VC205,
                               pctDanish=HC03_VC189))
  
  # Make lowercase for merging
  shapefile@data$state=tolower(shapefile@data$NAME)
  shapefile@data<-merge(shapefile@data,
                        trim.data,
                        by='state',
                        sort = F,
                        all.x=T)
  
  # Transform Shapefile into data.frame for plotting.
  states.plotting<-fortify(shapefile,region='state')
  
  # Merge Data with fortified shapefile
  state.dat<-merge(states.plotting,
                   trim.data,
                   by.x='id',
                   by.y = 'state',
                   sort = F,
                   all.x=T)   # Keep all coordinate lines
  
  state.dat=state.dat[order(state.dat$order),]
  
  # Return fortified map and shapefile with census data in @data.
  return(list(dat=state.dat,map=shapefile))
}

# Create text annotations
MakeTextLabels<-function(map,var){
  
  # Find centroid of each state
  trueCentroids = gCentroid(map,byid=TRUE)
  
  # Initialize annotations to be the centroid
  text.labels<-as.data.frame(trueCentroids@coords)
  
  # Add USPS ID codes to each state
  text.labels$id <- map@data$STUSPS
  
  # Format part of the annotation value 
  text.labels[,var] <- format(map@data[,var],digits=3)
  
  # Combine parts of what we want to be shown
  text.labels$lab<-paste(text.labels$id,text.labels[,var],sep='\n')
  
  # Lastly, we manually fix some states
  text.labels[,c('orig.x','orig.y')]=text.labels[,c('x','y')]
  
  # Shift some coordinates slightly
  text.labels[text.labels$id == 'CA',c('x')]<-c(-1750000)
  text.labels[text.labels$id == 'LA',c('x')]<-c(710000)
  text.labels[text.labels$id == 'FL',c('x')]<-c(1800000)
  text.labels[text.labels$id == 'KY',c('x')]<-c(1350224)
  text.labels[text.labels$id == 'NY',c('x')]<-c(1970000)
  text.labels[text.labels$id == 'SC',c('x')]<-c(1770000)
  text.labels[text.labels$id == 'VA',c('x')]<-c(1860000)
  text.labels[text.labels$id == 'WV',c('x')]<-c(1660000)
  text.labels[text.labels$id == 'MI',c('x','y')]<-c(1250000,-100000)
  
  # Move some labels to outside the state
  text.labels[text.labels$id == 'HI',c('x','y')]<-c(-388560.8,-2258325)
  text.labels[text.labels$id == 'NH',c('x','y')]<-c(1900000,600000)
  text.labels[text.labels$id == 'VT',c('x','y')]<-c(1900000,400000)
  
  text.labels[text.labels$id == 'MA',c('x','y')]<-c(2600000,292375.55)
  text.labels[text.labels$id == 'RI',c('x','y')]<-c(2600000,100000)
  text.labels[text.labels$id == 'CT',c('x','y')]<-c(2600000,-100000)
  text.labels[text.labels$id == 'NJ',c('x','y')]<-c(2600000,-300000)
  text.labels[text.labels$id == 'DE',c('x','y')]<-c(2600000,-500000)
  text.labels[text.labels$id == 'MD',c('x','y')]<-c(2600000,-700000)
  text.labels[text.labels$id == 'DC',c('x','y')]<-c(2600000,-900000)
  
  text.labels$outside=text.labels$id %in% c('MA','RI','CT','NJ','DE','MD','DC','NH','VT')
  
  # Draw a connecting line from label to state
  OFFSET=100000
  text.labels$init.horiz.x.start=text.labels$x
  text.labels$init.horiz.y.start=text.labels$y
  text.labels$init.horiz.x.end=text.labels$x-OFFSET
  text.labels[text.labels$id %in% c('VT','NH'),]$init.horiz.x.end=text.labels[text.labels$id %in% c('VT','NH'),]$x+1.1*OFFSET
  text.labels$init.horiz.y.end=text.labels$y
  text.labels$second.horiz.y.start=text.labels$y
  text.labels$second.horiz.x.start=text.labels$x-OFFSET
  text.labels[text.labels$id %in% c('VT','NH'),]$second.horiz.x.start=text.labels[text.labels$id %in% c('VT','NH'),]$x+1.1*OFFSET
  text.labels$second.horiz.x.end=text.labels$orig.x
  text.labels$second.horiz.y.end=text.labels$orig.y
  
  text.labels$size=6
  text.labels[text.labels$id %in% c('WV'),'size']=4.5
  text.labels$color='white'
  text.labels[text.labels$id %in% c('MA','RI','CT','NJ','DE','MD','DC','HI','NH','VT'),'color']='black'
  
  return(text.labels)
}

# Plot single variable choropleth of US States using ggplot
# Text annotations optional
make.plot<-function(map,shape,add.labels=FALSE,var,
                    label,
                    colors=(brewer.pal(9, "Blues")[4:9])){
  
  
  map$plotVar<-map[,var]
  map.plot<-ggplot(map)+
    theme_bw()+
    ggtitle(label)+
    theme(axis.ticks = element_blank(), 
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.title.x= element_blank(),
          axis.title.y= element_blank(),
          plot.background = element_blank(),
          plot.title = element_text(size=25, face="bold"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          legend.title=element_text(size=18,face='bold'),
          legend.text=element_text(size=18,face='bold'))+
          #legend.key.size = unit(x = 0.3,y=0.8,units = 'cm'))+
    geom_polygon(aes(x = long, y = lat,group = group,fill=plotVar))+
    geom_polygon(aes(x = long, y = lat,group = group),fill=NA,colour='black',size=0.1)+
    scale_fill_gradientn('%',colours=colors,guide =  guide_legend(title='%',
                                                                  title.position='top',
                                                                  label.position='right',
                                                                  title.hjust=0.2,
                                                                  label.hjust=0.5,
                                                                  keywidth=2,
                                                                  keyheight=3))
  if(add.labels){
    text.labels<-MakeTextLabels(shape,var)
    map.plot<-map.plot+geom_text(data = text.labels,
                                 aes(x = x,y=y,label=lab),
                                 colour=text.labels$color,
                                 size=text.labels$size)+
      geom_segment(data=text.labels[text.labels$outside,],
                   aes(x = second.horiz.x.start,
                       y = second.horiz.y.start, 
                       xend = second.horiz.x.end, 
                       yend = second.horiz.y.end),size=1.3)
  }
  
  return(map.plot)
}



base_dir='/media/roberto/Main Storage/Documents/bertplot/code/R/Multiple State-wide Plots Of Ancestry/'

setwd(base_dir)
dsn='cb_2013_us_state_20m'
layer='cb_2013_us_state_20m'

map<-LoadAndCreateMap(base_dir,dsn,layer)
ret<-CombineDataWithMap(base_dir,map)
map.w.data<-ret$dat
map<-ret$map


irish<-make.plot(map=map.w.data,shape=map,add.labels=F,var='pctIrish',label='Irish')
german<-make.plot(map=map.w.data,shape=map,add.labels=F,var='pctGerman',label='German')
swedish<-make.plot(map=map.w.data,shape=map,add.labels=F,var='pctSwedish',label='Swedish')
italian<-make.plot(map=map.w.data,shape=map,add.labels=F,var='pctItalian',label='Italian')
english<-make.plot(map=map.w.data,shape=map,add.labels=F,var='pctEnglish',label='English')
arab<-make.plot(map=map.w.data,shape=map,add.labels=F,var='pctArab',label='Arab')
russian<-make.plot(map=map.w.data,shape=map,add.labels=F,var='pctRussian',label='Russian')
greek<-make.plot(map=map.w.data,shape=map,add.labels=F,var='pctGreek',label='Greek')
polish<-make.plot(map=map.w.data,shape=map,add.labels=F,var='pctPolish',label='Polish')
french<-make.plot(map=map.w.data,shape=map,add.labels=F,var='pctFrench',label='French')
dutch<-make.plot(map=map.w.data,shape=map,add.labels=F,var='pctDutch',label='Dutch')
danish<-make.plot(map=map.w.data,shape=map,add.labels=F,var='pctDanish',label='Danish')


png('full-grid.png',height=1200,width=1200)
grid.arrange(german,english,irish,italian,french,polish,swedish,dutch,danish,arab,russian,greek,ncol=3)
dev.off()

for(i in c('Irish','German','Swedish','Italian','English','Arab','Russian',"Greek",'Polish','French','Dutch','Danish')){
    png(paste(i,'.png',sep=''),height=1000,width=1000)
    print(make.plot(map=map.w.data,shape=map,add.labels=T,var=paste('pct',i,sep=''),label=i))
    dev.off()
}

