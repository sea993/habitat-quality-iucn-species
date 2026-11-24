# Rebuild the national panels at the same physical scale as the species panels.
# Statistical definitions and geographic layers follow the original R map.
suppressPackageStartupMessages({library(ggplot2); library(sf); library(dplyr); library(grid); library(ggmapcn); library(jsonlite)})
map_ink <- '#2D3D46'
map_grid <- readRDS(file.path(data,"China_10km_plot_grid.rds"))
tr <- read.csv(file.path(data,"Fig_1b_grid_trends.csv"))
ind <- read.csv(file.path(data,"Fig_3_grid_mean_indices.csv"))
j <- match(map_grid$TID,ind$TID)
k <- match(map_grid$TID,tr$grid_id)
stopifnot(!anyNA(j),!anyNA(k),!anyDuplicated(j),!anyDuplicated(k))
map_grid$EN <- tr$province_code[k]
map_grid$ANCI_mean <- ind$ANCI_mean[j]
map_grid$APCI_mean <- ind$APCI_mean[j]
map_df <- st_drop_geometry(map_grid)
map_metrics <- c('ANCI','APCI')
map_palettes <- c(ANCI='YlOrRd',APCI='YlGn')
province_stats <- lapply(map_metrics,function(m) {
  z <- map_df[!is.na(map_df$EN),]
  z %>% group_by(EN) %>% summarise(mean=mean(.data[[paste0(m,'_mean')]]),n=n(),.groups='drop') %>% arrange(mean)
}); names(province_stats)<-map_metrics
map_plot <- function(m) {
  ggplot() +
    geom_mapcn(fill='white',color=NA) +
    geom_buffer_cn(mainland_dist=50000,fill='#EFF1F1') +
    geom_buffer_cn(mainland_dist=20000,fill='white') +
    geom_sf(data=map_grid,aes(fill=.data[[paste0(m,'_mean')]]),colour=NA) +
    geom_boundary_cn(province_color='#929B9F',province_size=.13,mainland_color='#4F5B61',mainland_size=.3,
      coastline_color='#8799A8',coastline_size=.12,ten_segment_line_color='#6A777E',ten_segment_line_size=.15,
      SAR_boundary_color='#929B9F',SAR_boundary_size=.12,undefined_boundary_color='#6A777E',undefined_boundary_size=.15) +
    scale_fill_distiller(palette=map_palettes[m],direction=1,limits=c(0,100),guide='none') +
    scale_x_continuous(breaks=seq(80,130,10),labels=function(x)paste0(x,'\u00b0E')) +
    scale_y_continuous(breaks=seq(20,50,10),labels=function(x)paste0(x,'\u00b0N')) +
    theme_minimal(base_size=7,base_family='Helvetica') +
    theme(text=element_text(family='Helvetica',colour=map_ink),axis.text=element_text(size=6.5,colour=map_ink),
      axis.title=element_blank(),panel.grid=element_blank(),
      legend.position='none',plot.margin=margin(0,0,0,0),plot.background=element_rect(fill=NA,colour=NA))
}
radial_plot <- function(m) {
  s <- province_stats[[m]]
  dd <- map_df[!is.na(map_df$EN),];dd$EN<-factor(dd$EN,levels=s$EN)
  dd$province_mean<-s$mean[match(dd$EN,s$EN)]
  n<-nrow(s); labs<-data.frame(EN=factor(s$EN,levels=s$EN),label=s$EN,i=seq_len(n))
  angle<-90-360*(labs$i-.5)/n
  labs$hjust<-ifelse(angle < -90,1,0)
  labs$angle<-ifelse(angle < -90,angle+180,angle)
  labs$label_radius<-ifelse(labs$label=='NMG',122,112)
  ggplot(dd,aes(EN,.data[[paste0(m,'_mean')]],fill=province_mean)) +
    geom_boxplot(linewidth=.14,outlier.size=.06,outlier.alpha=.16,colour='#657078') +
    geom_text(data=labs,aes(x=EN,y=label_radius,label=label,angle=angle,hjust=hjust),inherit.aes=FALSE,
      family='Helvetica',size=5/ggplot2::.pt,colour=map_ink) +
    coord_polar(theta='x',start=0,clip='off') +
    scale_y_continuous(limits=c(-10,125),expand=c(0,0)) +
    scale_fill_distiller(palette=map_palettes[m],direction=1,limits=c(0,100),guide='none') +
    theme_void(base_family='Helvetica') + theme(plot.margin=margin(2,2,2,2,'mm'),plot.background=element_rect(fill=NA,colour=NA))
}
# Build once; no separate image rescaling enters the final composition.
pdf(NULL,width=183/25.4,height=106/25.4,family='Helvetica',useDingbats=FALSE)
national_builds <- lapply(map_metrics,function(m)ggplot_build(map_plot(m)))
national_grobs <- lapply(national_builds,function(b)black_text_grob(ggplot_gtable(b)))
map_colourbars <- lapply(national_builds,function(b)b$plot$scales$get_scales('fill')$map(seq(0,100,length.out=128)))
radial_grobs <- lapply(map_metrics,function(m)black_text_grob(ggplotGrob(radial_plot(m))))
invisible(dev.off())
draw_map_sections <- function(measure=FALSE,newpage=FALSE) {
  if(newpage)grid.newpage()
  base_path <- current.vpPath()
  measured<-list()
  txt<-function(s,x,y,fs=7,bold=FALSE,just='left',col=map_ink)grid.text(s,unit(x,'mm'),unit(y,'mm'),just=just,
    gp=gpar(fontfamily='Helvetica',fontsize=fs,col=black_text_colour(col),fontface=if(bold)'bold' else 'plain'))
  for (j in 1:2) {
    left<-c(3,95)[j]; m<-map_metrics[j]
    txt(c('c','d')[j],left,102,10,TRUE,col='black')
    txt(paste0(m,'  |  2000–2020 mean'),left+6,102,7,TRUE)
    pushViewport(viewport(x=unit(left,'mm'),y=unit(2,'mm'),width=unit(85,'mm'),height=unit(96,'mm'),
      just=c('left','bottom'),clip='off',name=paste0('map-',m)))
    grid.draw(national_grobs[[j]])
    if(measure) {
      grid.force(); v<-grid.ls(viewports=TRUE,grobs=FALSE,recursive=TRUE,print=FALSE)
      wh<-which(grepl('^panel[.]',v$name)&grepl(paste0('map-',m),v$vpPath,fixed=TRUE));stopifnot(length(wh)>=1)
      k<-tail(wh,1);full<-paste(v$vpPath[k],v$name[k],sep='::');path<-strsplit(full,'::',fixed=TRUE)[[1]]
      upViewport(0);downViewport(do.call(vpPath,as.list(path[path!='ROOT'])))
      lo<-deviceLoc(unit(0,'npc'),unit(0,'npc'),valueOnly=TRUE);hi<-deviceLoc(unit(1,'npc'),unit(1,'npc'),valueOnly=TRUE)
      measured[[j]]<-list(id=paste0('map-',m),bbox_pt=c(lo$x,lo$y,hi$x,hi$y)*72,grid_id='maps',row_start=0L,row_stop=1L,col_start=j-1L,col_stop=j)
      upViewport(0)
      if(!is.null(base_path))downViewport(base_path)
    } else upViewport()
    pushViewport(viewport(x=unit(left-2.5,'mm'),y=unit(3,'mm'),width=unit(48,'mm'),height=unit(48,'mm'),just=c('left','bottom'),clip='off'))
    grid.draw(radial_grobs[[j]]);upViewport()
    pal<-map_colourbars[[j]]
    xx<-left+15;yy<-94;ww<-30
    for(k in 1:128)grid.rect(x=unit(xx+(k-.5)/128*ww,'mm'),y=unit(yy,'mm'),width=unit(ww/128+.01,'mm'),height=unit(2,'mm'),
      gp=gpar(fill=pal[k],col=NA))
    for(t in c(0,50,100))txt(as.character(t),xx+t/100*ww,yy-3.4,6.5,just='center')
    txt('Index',xx-2,yy,6.5,just='right')
  }
  invisible(measured)
}
