if (.Platform$OS.type=="windows") try(Sys.setlocale("LC_CTYPE","English_United States.utf8"), silent=TRUE)
options(device=function(...) grDevices::pdf(NULL))
script <- normalizePath(sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)),winslash="/")
scripts <- dirname(script)
data <- file.path(dirname(scripts),"figure_data")
args <- commandArgs(trailingOnly=TRUE)
out <- if(length(args)) args[1] else file.path(dirname(dirname(dirname(scripts))),"review-output","china")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
suppressPackageStartupMessages({library(ggplot2);library(sf);library(dplyr);library(grid);library(ggmapcn)})
g <- readRDS(file.path(data,"China_10km_plot_grid.rds"))
tr <- read.csv(file.path(data,"Fig_1b_grid_trends.csv"))
j <- match(g$TID,tr$grid_id)
stopifnot(!anyNA(j),!anyDuplicated(j),length(j)==96649)
g$HQ_Mean <- tr$mean_hq[j]
g$Slope <- tr$slope_per_year[j]
g$P_Value <- tr$p_value[j]
g$EN <- tr$province_code[j]
dd <- st_drop_geometry(g)
prov <- read.csv(file.path(data,"Fig_1_province_means.csv"))
ink<-"#000000";font<-"Helvetica"
hq_palette<-c("#F5F5F0","#C7E9C0","#74C476","#238B45","#00441B")
make_map<-function(metric){
 p<-ggplot()+geom_mapcn(fill="white",color=NA)+
    geom_buffer_cn(mainland_dist=50000,fill="#EFF1F1")+
    geom_buffer_cn(mainland_dist=20000,fill="white")
 if(metric=="HQ"){
   p<-p+geom_sf(data=g,aes(fill=HQ_Mean),colour=NA)+
     scale_fill_gradientn(colours=hq_palette,limits=c(0,1),guide="none")
 }else{
   p<-p+geom_sf(data=g[is.na(g$P_Value),],aes(fill=Slope),colour=NA,alpha=.1)+
     geom_sf(data=g[!is.na(g$P_Value)&g$P_Value>=.05,],fill="#D9D9D9",colour=NA)+
     geom_sf(data=g[!is.na(g$P_Value)&g$P_Value<.05,],aes(fill=Slope),colour=NA)+
     scale_fill_gradient2(low="#D73027",mid="#FFFFBF",high="#4575B4",midpoint=0,
       limits=c(-.01,.01),oob=scales::squish,guide="none")
 }
 p+geom_boundary_cn(province_color="#929B9F",province_size=.13,
    mainland_color="#4F5B61",mainland_size=.3,coastline_color="#8799A8",coastline_size=.12,
    ten_segment_line_color="#6A777E",ten_segment_line_size=.15,
    SAR_boundary_color="#929B9F",SAR_boundary_size=.12,
    undefined_boundary_color="#6A777E",undefined_boundary_size=.15)+
   scale_x_continuous(breaks=seq(80,130,10),labels=function(x)paste0(x,"\u00b0E"))+
   scale_y_continuous(breaks=seq(20,50,10),labels=function(x)paste0(x,"\u00b0N"))+
   theme_minimal(base_size=7,base_family=font)+
   theme(text=element_text(family=font,colour=ink),axis.text=element_text(size=6.5,colour=ink),
     axis.title=element_blank(),panel.grid=element_blank(),legend.position="none",
     plot.margin=margin(0,0,0,0),plot.background=element_rect(fill=NA,colour=NA))
}
boxdata<-dd[!is.na(dd$EN),];boxdata$EN<-factor(boxdata$EN,levels=prov$EN)
boxdata$province_mean<-prov$mean[match(boxdata$EN,prov$EN)]
n<-nrow(prov);labs<-data.frame(EN=factor(prov$EN,levels=prov$EN),label=prov$EN,i=seq_len(n))
angle<-90-360*(labs$i-.5)/n
labs$hjust<-ifelse(angle < -90,1,0);labs$angle<-ifelse(angle < -90,angle+180,angle)
labs$radius<-ifelse(labs$label=="NMG",1.22,1.12)
p_radial<-ggplot(boxdata,aes(EN,HQ_Mean,fill=province_mean))+
 geom_boxplot(linewidth=.14,outlier.size=.06,outlier.alpha=.16,colour="#657078")+
 geom_text(data=labs,aes(x=EN,y=radius,label=label,angle=angle,hjust=hjust),inherit.aes=FALSE,
   family=font,size=5.2/ggplot2::.pt,colour=ink)+
 coord_polar(theta="x",start=0,clip="off")+
 scale_y_continuous(limits=c(-.10,1.25),expand=c(0,0))+
 scale_fill_gradientn(colours=hq_palette,limits=c(0,1),guide="none")+
 theme_void(base_family=font)+theme(plot.margin=margin(2,2,2,2,"mm"),plot.background=element_rect(fill=NA,colour=NA))
# R-only ridge inset; all density estimates and summary statistics are reused.
# Adjustable appearance settings.
RIDGE_COLOURS<-c('2000'='#BFD9C8','2010'='#78AE95','2020'='#2E775D')
RIDGE_HEIGHT<-.70
REFERENCE_COLOUR<-'#858585'
ridge_s<-read.csv(file.path(data,'Fig_1_ridge_summaries.csv'),colClasses=c(year='character'))
ridge_d<-read.csv(file.path(data,'Fig_1_ridge_density.csv'))
ridge_peak<-max(ridge_d$density)
ridge_polys<-do.call(rbind,lapply(seq_len(nrow(ridge_s)),function(i){
 d<-ridge_d[ridge_d$year==as.numeric(ridge_s$year[i]),];base<-ridge_s$ridge_base[i]
 data.frame(year=ridge_s$year[i],x=c(d$hq,rev(d$hq)),y=c(base+d$density/ridge_peak*RIDGE_HEIGHT,rep(base,nrow(d))))
}))
ridge_s$tip<-vapply(seq_len(nrow(ridge_s)),function(i){d<-ridge_d[ridge_d$year==as.numeric(ridge_s$year[i]),]
 ridge_s$ridge_base[i]+approx(d$hq,d$density,xout=ridge_s$mean[i])$y/ridge_peak*RIDGE_HEIGHT
},numeric(1))
stopifnot(nrow(ridge_s)==3L,all(ridge_s$n==96649L))
make_ridge_plot<-function(compact=TRUE){
 ggplot()+geom_polygon(data=ridge_polys,aes(x,y,group=year,fill=year),colour='black',linewidth=.20)+
  annotate('segment',x=ridge_s$mean[1],xend=ridge_s$mean[1],y=.94,yend=3.78,
    colour=REFERENCE_COLOUR,linewidth=.23,linetype='22')+
  geom_segment(data=ridge_s,aes(x=mean,xend=mean,y=ridge_base,yend=ridge_base+.65),
    colour='black',linewidth=.24,linetype='22')+
  geom_point(data=ridge_s,aes(mean,tip),shape=21,fill='black',colour='white',stroke=.18,size=if(compact)1.0 else 1.4)+
  geom_text(data=ridge_s,aes(x=.99,y=ridge_base+.76,label=sprintf('%.3f',mean)),
    hjust=1,vjust=.5,family='Helvetica',colour='black',size=(if(compact)5.4 else 7)/ggplot2::.pt)+
  scale_fill_manual(values=RIDGE_COLOURS)+
  scale_x_continuous(limits=c(0,1),breaks=c(0,.5,1),labels=function(x)sprintf('%.1f',x),expand=expansion(mult=c(.02,.02)))+
  scale_y_continuous(limits=c(.9,3.91),breaks=3:1,labels=c('2000','2010','2020'),expand=c(0,0))+
  labs(x='Habitat quality',y=NULL,title=NULL)+
  theme_classic(base_size=if(compact)5.5 else 7,base_family='Helvetica')+
  theme(text=element_text(colour='black'),axis.text=element_text(colour='black',size=if(compact)5.2 else 6.3),
    axis.title=element_text(size=if(compact)5.5 else 6.5),
    plot.title=element_text(size=if(compact)6.2 else 7.5,face='bold',margin=margin(b=3)),
    axis.line.y=element_blank(),axis.ticks.y=element_blank(),axis.line.x=element_line(colour='black',linewidth=.22),
    axis.ticks.x=element_line(colour='black',linewidth=.22),axis.ticks.length=unit(.8,'mm'),
    plot.margin=margin(1,1,1,1,'mm'),legend.position='none',
    panel.background=element_rect(fill=NA,colour=NA),plot.background=element_rect(fill=NA,colour=NA))
}

p_time<-make_ridge_plot(TRUE)
pdf(NULL,width=183/25.4,height=106/25.4,family=font)
builds<-lapply(c("HQ","Trend"),function(m)ggplot_build(make_map(m)))
maps<-lapply(builds,ggplot_gtable)
bars<-list(builds[[1]]$plot$scales$get_scales("fill")$map(seq(0,1,length.out=128)),
 builds[[2]]$plot$scales$get_scales("fill")$map(seq(-.01,.01,length.out=128)))
radial<-ggplotGrob(p_radial);timeline<-ggplotGrob(p_time)
dev.off()
draw_figure<-function(measure=FALSE){
 grid.newpage();measured<-list()
 txt<-function(s,x,y,fs=7,bold=FALSE,just="left",col=ink){
  grid.text(s,unit(x,"mm"),unit(y,"mm"),just=just,
    gp=gpar(fontfamily=font,fontsize=fs,col=col,fontface=if(bold)"bold" else "plain"))
 }
 for(j in 1:2){
  left<-c(3,95)[j]
  txt(c("a","b")[j],left,102,10,TRUE,col="black")
  txt(c("Mean habitat quality (2000\u20132020)","Annual trend (2000\u20132020)")[j],left+6,102,7,TRUE)
  pushViewport(viewport(x=unit(left,"mm"),y=unit(2,"mm"),width=unit(85,"mm"),height=unit(96,"mm"),
    just=c("left","bottom"),clip="off",name=paste0("map-",j)))
  grid.draw(maps[[j]])
  if(measure){
    grid.force();v<-grid.ls(viewports=TRUE,grobs=FALSE,recursive=TRUE,print=FALSE)
    wh<-which(grepl("^panel[.]",v$name)&grepl(paste0("map-",j),v$vpPath,fixed=TRUE));stopifnot(length(wh)>=1)
    k<-tail(wh,1);path<-strsplit(paste(v$vpPath[k],v$name[k],sep="::"),"::",fixed=TRUE)[[1]]
    upViewport(0);downViewport(do.call(vpPath,as.list(path[path!="ROOT"])))
    lo<-deviceLoc(unit(0,"npc"),unit(0,"npc"),valueOnly=TRUE)
    hi<-deviceLoc(unit(1,"npc"),unit(1,"npc"),valueOnly=TRUE)
    measured[[j]]<-list(id=paste0("map-",j),bbox_pt=c(lo$x,lo$y,hi$x,hi$y)*72,
      grid_id="maps",row_start=0L,row_stop=1L,col_start=j-1L,col_stop=j)
    upViewport(0)
  }else upViewport()
  if(j==1){
    pushViewport(viewport(x=unit(left-2.5,"mm"),y=unit(3,"mm"),width=unit(48,"mm"),height=unit(48,"mm"),
      just=c("left","bottom"),clip="off"))
    grid.draw(radial);upViewport()
  }else{
    pushViewport(viewport(x=unit(left+2,"mm"),y=unit(9,"mm"),width=unit(40,"mm"),height=unit(32,"mm"),
      just=c("left","bottom"),clip="off"))
    grid.draw(timeline);upViewport()
  }
  xx<-left+if(j==1)15 else 19;yy<-95.5;ww<-30
  for(k in 1:128)grid.rect(x=unit(xx+(k-.5)/128*ww,"mm"),y=unit(yy,"mm"),
    width=unit(ww/128+.01,"mm"),height=unit(2,"mm"),gp=gpar(fill=bars[[j]][k],col=NA))
  ticks<-if(j==1)c("0","0.5","1") else c("\u22120.01","0","+0.01")
  for(k in 1:3)txt(ticks[k],xx+(k-1)/2*ww,yy-3.4,6.5,just="center")
  txt(if(j==1)"HQ" else "Slope/yr",xx-2,yy,6.2,just="right")
  if(j==2){
    grid.rect(x=unit(xx+ww+6.5,"mm"),y=unit(yy,"mm"),width=unit(5,"mm"),height=unit(2,"mm"),gp=gpar(fill="#D9D9D9",col=NA))
    txt("p \u2265 0.05",xx+ww+10.5,yy,5.8)
  }
 }
 invisible(measured)
}
prefix <- file.path(out,"Fig_1")
cairo_pdf(paste0(prefix,".pdf"),width=183/25.4,height=106/25.4,family=font)
draw_figure();dev.off()
ragg::agg_png(paste0(prefix,".png"),width=183,height=106,units="mm",res=300,background="white")
draw_figure();dev.off()
