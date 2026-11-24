if (.Platform$OS.type=="windows") try(Sys.setlocale("LC_CTYPE","English_United States.utf8"), silent=TRUE)
options(device=function(...) grDevices::pdf(NULL))
script <- normalizePath(sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)),winslash="/")
scripts <- dirname(script)
data <- file.path(dirname(scripts),"figure_data")
args <- commandArgs(trailingOnly=TRUE)
out <- if(length(args)) args[1] else file.path(dirname(dirname(dirname(scripts))),"review-output","china")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
suppressPackageStartupMessages({library(ggplot2);library(dplyr);library(grid)})
source(file.path(scripts,'figure_text.R'),encoding='UTF-8')
read_source <- function(n)read.csv(file.path(data,paste0('Fig_2_',n,'.csv')),check.names=FALSE)
annual<-read_source('annual_all_species'); groups<-read_source('annual_group_summaries')
period<-read_source('species_period_means');national<-read.csv(file.path(dirname(data),'statistics','annual_national.csv'))
heat<-read_source('crossclassified_means');directions<-read_source('group_endpoint_proportions')
cats<-c('CR','EN','VU','NT','LC');taxa<-c('Amphibians','Reptiles','Mammals','Plants')
metrics<-c('Q_mean','Q_5pct','R_low');all_groups<-c(cats,taxa)
direction_col<-c(Increasing='#3B75AF',Decreasing='#C95849',Unchanged='#C5C9CC')
cat_col<-c(CR='#A84343',EN='#CF7146',VU='#C7A54A',NT='#91A75A',LC='#4F957B')
tax_col<-c(Amphibians='#467FAD',Reptiles='#C89449',Mammals='#B45B68',Plants='#59947A')
blue<-'#387F9B';ink<-'#2D3D46';gray<-'#6B7379'
theme_set(theme_classic(base_size=7,base_family='Helvetica')+theme(
 axis.text=element_text(size=6.4,color=ink),axis.title=element_blank(),axis.line=element_line(linewidth=.23,color='#7D8990'),
 axis.ticks=element_line(linewidth=.2,color='#7D8990'),axis.ticks.length=unit(1.5,'pt'),
 legend.position='none',panel.grid=element_blank(),plot.margin=margin(1.5,1.5,1.5,1.5),plot.background=element_rect(fill=NA,color=NA)))
fmt<-function(x,m)if(m=='R_low')sprintf('%.2f',x) else sprintf('%.3f',x)
annual_change<-annual%>%group_by(metric)%>%arrange(year,.by_group=TRUE)%>%mutate(change=100*(mean/mean[year==2000]-1))%>%ungroup()
national_change<-national%>%arrange(year)%>%mutate(change=100*(mean/mean[year==2000]-1))
stopifnot(nrow(annual_change)==63,all(annual$n==1820),nrow(period)==5460)
change_plot<-function(m){
 dd<-annual_change[annual_change$metric==m,];last<-dd[dd$year==2020,]
 lim<-switch(m,Q_mean=c(-6,3.1),Q_5pct=c(-14,3.2),R_low=c(-4,36))
 ticks<-switch(m,Q_mean=c(-6,-3,0),Q_5pct=c(-10,-5,0),R_low=c(0,15,30))
 p<-ggplot(dd,aes(year,change))+geom_ribbon(aes(ymin=pmin(change,0),ymax=pmax(change,0)),fill=blue,alpha=.16)+
  geom_hline(yintercept=0,linewidth=.25,color='#9AA4A8')+geom_line(color=blue,linewidth=.6)+
  geom_point(data=last,color=blue,size=1.4)+
  geom_text(data=last,aes(label=sprintf('%+.1f%%',change)),hjust=1.1,vjust=if(m=='R_low')-.7 else 1.5,size=6.8/ggplot2::.pt,color=blue,fontface='bold')+
  annotate('text',x=2000,y=lim[2],label=paste0(if(m=='Q_mean')'Species: ' else '',fmt(dd$mean[dd$year==2000],m),' to ',fmt(last$mean,m)),hjust=0,vjust=1,size=6.3/ggplot2::.pt,color=ink)+
  scale_x_continuous(breaks=c(2000,2010,2020),limits=c(1999.5,2020.5),expand=c(0,0))+
  scale_y_continuous(limits=lim,breaks=ticks,labels=function(x)paste0(x,'%'),expand=c(0,0))
 if(m=='Q_mean')p<-p+annotate('text',x=2000,y=2,label=sprintf('National: %.3f to %.3f',national$mean[national$year==2000],national$mean[national$year==2020]),hjust=0,vjust=1,size=6.1/ggplot2::.pt,color=gray)+
  geom_line(data=national_change,aes(year,change),color=gray,linewidth=.45,linetype='22')+
  geom_text(data=national_change[national_change$year==2020,],aes(label=sprintf('%+.1f%%',change)),hjust=1.1,vjust=-2,size=6.1/ggplot2::.pt,color=gray)
 p
}
row_data<-function(type,m){
 lev<-if(type=='category')cats else taxa
 dd<-groups[groups$group_type==type & groups$metric==m,]
 dd$y<-length(lev)+1-match(dd$group,lev);dd
}
endpoints<-function(type,m,first=FALSE){
 lev<-if(type=='category')cats else taxa;cols<-if(type=='category')cat_col else tax_col
 dd<-row_data(type,m);e0<-dd[dd$year==2000,];e1<-dd[dd$year==2020,]
 pair<-merge(e0,e1,by=c('group','y'),suffixes=c('0','1'));pair$change<-100*(pair$mean1/pair$mean0-1)
 est<-dd[dd$year %in% c(2000,2020),];est$yp<-est$y+ifelse(est$year==2000,.045,-.045)
 xlim<-if(m=='Q_mean')c(.62,.86) else c(.29,.75)
 stopifnot(all(est$mean-est$sem>xlim[1]),all(est$mean+est$sem<xlim[2]))
 p<-ggplot()+
  geom_segment(data=pair,aes(x=mean0,xend=mean1,y=y+.045,yend=y-.045,color=group),linewidth=.8,alpha=.65)+
  geom_segment(data=est,aes(x=mean-sem,xend=mean+sem,y=yp,yend=yp,color=group),linewidth=.3)+
  geom_point(data=est[est$year==2000,],aes(mean,yp,color=group),shape=21,fill='white',size=2,stroke=.5)+
  geom_point(data=est[est$year==2020,],aes(mean,yp,color=group,fill=group),shape=21,size=2,stroke=.5)+
  geom_text(data=pair,aes(x=xlim[2]+diff(xlim)*.025,y=y,label=sprintf('%+.1f%%',change)),hjust=0,size=6.1/ggplot2::.pt,color=ink)+
  scale_color_manual(values=cols)+scale_fill_manual(values=cols)+
  scale_x_continuous(breaks=if(m=='Q_mean')seq(.65,.85,.05) else seq(.3,.7,.1),labels=function(x)sprintf('%.2f',x),expand=c(0,0))+
  coord_cartesian(xlim=xlim,clip='off')+
  scale_y_continuous(limits=c(.45,length(lev)+.65),breaks=rev(seq_along(lev)),labels=if(first)lev else rep('',length(lev)),expand=c(0,0))+
  theme(axis.line.y=element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_text(size=6.4))
 p
}
cloud_geometry<-list();cloud_summary<-list()
raincloud<-function(type){
 lev<-if(type=='category')cats else taxa;cols<-if(type=='category')cat_col else tax_col
 field<-if(type=='category')'category' else 'class'
 dd<-period[period$metric=='R_low',];dd$group<-dd[[field]];dd$y<-length(lev)+1-match(dd$group,lev)
 dd$x<-log1p(dd$value);stopifnot(all(dd$value>=0 & dd$value<=100),nrow(dd)==1820)
 set.seed(20261023);dd$yp<-dd$y-runif(nrow(dd),.10,.29)
 bandwidth<-bw.nrd0(dd$x);L<-log1p(100);polys<-list();sums<-list()
 for(g in lev){
  v<-dd[dd$group==g,];den<-density(c(v$x,-v$x,2*L-v$x),bw=bandwidth,from=0,to=L,n=256)
  den$y<-3*den$y;norm<-sum(diff(den$x)*(head(den$y,-1)+tail(den$y,-1))/2);den$y<-den$y/norm
  yy<-unique(v$y);polys[[g]]<-data.frame(group=g,x=c(den$x,rev(den$x)),y=c(yy+.04+.35*den$y/max(den$y),rep(yy+.04,256)))
  q<-quantile(v$value,c(.25,.5,.75),names=FALSE)
  sums[[g]]<-data.frame(group=g,y=yy,n=nrow(v),zero=sum(v$value==0),q1=q[1],median=q[2],q3=q[3],bandwidth=bandwidth)
 }
 poly<-do.call(rbind,polys);ss<-do.call(rbind,sums)
 cloud_geometry[[type]]<<-list(n_points=nrow(dd),n_zero=sum(dd$value==0),density_domain=c(0,L),bandwidth=bandwidth)
 cloud_summary[[type]]<<-ss
 ggplot()+geom_vline(xintercept=log1p(c(0,1,5,20,100)),linewidth=.2,color='#E5E9EC')+
  geom_polygon(data=poly,aes(x,y,group=group,fill=group),alpha=.42,color=NA)+
  geom_point(data=dd,aes(x,yp,color=group),size=.35,alpha=.32,stroke=0)+
  geom_segment(data=ss,aes(x=log1p(q1),xend=log1p(q3),y=y,yend=y),linewidth=1.1,color='#35434B')+
  geom_point(data=ss,aes(x=log1p(median),y=y),shape=21,fill='white',color='#35434B',size=1.25,stroke=.35)+
  scale_color_manual(values=cols)+scale_fill_manual(values=cols)+
  scale_x_continuous(limits=c(-.08,L+.12),breaks=log1p(c(0,1,5,20,100)),labels=c('0','1','5','20','100'),expand=c(0,0))+
  scale_y_continuous(limits=c(.45,length(lev)+.65),breaks=rev(seq_along(lev)),labels=rep('',length(lev)),expand=c(0,0))+
  theme(axis.line.y=element_blank(),axis.ticks.y=element_blank())
}
heat_plot <- function(m,first=FALSE){
 dd<-heat[heat$metric==m,];dd$category<-factor(dd$category,levels=cats);dd$class<-factor(dd$class,levels=rev(taxa))
 dd$label<-if(m=='R_low')sprintf('%.1f',dd$value) else sprintf('%.2f',dd$value)
 pal<-if(m=='R_low')c('#FFF1DF','#EEB68D','#CB7858','#99423C') else c('#F0F7F5','#BDDCD7','#80BDB4','#418E8F','#245D74')
 tile_cols<-scales::gradient_n_pal(pal)(dd$value/if(m=='R_low')12 else 1)
 rgb<-grDevices::col2rgb(tile_cols)/255;linear<-ifelse(rgb<=.04045,rgb/12.92,((rgb+.055)/1.055)^2.4)
 lum<-colSums(linear*c(.2126,.7152,.0722));dd$text_col<-ifelse(lum>.179,'black','white')
 p<-ggplot(dd,aes(category,class,fill=value))+geom_tile(color='white',linewidth=.4)+geom_text(aes(label=label,color=text_col),size=6.5/ggplot2::.pt)+
  scale_color_identity()+scale_x_discrete(expand=c(0,0))+scale_y_discrete(expand=c(0,0))+
  theme(axis.line=element_blank(),axis.ticks=element_blank(),panel.grid=element_blank(),panel.grid.major.y=element_blank(),
        axis.text.y=if(first)element_text(size=6.3) else element_blank(),legend.position='bottom',legend.title=element_blank(),
        legend.text=element_text(size=6),legend.margin=margin(0),legend.key.height=unit(2,'mm'),legend.box.spacing=unit(.5,'mm'))+
  guides(fill=guide_colorbar(barwidth=unit(37,'mm'),barheight=unit(1.8,'mm'),ticks=FALSE,frame.colour=NA))
 if(m=='R_low')p+scale_fill_gradientn(colors=c('#FFF1DF','#EEB68D','#CB7858','#99423C'),limits=c(0,12),breaks=c(0,6,12))
 else p+scale_fill_gradientn(colors=c('#F0F7F5','#BDDCD7','#80BDB4','#418E8F','#245D74'),limits=c(0,1),breaks=c(0,.5,1))
}
y_pos<-setNames(c(10,9,8,7,6,4.7,3.7,2.7,1.7),all_groups)
bar_data<-directions[directions$group_type!='all',]
bar_data$group<-factor(bar_data$group,levels=all_groups);bar_data$direction<-factor(bar_data$direction,levels=names(direction_col))
bar_data<-bar_data%>%group_by(group,metric)%>%arrange(direction,.by_group=TRUE)%>%mutate(xmax=cumsum(count)/sum(count)*100,xmin=lag(xmax,default=0),xmid=(xmin+xmax)/2)%>%ungroup()
bar_data$y<-unname(y_pos[as.character(bar_data$group)])
bar_data$label<-ifelse(bar_data$pct>=12,sprintf('%.0f',bar_data$pct),'')
bar_data$text_col<-ifelse(bar_data$direction=='Unchanged','#344047','white')
counts<-sapply(all_groups,function(g)unique(bar_data$n[as.character(bar_data$group)==g]))
group_labels<-paste0(all_groups,'  ',counts)
stacked<-function(m,first=FALSE){
 dd<-bar_data[bar_data$metric==m,]
 ggplot(dd)+geom_rect(aes(xmin=xmin,xmax=xmax,ymin=y-.34,ymax=y+.34,fill=direction),linewidth=0)+
  geom_text(aes(x=xmid,y=y,label=label,color=text_col),size=6.4/ggplot2::.pt,fontface='bold')+
  scale_fill_manual(values=direction_col,drop=FALSE)+scale_color_identity()+
  scale_x_continuous(breaks=c(0,50,100),labels=c('0','50','100%'),expand=c(0,0))+coord_cartesian(xlim=c(0,100),clip='off')+
  scale_y_continuous(limits=c(.9,10.6),breaks=unname(y_pos),labels=if(first)group_labels else rep('',9),expand=c(0,0))+
  theme(axis.line.y=element_blank(),axis.ticks.y=element_blank(),panel.grid=element_blank(),panel.grid.major.y=element_blank(),axis.text.y=element_text(size=6.2))
}

rain_b<-raincloud('category');rain_c<-raincloud('taxon')
plots<-list(change_plot('Q_mean'),change_plot('Q_5pct'),change_plot('R_low'),
 endpoints('category','Q_mean',TRUE),endpoints('category','Q_5pct'),rain_b,
 endpoints('taxon','Q_mean',TRUE),endpoints('taxon','Q_5pct'),rain_c,
 heat_plot('Q_mean',TRUE),heat_plot('Q_5pct'),heat_plot('R_low'),
 stacked('Q_mean',TRUE),stacked('Q_5pct'),stacked('R_low'))
ids<-as.vector(t(outer(letters[1:5],c('mean','q5','low'),paste,sep='-')))
width_mm<-183;height_mm<-225;cell_x<-c(16,72,128);cell_w<-54
row_bottom<-c(186,143,103,62,8);row_h<-c(25,35,32,32,41)
prefix<-file.path(out,'Fig_2')
grobs<-lapply(seq_along(plots),function(i){
 g<-ggplotGrob(plots[[i]]);p<-g$layout[g$layout$name=='panel',][1,]
 left<-seq_len(p$l-1);right<-seq.int(p$r+1,length(g$widths));
 # Padding allocations are physical and independent of label string length.
 g$widths[left]<-unit(rep(0,length(left)),'mm');g$widths[p$l-1]<-unit(6,'mm')
 g$widths[right]<-unit(rep(0,length(right)),'mm');g$widths[p$r+1]<-unit(4,'mm')
 # All plot labels are permitted to use the unused outer gutter.
 g$layout$clip[g$layout$name!='panel']<-'off'
 black_text_grob(g)
})
draw_figure<-function(measure=FALSE){
 grid.newpage();grid.rect(gp=gpar(fill='white',col=NA));measured<-list()
 txt<-function(s,x,y,fs=7,col=ink,bold=FALSE,just='left')grid.text(s,unit(x,'mm'),unit(y,'mm'),just=just,gp=gpar(fontfamily='Helvetica',fontsize=fs,col=black_text_colour(col),fontface=if(bold)'bold' else 'plain'))
 for(j in 1:3)txt(c('Mean habitat quality','5th-percentile habitat quality','Proportion of low-quality habitat (%)')[j],cell_x[j]+28,221,8,bold=TRUE,just='center')
 titles<-c('Change from 2000','IUCN Red List Categories','Taxonomic groups','Category-by-taxon means','Direction of change: 2020 versus 2000');yy<-c(215,183,140,96,56)
 for(i in 1:5){txt(letters[i],3,yy[i],10,bold=TRUE);txt(titles[i],8,yy[i],7.6,bold=TRUE)}
 grid.segments(unit(112,'mm'),unit(yy[1],'mm'),unit(115,'mm'),unit(yy[1],'mm'),gp=gpar(col=blue,lwd=1.4));txt('Species',116,yy[1],6.2)
 grid.segments(unit(145,'mm'),unit(yy[1],'mm'),unit(148,'mm'),unit(yy[1],'mm'),gp=gpar(col=gray,lwd=1.2,lty='22'));txt('National mean',149,yy[1],6.2)

 grid.points(unit(c(81,102),'mm'),unit(rep(yy[2],2),'mm'),pch=c(21,21),size=unit(1.7,'mm'),gp=gpar(col=ink,fill=c('white',ink)))
 txt('2000',84,yy[2],6.3);txt('2020',105,yy[2],6.3)
 txt('Mean +/- SEM',82,yy[3],6.2);txt('Change (%)',115,yy[3],6.2)
 for(k in 1:3){xx<-87+(k-1)*30
  grid.rect(unit(xx+1.2,'mm'),unit(56,'mm'),width=unit(2.3,'mm'),height=unit(1.8,'mm'),gp=gpar(col=NA,fill=direction_col[k]))
  txt(names(direction_col)[k],xx+4.2,56,6.5,col='#47535B')}
 txt('Group  n',3,49,6.3,col='#657078')
 txt('Raincloud | log1p axis',137,yy[2],6.1);txt('Raincloud | log1p axis',137,yy[3],6.1)
 for(i in seq_along(grobs)){
  r<-(i-1)%/%3+1;c<-(i-1)%%3+1
  pushViewport(viewport(x=unit(cell_x[c],'mm'),y=unit(row_bottom[r],'mm'),width=unit(cell_w,'mm'),height=unit(row_h[r],'mm'),just=c('left','bottom'),name=paste0('cell-',ids[i])))
  grid.draw(grobs[[i]])
  if(measure){grid.force();vs<-grid.ls(viewports=TRUE,grobs=FALSE,recursive=TRUE,print=FALSE);wh<-which(grepl('^panel[.]',vs$name)&grepl(paste0('cell-',ids[i]),vs$vpPath,fixed=TRUE))
   stopifnot(length(wh)>=1);k<-tail(wh,1);full<-paste(vs$vpPath[k],vs$name[k],sep='::');upViewport(0);parts<-strsplit(full,'::',fixed=TRUE)[[1]];downViewport(do.call(vpPath,as.list(parts[parts!='ROOT'])))
   lo<-deviceLoc(unit(0,'npc'),unit(0,'npc'),valueOnly=TRUE);hi<-deviceLoc(unit(1,'npc'),unit(1,'npc'),valueOnly=TRUE)
   measured[[i]]<-list(id=ids[i],bbox_pt=c(lo$x,lo$y,hi$x,hi$y)*72,grid_id='figure4',row_start=r-1L,row_stop=r,col_start=c-1L,col_stop=c);upViewport(0)
  }else upViewport()
 }
 if(measure)return(measured)
}
pdf(paste0(prefix,'.pdf'),width=width_mm/25.4,height=height_mm/25.4,family='Helvetica',useDingbats=FALSE);draw_figure();dev.off()
svglite::svglite(paste0(prefix,'.svg'),width=width_mm/25.4,height=height_mm/25.4);draw_figure();dev.off()
ragg::agg_png(paste0(prefix,'.png'),width=width_mm,height=height_mm,units='mm',res=600,background='white');draw_figure();dev.off()
