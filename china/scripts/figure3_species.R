# Upper Figure 5 module: species heterogeneity (a) and paired group means (b).
# All observations and uncertainty summaries are retained from the verified data.
suppressPackageStartupMessages({library(ggplot2);library(dplyr);library(grid);library(jsonlite)})
.species_sections <- local({
  d <- read.csv(file.path(data,'Figure_S11_paired_species.csv'))
  ordering <- read.csv(file.path(data,'Fig_3_species_order.csv'))
  d <- merge(d,ordering,by='Species_ID')
  d <- d[order(d$x),]
  d$taxon <- d$Taxon
  d$ANCI_direction <- c(Increase='Increasing',Decrease='Decreasing',Unchanged='Unchanged')[d$direction]
  d$APCI_direction <- c(Increase="Increasing",Decrease="Decreasing",Unchanged="Unchanged")[d$APCI_direction]
  est <- read.csv(file.path(data,'Fig_3_group_mean_sem.csv'))
  metrics <- c('ANCI','APCI'); cats <- c('CR','EN','VU','NT','LC')
  taxa <- c('Amphibians','Reptiles','Mammals','Plants')
  tax_col <- c(Amphibians='#467FAD',Reptiles='#C89449',Mammals='#B45B68',Plants='#59947A')
  cat_col <- c(CR='#A84343',EN='#CF7146',VU='#C7A54A',NT='#91A75A',LC='#4F957B')
  direction_col <- c(Increasing='#3B75AF',Decreasing='#C95849',Unchanged='#C5C9CC')
  ink <- '#2D3D46'; faint <- '#E7ECEF'
  stopifnot(nrow(d)==1820L,!anyDuplicated(d$Species_ID),identical(as.integer(d$x),seq_len(1820L)),nrow(est)==80L)
  stopifnot(all(is.finite(as.matrix(d[,c('ANCI_2000','ANCI_2020','APCI_2000','APCI_2020')]))),all(is.finite(est$sem)))
  counts <- setNames(lapply(metrics,function(m)table(factor(d[[paste0(m,'_direction')]],levels=names(direction_col)))),metrics)
  stopifnot(identical(as.integer(counts$ANCI),c(1407L,371L,42L)),identical(as.integer(counts$APCI),c(1289L,529L,2L)))
  stopifnot(max(d$ANCI_2020,d$APCI_2020)<=80,all(est$mean-est$sem>=0),all(est$mean+est$sem<=30))
  for(m in metrics)for(y in c(2000,2020))for(t in taxa)for(ca in cats){
    values <- d[d$taxon==t & d$IUCN_clean==ca,paste0(m,'_',y)]
    e <- est[est$metric==m & est$year==y & est$taxon==t & est$IUCN_clean==ca,]
    stopifnot(nrow(e)==1L,e$n==length(values),abs(e$mean-mean(values))<1e-10,abs(e$sem-sd(values)/sqrt(length(values)))<1e-10)
  }
  segments <- d %>% group_by(taxon,IUCN_clean) %>% summarise(n=n(),start=min(x)-.5,end=max(x)+.5,.groups='drop')
  tax_segments <- d %>% group_by(taxon) %>% summarise(n=n(),start=min(x)-.5,end=max(x)+.5,.groups='drop')
  tax_segments <- tax_segments[match(taxa,tax_segments$taxon),]
  boundaries <- head(tax_segments$end,-1)
  theme_local <- theme_classic(base_size=7,base_family='Helvetica') + theme(
    axis.title=element_blank(),axis.text=element_text(size=6.5,color=ink),
    axis.line=element_line(linewidth=.23,color='#A8B2B9'),axis.ticks=element_line(linewidth=.23,color='#A8B2B9'),
    axis.ticks.length=unit(.8,'mm'),legend.position='none',panel.grid=element_blank(),
    plot.margin=margin(0,0,0,0),plot.background=element_blank())
  track <- function(m){
    dd <- d;dd$value <- dd[[paste0(m,'_2020')]];dd$direction_display <- dd[[paste0(m,'_direction')]]
    ggplot(dd,aes(x,value)) + geom_hline(yintercept=seq(20,80,20),color=faint,linewidth=.18) +
      geom_vline(xintercept=boundaries,color='#B4BEC4',linewidth=.25,linetype='22') +
      geom_segment(aes(xend=x,y=0,yend=value,color=direction_display),linewidth=.13,alpha=.7) +
      geom_point(aes(color=direction_display),size=.20,stroke=0) + scale_color_manual(values=direction_col) +
      scale_x_continuous(limits=c(.5,1820.5),expand=c(0,0)) + scale_y_continuous(limits=c(0,80),breaks=seq(0,80,20),expand=c(0,0)) +
      theme_local + theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())
  }
  dumbbell <- function(m,t){
    dd <- est[est$metric==m & est$taxon==t,];p0 <- dd[dd$year==2000,];p1 <- dd[dd$year==2020,]
    pp <- merge(p0,p1,by=c('taxon','IUCN_clean','y'),suffixes=c('0','1'));pp <- pp[match(cats,pp$IUCN_clean),]
    labels <- if(m=='ANCI')paste0(cats,' (',pp$n0,')') else cats
    ggplot() + geom_vline(xintercept=c(10,20,30),linewidth=.18,color=faint) +
      geom_segment(data=pp,aes(x=mean0,xend=mean1,y=yp0,yend=yp1,color=IUCN_clean),linewidth=.8,alpha=.65) +
      geom_segment(data=dd,aes(x=mean-sem,xend=mean+sem,y=yp,yend=yp,color=IUCN_clean),linewidth=.35) +
      geom_point(data=p0,aes(mean,yp,color=IUCN_clean),shape=21,fill='white',size=1.9,stroke=.5) +
      geom_point(data=p1,aes(mean,yp,color=IUCN_clean,fill=IUCN_clean),shape=21,size=1.9,stroke=.5) +
      scale_color_manual(values=cat_col) + scale_fill_manual(values=cat_col) +
      scale_x_continuous(limits=c(0,30),breaks=c(0,10,20,30),expand=c(0,0)) +
      scale_y_continuous(limits=c(.48,5.55),breaks=5:1,labels=labels,expand=c(0,0)) +
      theme_local + theme(axis.line.y=element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_text(size=6.5))
  }
  plots <- c(lapply(metrics,track),lapply(taxa,function(t)dumbbell('ANCI',t)),lapply(taxa,function(t)dumbbell('APCI',t)))
  ids <- c('track-ANCI','track-APCI',paste0('ANCI-',tolower(taxa)),paste0('APCI-',tolower(taxa)))
  cell_x <- c(9,9,rep(c(9,52,95,138),2));cell_y <- c(103,78,rep(32,4),rep(3,4))
  cell_w <- c(171,171,rep(42,8));cell_h <- c(15,15,rep(24,8))
  pdf(NULL,width=183/25.4,height=134/25.4,family='Helvetica',useDingbats=FALSE)
  grobs <- lapply(seq_along(plots),function(i){
    g <- ggplotGrob(plots[[i]]);p <- g$layout[g$layout$name=='panel',][1,]
    ll <- seq_len(p$l-1);rr <- seq.int(p$r+1,length(g$widths))
    g$widths[ll] <- unit(rep(0,length(ll)),'mm');g$widths[rr] <- unit(rep(0,length(rr)),'mm')
    g$widths[p$l-1] <- unit(if(i<=2)6 else 11,'mm');g$widths[p$r+1] <- unit(1,'mm')
    if(i<=2){g$heights <- unit(rep(0,length(g$heights)),'mm');g$heights[p$t] <- unit(15,'mm')}
    g$layout$clip[g$layout$name!='panel'] <- 'off';black_text_grob(g)
  })
  invisible(dev.off())
  draw <- function(measure=FALSE,newpage=FALSE){
    if(newpage)grid.newpage()
    base_path <- current.vpPath()
    grid.rect(gp=gpar(fill='white',col=NA));measured <- list()
    txt <- function(s,x,y,fs=7,col=ink,bold=FALSE,just='left')grid.text(s,unit(x,'mm'),unit(y,'mm'),just=just,gp=gpar(fontfamily='Helvetica',fontsize=fs,col=black_text_colour(col),fontface=if(bold)'bold' else 'plain'))
    txt('a',3,130,10,col='black',bold=TRUE)
    txt('2020',3,122.8,7,bold=TRUE)
    for(k in 1:5){xx <- 132+(k-1)*9.5;grid.rect(unit(xx,'mm'),unit(130,'mm'),width=unit(2,'mm'),height=unit(1.5,'mm'),gp=gpar(fill=cat_col[k],col=NA));txt(cats[k],xx+1.7,130,6.5)}
    xpos <- function(v)15+(v-.5)/1820*164
    for(i in seq_len(nrow(segments))){s <- segments[i,];grid.rect(unit((xpos(s$start)+xpos(s$end))/2,'mm'),unit(120.5,'mm'),width=unit(xpos(s$end)-xpos(s$start),'mm'),height=unit(1.5,'mm'),gp=gpar(fill=cat_col[s$IUCN_clean],col='white',lwd=.2))}
    for(i in seq_len(nrow(tax_segments))){s <- tax_segments[i,];txt(sprintf('%s (n = %s)',s$taxon,s$n),xpos((s$start+s$end)/2),124.6,6.7,col=tax_col[s$taxon],bold=TRUE,just='center')}
    for(j in 1:2){
      grid.text(metrics[j],unit(6.7,'mm'),unit(c(110.5,85.5)[j],'mm'),rot=90,gp=gpar(fontfamily='Helvetica',fontsize=7,col=black_text_colour(ink)))
      for(k in 1:3){xx <- 18+(k-1)*55;yy <- c(99,74)[j];lab <- sprintf('%s  %s (%.1f%%)',names(direction_col)[k],format(counts[[j]][k],big.mark=',',trim=TRUE),100*counts[[j]][k]/nrow(d))
        grid.rect(unit(xx,'mm'),unit(yy,'mm'),width=unit(2.7,'mm'),height=unit(1.8,'mm'),gp=gpar(fill=direction_col[k],col=NA));txt(lab,xx+3,yy,6.5)
      }
    }
    txt('b',3,65,10,col='black',bold=TRUE)
    txt('Mean \u00b1 SEM',98,65,6.5)
    grid.points(unit(c(137,159),'mm'),unit(c(65,65),'mm'),pch=21,size=unit(1.7,'mm'),gp=gpar(col=ink,fill=c('white',ink)))
    txt('2000',140,65,6.5);txt('2020',162,65,6.5)
    for(k in 1:4)txt(taxa[k],cell_x[k+2]+26,59,7,col=tax_col[k],bold=TRUE,just='center')
    for(j in 1:2)grid.text(metrics[j],unit(4.3,'mm'),unit(c(45.2,16.2)[j],'mm'),rot=90,gp=gpar(fontfamily='Helvetica',fontsize=7,col=black_text_colour(ink)))
    for(i in seq_along(grobs)){
      pushViewport(viewport(x=unit(cell_x[i],'mm'),y=unit(cell_y[i],'mm'),width=unit(cell_w[i],'mm'),height=unit(cell_h[i],'mm'),just=c('left','bottom'),name=paste0('species-cell-',ids[i])))
      grid.draw(grobs[[i]])
      if(measure){
        grid.force();vs <- grid.ls(viewports=TRUE,grobs=FALSE,recursive=TRUE,print=FALSE)
        wh <- which(grepl('^panel[.]',vs$name)&grepl(paste0('species-cell-',ids[i]),vs$vpPath,fixed=TRUE));stopifnot(length(wh)>=1)
        k <- tail(wh,1);full <- paste(vs$vpPath[k],vs$name[k],sep='::');upViewport(0)
        parts <- strsplit(full,'::',fixed=TRUE)[[1]];downViewport(do.call(vpPath,as.list(parts[parts!='ROOT'])))
        lo <- deviceLoc(unit(0,'npc'),unit(0,'npc'),valueOnly=TRUE);hi <- deviceLoc(unit(1,'npc'),unit(1,'npc'),valueOnly=TRUE)
        measured[[i]] <- list(id=ids[i],bbox_pt=c(lo$x,lo$y,hi$x,hi$y)*72,grid_id=if(i<=2)'tracks' else 'summaries',
          row_start=if(i<=2)i-1L else (i-3L)%/%4L,row_stop=if(i<=2)i else (i-3L)%/%4L+1L,
          col_start=if(i<=2)0L else (i-3L)%%4L,col_stop=if(i<=2)1L else (i-3L)%%4L+1L)
        upViewport(0);if(!is.null(base_path))downViewport(base_path)
      }else upViewport()
    }
    if(measure)return(measured)
    invisible(NULL)
  }
  list(draw=draw)
})
draw_species_sections <- function(measure=FALSE,newpage=FALSE).species_sections$draw(measure,newpage)

