# Change text only; preserve all data-mark and explicitly chromatic text colours.
neutral_text_colours <- c("#2D3D46", "#344047", "#47535B", "#657078")
black_text_colour <- function(col) {
  if(is.null(col)) return(col)
  hit <- !is.na(col) & toupper(substr(col,1,7)) %in% neutral_text_colours
  col[hit] <- ifelse(nchar(col[hit])==9,paste0("#000000",substr(col[hit],8,9)),"#000000")
  col
}
black_text_grob <- function(g) {
  if(inherits(g,"text") && !is.null(g$gp$col)) g$gp$col <- black_text_colour(g$gp$col)
  if(!is.null(g$grobs)) g$grobs <- lapply(g$grobs,black_text_grob)
  if(!is.null(g$children)) for(k in seq_along(g$children)) g$children[[k]] <- black_text_grob(g$children[[k]])
  g
}
