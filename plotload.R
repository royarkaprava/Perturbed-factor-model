library(IMIFA)

plotload <- function(x, na.col = "#808080FF", ptype = c("image", "points"), 
                     border.col = "#808080FF", dlabels = NULL, rlabels = FALSE, 
                     clabels = FALSE, pch = 15, cex = 3, label.cex = 0.6){
  cmat <- mat2cols(t(x))
  N <- nrow(cmat)
  P <- ncol(cmat)
  cmat <- replace(cmat, is.na(cmat), na.col)
  levels <- sort(unique(as.vector(cmat)))
  z <- matrix(unclass(factor(cmat, levels = levels, labels = seq_along(levels))), 
              nrow = N, ncol = P)
  info <- list(x = seq_len(P), y = seq_len(N), z = t(z), 
               col = levels)
  return(info)
}

Infolam <- plotload(lambda)
graphics::image(Infolam$x, Infolam$y, Infolam$z, col = Infolam$col, main = "", xlab="", ylab="")
fields::image.plot(Infolam$x, Infolam$y, Infolam$z, col = Infolam$col, main = "", xlab="", ylab="")
