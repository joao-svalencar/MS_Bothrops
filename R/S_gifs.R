library(raster)
library(magick)

rasters <- list(pres_consenso,
                fut1.1_consenso,
                fut1.2_consenso,
                fut2.1_consenso,
                fut2.2_consenso)

rasters_bin <- list(binario,
                binario_fut1.1,
                binario_fut1.2,
                binario_fut2.1,
                binario_fut2.2)

# Titles for each frame
titles <- c("Present", "Futuro 1.1", "Futuro 1.2", "Futuro 2.1", "Futuro 2.2")

imgs <- mapply(function(r, title) {
  tmp <- tempfile(fileext = ".png")
  png(tmp, width = 800, height = 600)
  
  # Optional styling
  par(cex.main = 1.5, font.main = 2)
  
  # Each plot gets its own scale automatically
  plot(r, main = title)
  
  dev.off()
  image_read(tmp)
}, rasters, titles, SIMPLIFY = FALSE) #trocar para rasters_bin para binarios


gif <- image_animate(image_join(imgs), fps = 1)
image_write(gif, here::here("outputs","gifs","fonsecai.gif"))
