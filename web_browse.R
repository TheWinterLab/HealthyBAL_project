web_browse <- function(x, ...) {
  UseMethod("web_browse")
}

web_browse.data.frame <- function(x, ...) {
  file = tempfile(paste0(deparse(substitute(x)), "_"), fileext = ".csv")
  readr::write_csv(x, path = file)
  browseURL(file)
  
}

web_browse.gg <- function(x,
                          width = 20,
                          height = 20,
                          save = FALSE,
                          file_name = NULL,
                          ...) {
  if (isTRUE(save)) {
    if (is.null(file_name)) {
      file = paste0(gsub("-", "", Sys.Date()), "_image_output.png")
    } else {
      file = file_name
    }
    
    ggplot2::ggsave(
      plot = x,
      filename = file,
      width = width,
      height = height,
      
    )
  } else {
    file = tempfile(paste0(deparse(substitute(x)), "_"), fileext = ".png")
    ggplot2::ggsave(
      plot = x,
      filename = file,
      width = width,
      height = height
    )
  }
  
  browseURL(file)
}


web_browse.pheatmap <- function(x,
                                width = 20,
                                height = 20,
                                save = FALSE,
                                file_name = NULL,
                                ...) {
  file = tempfile(paste0(deparse(substitute(x)), "_"),
                  fileext = ".png")
  
  if (isTRUE(save)) {
    if (is.null(file_name)) {
      file = paste0(gsub("-", "", Sys.Date()), "_image_output.png")
    } else {
      file = file_name
    }
  }
  
  png(
    file,
    width = width,
    height = height,
    units = "in",
    res = 300
  )
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
  browseURL(file)
}
