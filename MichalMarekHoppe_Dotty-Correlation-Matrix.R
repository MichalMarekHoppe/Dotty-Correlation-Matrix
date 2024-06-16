# Dotty Correlation Matrix
# by Dr Michal Marek Hoppe

library(ggplot2)

# load data in WIDE format
data <- read.csv("data.CSV", row.names = 1)

## Dotty Correlation Matrix ----
dotty_corr_matrix <- 
  function(data,
           total_n = TRUE,
           max_colour = "#344470",
           min_colour = "#CA6D88",
           method = "spearman",
           title = "Dotty Correlation Matrix",
           box_p = 0.01,
           box_r = 0.5,
           params = 
             list("widths" = list("cell" = 10,
                                  "minor" = 2,
                                  "major" = 12,
                                  "margin" = 4),
                  "label_dist" = 2,
                  "label_size" = 2,
                  "label_colour" = "grey33",
                  "box_fill" = "grey94",
                  "box_nan_fill" = "grey98",
                  "box_X_colour" = "grey75",
                  "sig_colour" = "#CA6D88",
                  "circle_scale" = 0.95)) {
    # functions
    # gradient palette
    colfunc <- colorRampPalette(c(min_colour, "grey75", max_colour))
    params$colours = colfunc(201)
    names(params$colours) <-
      gsub(" ", "", as.character(format(seq(-1, 1, by = 0.01), digits = 2)))
    
    # matrix size calculator
    corr_size_calculator <- function (x) {
      width <-
        sum(sum(!grepl(paste(c("minor",
                               "major"),
                             collapse = "|"),
                       colnames(x))) * params$widths$cell,
            sum(grepl("minor", colnames(x))) * params$widths$minor,
            sum(grepl("minor", colnames(x))) * params$widths$major)
      return(width)
    }
    
    # cell locator
    corr_loc <- function(i, forward = T) {
      nuggets <- colnames(data)[1:i]
      if (forward) {
        width <- 0
        for (n in seq_along(nuggets)) {
          width <-
            width + ifelse(grepl("minor", nuggets[n]),
                           params$widths$minor,
                           ifelse(grepl("major", nuggets[n]),
                                  params$widths$major,
                                  params$widths$cell))
        }
      } else {
        width <- corr_size_calculator(data)
        for (n in seq_along(nuggets)) {
          width <- width - ifelse(grepl("minor", nuggets[n]),
                                  params$widths$minor,
                                  ifelse(grepl("major", nuggets[n]),
                                         params$widths$major,
                                         params$widths$cell))
        }
      }
      width <- 
        ifelse(forward,
               width - 0.5 * params$widths$cell,
               width + 0.5 * params$widths$cell)
      return(width)
    }
    
    # generate circle points
    gen_cir <- function(center_x, center_y, radius, n_points = 100) {
      theta <- seq(0, 2 * pi, length.out = n_points)
      x <- center_x + radius * cos(theta)
      y <- center_y + radius * sin(theta)
      return(data.frame(x = x, y = y))
    }
    
    # generate matrix
    if (method == "spearman") {
      corr <- suppressWarnings(Hmisc::rcorr(as.matrix(data), type = "spearman")) 
    } else if (method == "pearson") {
      corr <- suppressWarnings(Hmisc::rcorr(as.matrix(data), type = "pearson")) 
    }
    
    corr$r[upper.tri(corr$r, diag = TRUE)] <- NA
    corr$p_adj <-
      ifelse(corr[["P"]]*sum(!is.na(corr[["r"]])) < 1, # account only for unqique comparisons
             corr[["P"]]*sum(!is.na(corr[["r"]])), 1)
    
    
    
    # generate plot
    p <- ggplot() + 
      scale_x_continuous(limits = c(-params$widths$margin * params$widths$cell,
                                    corr_size_calculator(data)),
                         expand = c(0.01, 0.01)) +
      scale_y_continuous(limits = c(-params$widths$margin * params$widths$cell,
                                    corr_size_calculator(data)),
                         expand = c(0.01, 0.01)) +
      theme(plot.title = element_text(size = 8, hjust = 0.5,
                                      margin = margin(b = 2)),
            plot.subtitle = element_text(size = 6, hjust = 0.5,
                                         colour = "grey66",
                                         margin = margin(b = 1)),
            panel.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank()) +
      labs(title = title,
           subtitle = paste0("as of ", format(Sys.time(), format = "%Y-%m-%d")))
    
    
    # add lables
    for (i in seq_along(colnames(data))) {
      if (!grepl(paste(c("minor", "major"), collapse = "|"), 
                 colnames(data)[i])) {
        p[["layers"]][[length(p[["layers"]])+1]] <- 
          geom_text(aes(x = !!corr_loc(i),
                        y = !!(0 - params$label_dist)),
                    label = colnames(data)[i],
                    vjust = 0.5, hjust = 1, angle = 90,
                    color = params$label_colour,
                    size = params$label_size)
        p[["layers"]][[length(p[["layers"]])+1]] <- 
          geom_text(aes(x = !!(0 - params$label_dist),
                        y = !!corr_loc(i, F)),
                    label = colnames(data)[i],
                    vjust = 0.5, hjust = 1, angle = 0,
                    color = params$label_colour,
                    size = params$label_size)
      }
    }
    # add data
    for (x in seq_along(colnames(data))) {
      for (y in seq_along(colnames(data))) {
        value <- round(corr$r[y, x],2)
        if (is.nan(value)) {
          p[["layers"]][[length(p[["layers"]])+1]] <- 
            geom_rect(aes(xmin = !!(corr_loc(x) - 0.5 * params$widths$cell),
                          xmax = !!(corr_loc(x) + 0.5 * params$widths$cell),
                          ymin = !!(corr_loc(y, F) - 0.5 * params$widths$cell), 
                          ymax = !!(corr_loc(y, F) + 0.5 * params$widths$cell)), 
                      fill = params$box_nan_fill, color = "white")
          p[["layers"]][[length(p[["layers"]])+1]] <- 
            geom_text(aes(x = !!corr_loc(x),
                          y = !!corr_loc(y, F)),
                      label = "X",
                      vjust = 0.5, hjust = 0.5,
                      color = params$box_X_colour,
                      size = params$label_size * 0.75)
        } else {
          if (!is.na(value)) {
            # background square
            p[["layers"]][[length(p[["layers"]])+1]] <- 
              geom_rect(aes(xmin = !!(corr_loc(x) - 0.5 * params$widths$cell),
                            xmax = !!(corr_loc(x) + 0.5 * params$widths$cell),
                            ymin = !!(corr_loc(y, F) - 0.5 * params$widths$cell), 
                            ymax = !!(corr_loc(y, F) + 0.5 * params$widths$cell)), 
                        fill = params$box_fill, color = "white")
          }
        }
      }
    }
    # add significance square
    for (x in seq_along(colnames(data))) {
      for (y in seq_along(colnames(data))) {
        value <- round(corr$r[y, x],2)
        if (!is.na(value)) {
          if (corr$p_adj[y, x] <= box_p & abs(corr$r[y, x]) >= box_r ) { # define significance level
            p[["layers"]][[length(p[["layers"]])+1]] <- 
              geom_rect(aes(xmin = !!(corr_loc(x) - 0.5 * params$widths$cell),
                            xmax = !!(corr_loc(x) + 0.5 * params$widths$cell),
                            ymin = !!(corr_loc(y, F) - 0.5 * params$widths$cell), 
                            ymax = !!(corr_loc(y, F) + 0.5 * params$widths$cell)), 
                        fill = NA, color = params$sig_colour,
                        linewidth = 0.25)
          }
        }
      }
    }
    # add circles
    for (x in seq_along(colnames(data))) {
      for (y in seq_along(colnames(data))) {
        value <- round(corr$r[y, x],2)
        if (!is.na(value)) {
          circle_data <- 
            gen_cir("center_x" = corr_loc(x),
                    "center_y" = corr_loc(y, F),
                    "radius" = ((params$widths$cell/2) * params$circle_scale) * abs(value))
          
          p[["layers"]][[length(p[["layers"]])+1]] <- 
            annotate("polygon", x = circle_data$x, y = circle_data$y,
                     fill = params$colours[names(params$colours) == format(value, 
                                                                           nsmall = 2)])
        }
      }
    }
    if (total_n == TRUE) {
      p[["layers"]][[length(p[["layers"]])+1]] <- 
        geom_text(aes(x = !!(corr_size_calculator(data) * 0.66),
                      y = !!(corr_size_calculator(data) * 0.66)),
                  label = paste0("n = ", nrow(data)),
                  vjust = 0.5, hjust = 0,
                  color = params$label_colour,
                  size = params$label_size)
    }
    return(list(graph = p,
                matrix = corr))
  }

## print graphs ----
print_graphs <- function (p, width, height, res, name, folder = "") {
  folder <- ifelse(folder == "",
                   "",
                   paste0(folder, "/"))
  
  if (!dir.exists(folder)) {
    dir.create(folder, recursive = TRUE, showWarnings = FALSE)
  }
  ggsave(paste0(folder, name, ".pdf"),
         plot = p,
         width = width,
         height = height,
         dpi = res,
         units = "px",
         bg = "white")
  png(filename = paste0(folder, name, ".png"),
      width = width,
      height = height,
      res = res,
      units = "px",
      bg = "white")
  options(warn = -1)
  print(p)
  invisible(dev.off())
  options(warn = 0)
  return(p)
}

## generate Dotty object ----
dotty <-
  dotty_corr_matrix(data, # data in WIDE format
                    total_n = TRUE, # show total number of comparisons
                    max_colour = "#344470", # colour for r = 1
                    min_colour = "#CA6D88", # colour for r = -1
                    method = "spearman", # correlation method: spearman/pearson
                    title = "Dotty Correlation Matrix", # plot title
                    box_p = 0.001, # max adjusted p-value for significance box 
                    box_r = 0.5, # minimun abs(r) for significance box
                    params = 
                      list("widths" = list("cell" = 10, # cell size
                                           "minor" = 2, # minor divider
                                           "major" = 12, # major divider
                                           "margin" = 4), # margin size
                           "label_dist" = 2, # distance of labels from matrix
                           "label_size" = 2, # label size
                           "label_colour" = "grey33", # label text colour 
                           "box_fill" = "grey94",# box fill colour
                           "box_nan_fill" = "grey98", # box fill of missing correlations
                           "box_X_colour" = "grey75", # X colour for missing correlations
                           "sig_colour" = "#CA6D88", # significance box colour
                           "circle_scale" = 0.95)) # Dotty circle size relative to box size

## print Dotty Correlation Matrix and save matrix RDS----
print_graphs(dotty$graph,
             1200, # width
             1200, # height
             300, # resolution
             "dotty_corr_matrix" # file name
); saveRDS(dotty$matrix, "dotty.RDS")

