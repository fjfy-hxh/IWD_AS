#' Filter a Seurat object based on a CSV file and add metadata
#'
#' @param csv_file Path to a CSV file with cell information; row names must be cell IDs
#' @param seurat_obj Input Seurat object
#' @param metadata_cols Vector of column names to add from the CSV to the Seurat object
#' 
#' @return A filtered Seurat object with metadata added
#' @export
#' By default the first column of the CSV is the cell barcode and becomes the row names
#' @examples
#' # af.ST_radius <- add_metadata_from_csv("neighbor_data.csv", af.ST_both, c("radius", "distance"))
add_metadata_from_csv <- function(csv_file, seurat_obj, metadata_cols = NULL) {
  # Read the CSV file (ensure the first column holds cell names)
  neighbor_data <- read.csv(csv_file, row.names = 1, stringsAsFactors = FALSE)
  
  # Get common cells (strict matching)
  common_cells <- intersect(rownames(neighbor_data), colnames(seurat_obj))
  
  # Check whether there are common cells
  if (length(common_cells) == 0) {
    stop("Error: no common cell names between the Seurat object and the CSV file")
  }
  
  # Filter the Seurat object (strict filtering)
  seurat_filtered <- subset(seurat_obj, cells = common_cells)
  
  # Filter the CSV data and strictly align the order
  neighbor_data <- neighbor_data[common_cells, , drop = FALSE]
  neighbor_data <- neighbor_data[colnames(seurat_filtered), , drop = FALSE] # Confirm the order again
  
  # Validate the cell order
  if (!identical(rownames(neighbor_data), colnames(seurat_filtered))) {
    stop("Final cell order matching failed")
  }
  
  # Determine the columns to add
  if (is.null(metadata_cols)) {
    metadata_cols <- colnames(neighbor_data) # Add all columns by default
  } else {
    metadata_cols <- intersect(metadata_cols, colnames(neighbor_data))
  }
  
  # Add metadata
  for (col in metadata_cols) {
    seurat_filtered[[col]] <- neighbor_data[[col]]
  }
  
  # Print a processing summary
  message(
    sprintf("Successfully matched %d cells\nAdded %d metadata columns from CSV:\n%s",
            ncol(seurat_filtered),
            length(metadata_cols),
            paste(metadata_cols, collapse = ", "))
  )
  
  return(seurat_filtered)
}


#' Compute and classify the minimum distance from each cell to the specified target cells
#'
#' @param distance_data Input data frame; must contain columns: group, Celltype, x_avg, y_avg
#' @param target_celltype Target cell type (default "FB C4")
#' @param distance_breaks Vector of distance breakpoints (default c(-Inf, 0, 50, 100, 200, Inf))
#' @param distance_labels Distance labels (generated automatically by default)
#'
#' @return A data frame with min_distance and radius_range columns added
#'
#' @examples
#' result <- calculate_min_distance_to_target(
#'   distance_data,
#'   target_celltype = "FB C4",
#'   distance_breaks = c(-Inf, 0, 10, 50, 100, Inf)
#' )
calculate_min_distance_to_target <- function(
    distance_data,
    target_celltype = "FB C4",
    distance_breaks = c(-Inf, 0, 50, 100, 200, Inf),
    distance_labels = NULL
) {
  # Parameter validation
  required_cols <- c("group", "Celltype", "x_avg", "y_avg")
  if (!all(required_cols %in% colnames(distance_data))) {
    stop("Input data frame must contain columns: ", paste(required_cols, collapse = ", "))
  }
  
  # Auto-generate labels (if not provided)
  if (is.null(distance_labels)) {
    distance_labels <- sapply(1:(length(distance_breaks)-1), function(i) {
      if (i == 1) {
        if (distance_breaks[i] == -Inf) {
          return(paste0("0-", distance_breaks[i+1]))
        } else {
          return(paste0("<=", distance_breaks[i+1]))
        }
      } else if (i == length(distance_breaks)-1) {
        if (distance_breaks[i+1] == Inf) {
          return(paste0(">", distance_breaks[i]))
        } else {
          return(paste0(distance_breaks[i], "-", distance_breaks[i+1]))
        }
      } else {
        return(paste0(distance_breaks[i], "-", distance_breaks[i+1]))
      }
    })
  }
  
  # Initialize the result
  result_df <- data.frame()
  groups <- unique(distance_data$group)
  
  for (current_group in groups) {
    group_data <- distance_data[distance_data$group == current_group, ]
    target_cells <- group_data[group_data$Celltype == target_celltype, ]
    
    if (nrow(target_cells) == 0) {
      warning("Group ", current_group, " does not contain target cell type '", target_celltype, "', skipping")
      next
    }
    
    # Compute the minimum distance
    group_data$min_distance_to_target <- sapply(1:nrow(group_data), function(i) {
      if (group_data$Celltype[i] == target_celltype) {
        return(0)
      }
      distances <- sqrt((group_data$x_avg[i] - target_cells$x_avg)^2 + 
                          (group_data$y_avg[i] - target_cells$y_avg)^2)
      min(distances)
    })
    
    # Classify distance ranges
    group_data$radius_range <- cut(
      group_data$min_distance_to_target,
      breaks = distance_breaks,
      labels = distance_labels,
      include.lowest = TRUE
    )
    
    result_df <- rbind(result_df, group_data)
  }
  
  return(result_df)
}








###################################
library(dplyr)

calculate_grouped_euclidean_distance <- function(
    distance_data,
    source_niche = "niche_3",
    target_niches = c("niche_1", "niche_4", "niche_6"),
    distance_breaks = c(0, 50, 100, 200, Inf),
    distance_labels = c("0-50", "50-100", "100-200", ">200")
) {
  # 1. Initialize result columns
  distance_data$euclidean_distance <- NA
  distance_data$distance_range <- NA
  
  # 2. Check whether required columns exist
  required_cols <- c("cell", "group", "banksy_niche", "x_avg", "y_avg")
  missing_cols <- setdiff(required_cols, colnames(distance_data))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  # 3. Check that coordinate columns are numeric
  if (!is.numeric(distance_data$x_avg) || !is.numeric(distance_data$y_avg)) {
    stop("Coordinate columns (x_avg/y_avg) must be numeric")
  }
  
  # 4. Process by group
  groups <- unique(distance_data$group)
  for (current_group in groups) {
    group_data <- distance_data[distance_data$group == current_group, ]
    
    # Check the source region
    if (!source_niche %in% group_data$banksy_niche) {
      message("Skipping group ", current_group, ": missing source region '", source_niche, "'")
      next
    }
    
    # Check the target regions
    valid_targets <- intersect(target_niches, unique(group_data$banksy_niche))
    if (length(valid_targets) == 0) {
      message("Skipping group ", current_group, ": missing target regions")
      next
    }
    
    # Extract cells
    source_cells <- group_data[group_data$banksy_niche == source_niche, ]
    target_cells <- group_data[group_data$banksy_niche %in% valid_targets, ]
    
    # Safely compute distances (key fix)
    if (nrow(source_cells) == 0 || nrow(target_cells) == 0) {
      message("Skipping group ", current_group, ": source or target cell count is 0")
      next
    }
    
    # Compute Euclidean distances (vectorized)
    source_coords <- as.matrix(source_cells[, c("x_avg", "y_avg")])
    target_coords <- as.matrix(target_cells[, c("x_avg", "y_avg")])
    
    dist_matrix <- sqrt(
      outer(source_coords[,1], target_coords[,1], "-")^2 +
        outer(source_coords[,2], target_coords[,2], "-")^2
    )
    
    source_cells$euclidean_distance <- apply(dist_matrix, 1, min)
    
    # Classify distance ranges
    source_cells$distance_range <- cut(
      source_cells$euclidean_distance,
      breaks = distance_breaks,
      labels = distance_labels,
      include.lowest = TRUE
    )
    
    # Update the results
    distance_data[distance_data$cell %in% source_cells$cell, "euclidean_distance"] <-
      source_cells$euclidean_distance
    distance_data[distance_data$cell %in% source_cells$cell, "distance_range"] <-
      source_cells$distance_range
    
    message("Group ", current_group, " done (", nrow(source_cells), " cells)")
  }
  
  return(distance_data)
}





#' Plot a percentage stacked area chart as distance changes
#'
#' @param seurat_obj Seurat object
#' @param radius_range_col Name of the radius range column, default "radius_range"
#' @param domain_col Name of the domain column, default "domain"
#' @param x_labels X-axis labels, default c("0-20", "20-100", "100-200", "200-500")
#' @param color_palette Color palette; uses the viridis plasma palette by default
#' @param title Plot title
#' @param xlab X-axis label
#' @param ylab Y-axis label
#' @param legend_title Legend title
#' @param alpha Transparency, default 0.8
#'
#' @return A ggplot2 object
plot_domain_percentage <- function(seurat_obj,
                                   radius_range_col = "radius_range",
                                   domain_col = "domain",
                                   x_labels = c("0-20", "20-100", "100-200", "200-500"),
                                   color_palette = NULL,
                                   title = NULL,
                                   xlab = "Radius Range",
                                   ylab = "Percentage (%)",
                                   legend_title = "Domain",
                                   alpha = 0.8) {
  
  # Load required packages
  if (!require(ggplot2)) install.packages("ggplot2")
  if (!require(dplyr)) install.packages("dplyr")
  if (!require(viridis)) install.packages("viridis")
  
  library(ggplot2)
  library(dplyr)
  library(viridis)
  
  # 1. Extract the data
  plot_data <- data.frame(
    radius_range = factor(seurat_obj@meta.data[[radius_range_col]]),
    domain = factor(seurat_obj@meta.data[[domain_col]])
  ) %>% 
    na.omit()
  
  # Check the data
  if (nrow(plot_data) == 0) {
    stop("Data extraction failed, please check whether the column names are correct")
  }
  
  # 2. Count and compute the proportion of each domain within each radius_range
  stack_data <- plot_data %>%
    dplyr::count(radius_range, domain) %>%              # Count
    dplyr::group_by(radius_range) %>%                   # Group by radius_range
    dplyr::mutate(proportion = n / sum(n) * 100) %>%   # Compute percentages
    dplyr::ungroup()                                    # Ungroup
  
  # Ensure the order of x_labels matches the levels of radius_range
  plot_data$radius_range <- factor(
    plot_data$radius_range,
    levels = x_labels,
    ordered = TRUE
  )
  
  stack_data$radius_range <- factor(
    stack_data$radius_range,
    levels = x_labels,
    ordered = TRUE
  )
  
  # 3. Create the color palette
  domain_levels <- levels(stack_data$domain)
  n_domains <- length(domain_levels)
  
  if (is.null(color_palette)) {
    # Use the viridis plasma palette by default
    colors_domain <- viridis::plasma(n_domains)
    names(colors_domain) <- domain_levels
  } else if (length(color_palette) >= n_domains) {
    # Use user-provided colors
    colors_domain <- color_palette[1:n_domains]
    names(colors_domain) <- domain_levels
  } else {
    # Not enough colors, complete with viridis
    warning(paste("Not enough colors provided, completing with the viridis palette. Required", 
                  n_domains, "colors, but only", length(color_palette), "provided"))
    colors_domain <- viridis::plasma(n_domains)
    names(colors_domain) <- domain_levels
  }
  
  # 4. Draw the stacked area chart
  p <- ggplot(stack_data, aes(x = as.numeric(factor(radius_range, levels = x_labels)), 
                              y = proportion, 
                              fill = domain)) +
    geom_area(position = "stack", alpha = alpha) +
    scale_x_continuous(
      breaks = 1:length(x_labels),
      labels = x_labels
    ) +
    labs(x = xlab, y = ylab, title = title) +
    theme_classic() +
    theme(
      legend.position = "right",
      axis.title = element_text(size = 14, face = "bold"),
      axis.text = element_text(size = 12),
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 12, face = "bold"),
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
    )
  
  # Apply the manual color settings
  p <- p + scale_fill_manual(
    values = colors_domain,
    name = legend_title
  )
  
  return(p)
}



#' Plot a percentage stacked bar chart as distance changes
#'
#' @param seurat_obj Seurat object
#' @param radius_range_col Name of the radius range column, default "radius_range"
#' @param domain_col Name of the domain column, default "domain"
#' @param x_labels X-axis labels, default c("0-20", "20-100", "100-200", "200-500")
#' @param color_palette Color palette; uses the viridis plasma palette by default
#' @param title Plot title
#' @param xlab X-axis label
#' @param ylab Y-axis label
#' @param legend_title Legend title
#' @param alpha Transparency, default 0.8
#' @param width Bar width, default 0.7
#' @param show_value Whether to display percentage values, default FALSE
#' @param angle Rotation angle of X-axis labels, default 0
#'
#' @return A ggplot2 object
plot_domain_percentage_bar <- function(seurat_obj,
                                       radius_range_col = "radius_range",
                                       domain_col = "domain",
                                       x_labels = c("0-20", "20-100", "100-200", "200-500"),
                                       color_palette = NULL,
                                       title = NULL,
                                       xlab = "Radius Range",
                                       ylab = "Percentage (%)",
                                       legend_title = "Domain",
                                       alpha = 0.8,
                                       width = 0.7,
                                       show_value = FALSE,
                                       angle = 0) {
  
  # Load required packages
  if (!require(ggplot2)) install.packages("ggplot2")
  if (!require(dplyr)) install.packages("dplyr")
  if (!require(viridis)) install.packages("viridis")
  
  library(ggplot2)
  library(dplyr)
  library(viridis)
  
  # 1. Extract the data
  plot_data <- data.frame(
    radius_range = factor(seurat_obj@meta.data[[radius_range_col]]),
    domain = factor(seurat_obj@meta.data[[domain_col]])
  ) %>% 
    na.omit()
  
  # Check the data
  if (nrow(plot_data) == 0) {
    stop("Data extraction failed, please check whether the column names are correct")
  }
  
  # 2. Count and compute the proportion of each domain within each radius_range
  stack_data <- plot_data %>%
    dplyr::count(radius_range, domain) %>%              # Count
    dplyr::group_by(radius_range) %>%                   # Group by radius_range
    dplyr::mutate(proportion = n / sum(n) * 100) %>%   # Compute percentages
    dplyr::ungroup()                                    # Ungroup
  
  # Ensure the order of x_labels matches the levels of radius_range
  if (!is.null(x_labels)) {
    stack_data$radius_range <- factor(
      stack_data$radius_range,
      levels = x_labels,
      ordered = TRUE
    )
  }
  
  # 3. Create the color palette
  domain_levels <- levels(stack_data$domain)
  n_domains <- length(domain_levels)
  
  if (is.null(color_palette)) {
    # Use the viridis plasma palette by default
    colors_domain <- viridis::plasma(n_domains)
    names(colors_domain) <- domain_levels
  } else if (length(color_palette) >= n_domains) {
    # Use user-provided colors
    colors_domain <- color_palette[1:n_domains]
    names(colors_domain) <- domain_levels
  } else {
    # Not enough colors, complete with viridis
    warning(paste("Not enough colors provided, completing with the viridis palette. Required", 
                  n_domains, "colors, but only", length(color_palette), "provided"))
    colors_domain <- viridis::plasma(n_domains)
    names(colors_domain) <- domain_levels
  }
  
  # 4. Compute cumulative percentages for label positioning
  stack_data <- stack_data %>%
    arrange(radius_range, domain) %>%
    group_by(radius_range) %>%
    mutate(
      cumulative = cumsum(proportion),
      position = cumulative - proportion/2
    ) %>%
    ungroup()
  
  # 5. Draw the stacked bar chart
  p <- ggplot(stack_data, aes(x = radius_range, y = proportion, fill = domain)) +
    geom_bar(stat = "identity", position = "stack", 
             alpha = alpha, width = width) +
    labs(x = xlab, y = ylab, title = title) +
    theme_classic() +
    theme(
      legend.position = "right",
      axis.title = element_text(size = 14, face = "bold"),
      axis.text = element_text(size = 12),
      axis.text.x = element_text(angle = angle, hjust = ifelse(angle == 0, 0.5, 1)),
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 12, face = "bold"),
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
    )
  
  # Apply the manual color settings
  p <- p + scale_fill_manual(
    values = colors_domain,
    name = legend_title
  )
  
  # 6. Whether to display percentage values
  if (show_value) {
    p <- p + geom_text(
      aes(x = radius_range, y = position, 
          label = paste0(round(proportion, 1), "%")),
      size = 3,
      color = "white",
      fontface = "bold"
    )
  }
  
  return(p)
}

###################################


library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)

plot_cell_composition_by_distance <- function(seurat_obj, 
                                              distance_col = "vertical_distance", 
                                              cell_col = "Cell_Anno",
                                              bin_num = 50,          
                                              plot_type = "proportion", 
                                              colors = NULL) {
  
  meta <- as.data.frame(seurat_obj@meta.data)
  
  if (!(distance_col %in% colnames(meta)) || !(cell_col %in% colnames(meta))) {
    stop("Error: specified distance or cell annotation columns not found in meta.data!")
  }
  
  # 💡[Key fix 1: remember the original order]
  # Extract your original cell order before force-unwrapping
  orig_cell_data <- meta[[cell_col]]
  if (is.factor(orig_cell_data)) {
    # If it is already a factor, inherit its levels directly
    cell_order <- levels(orig_cell_data)
  } else {
    # Otherwise record the unique order in which they appear in your table
    cell_order <- unique(as.character(orig_cell_data))
  }
  
  # Force-unwrap to avoid errors
  raw_dist <- unlist(meta[[distance_col]], use.names = FALSE)
  raw_dist <- as.numeric(as.character(raw_dist))
  raw_cell <- unlist(meta[[cell_col]], use.names = FALSE)
  raw_cell <- as.character(raw_cell)
  
  clean_meta <- data.frame(
    Dist = raw_dist,
    Cell = raw_cell,
    stringsAsFactors = FALSE
  ) %>% filter(!is.na(Dist), !is.na(Cell))
  
  if (nrow(clean_meta) == 0) {
    stop("Error: 0 valid rows after cleaning.")
  }
  
  # Numeric-only binning
  breaks <- seq(min(clean_meta$Dist), max(clean_meta$Dist), length.out = bin_num + 1)
  bin_centers <- (breaks[-1] + breaks[-length(breaks)]) / 2
  
  message("--> [Safety fallback] Using base R table matrices for automatic zero-padding alignment...")
  
  # Base R matrix alignment
  clean_meta$Dist_Bin <- factor(
    cut(clean_meta$Dist, breaks = breaks, include.lowest = TRUE, labels = FALSE),
    levels = 1:bin_num
  )
  
  tab_matrix <- table(Dist_Bin = clean_meta$Dist_Bin, Cell = clean_meta$Cell)
  plot_df <- as.data.frame(tab_matrix, stringsAsFactors = FALSE)
  colnames(plot_df)[3] <- "count"
  
  # Compute proportions
  plot_df <- plot_df %>%
    mutate(
      Dist_Bin = as.integer(Dist_Bin),
      Bin_Center = bin_centers[Dist_Bin]
    ) %>%
    group_by(Dist_Bin) %>%
    mutate(
      total_in_bin = sum(count),
      Proportion = ifelse(total_in_bin == 0, 0, count / total_in_bin)
    ) %>%
    ungroup()
  
  if (plot_type == "proportion") {
    plot_df$Y_Value <- plot_df$Proportion
    y_label <- "Cell Type Proportion"
    title_suffix <- "(Relative Proportion)"
  } else {
    plot_df$Y_Value <- plot_df$count
    y_label <- "Cell Count / Density per Bin"
    title_suffix <- "(Absolute Density)"
  }
  
  # 💡[Key fix 2: force-reset the factor order before plotting]
  # Drop redundant labels that do not exist and lock in your original ordering!
  valid_levels <- cell_order[cell_order %in% unique(plot_df$Cell)]
  plot_df$Cell <- factor(plot_df$Cell, levels = valid_levels)
  
  # 💡[Key fix 3: absolute binding of color key-value pairs]
  if (is.null(colors)) {
    n_cells <- length(valid_levels)
    built_in_colors <- c("#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F", 
                         "#8491B4", "#91D1C2", "#DC0000", "#7E6148", "#B09C85",
                         "#8C510A", "#BF812D", "#DFC27D", "#80CDC1", "#35978F",
                         "#C51B7D", "#E9A3C9", "#4D9221", "#A1D76A", "#1B7837")
    if (n_cells <= length(built_in_colors)) {
      colors <- built_in_colors[1:n_cells]
    } else {
      colors <- scales::hue_pal()(n_cells)
    }
    # This step is critical: anchor each color precisely to its name!
    names(colors) <- valid_levels 
  } else {
    # If you supply colors yourself (e.g. a vector of 10 colors), also assign the matching name labels
    if (is.null(names(colors)) && length(colors) >= length(valid_levels)) {
      colors <- colors[1:length(valid_levels)]
      names(colors) <- valid_levels
    }
  }
  
  message("--> Running LOESS locally weighted regression fit...")
  p <- ggplot(plot_df, aes(x = Bin_Center, y = Y_Value, color = Cell, fill = Cell)) +
    geom_smooth(method = "loess", span = 0.4, se = TRUE, alpha = 0.1, size = 1.2) +
    theme_bw() +
    scale_color_manual(values = colors) +
    scale_fill_manual(values = colors) +
    labs(
      title = paste("Cell Composition Dynamics", title_suffix),
      x = "Vertical Distance",
      y = y_label,
      color = "Cell Type",
      fill = "Cell Type"
    ) +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
      axis.title = element_text(size = 11, face = "bold"),
      axis.text = element_text(size = 10, color = "black"),
      legend.text = element_text(size = 9),
      strip.background = element_blank()
    )
  
  return(p)
}
######################################
add_spatial_centroid_reduction <- function(seurat_obj, 
                                           x_col = "x_centroid", 
                                           y_col = "y_centroid", 
                                           reduction_name = "spatial_plot",
                                           key_prefix = "spatial_") {
  
  # 1. Check whether the object is empty
  if (ncol(seurat_obj) == 0) {
    stop("Error: the input Seurat object has 0 cells/points, please check the subset result.")
  }
  
  # 2. Dynamically check whether the specified coordinate columns exist in the metadata
  meta_cols <- colnames(seurat_obj@meta.data)
  if (!(x_col %in% meta_cols) || !(y_col %in% meta_cols)) {
    stop(paste0("Error: specified coordinate columns not found in meta.data!\n",
                "Specified column names: x = '", x_col, "', y = '", y_col, "'\n",
                "Available column names: ", paste(meta_cols, collapse = ", ")))
  }
  
  message(sprintf("--> Building spatial physical coordinate matrix from [%s] and [%s]...", x_col, y_col))
  
  # 3. Build the physical coordinate matrix
  spatial_coords <- cbind(
    x = as.numeric(seurat_obj@meta.data[[x_col]]),
    y = as.numeric(seurat_obj@meta.data[[y_col]])
  )
  
  # Strictly align cell/point names (barcodes)
  rownames(spatial_coords) <- colnames(seurat_obj)
  colnames(spatial_coords) <- c(paste0(key_prefix, 1), paste0(key_prefix, 2))
  
  # 4. Get the current default assay
  current_assay <- DefaultAssay(seurat_obj)
  
  # 5. Create the standard DimReduc object
  spatial_reduction <- Seurat::CreateDimReducObject(
    embeddings = spatial_coords,
    key        = key_prefix,
    assay      = current_assay
  )
  
  # 6. Write back to the original Seurat object
  seurat_obj[[reduction_name]] <- spatial_reduction
  
  message(sprintf("✅ Spatial coordinates successfully written into reduction layer: '%s'", reduction_name))
  return(seurat_obj)
}
######################
## Plot gene maps across samples
ks_plot_spatial_features_combined <- function(seurat_obj,
                                              features = c("SCGB1A1", "CAPS"), 
                                              group_id = "A9",
                                              group_col = "group",
                                              x_col = "x_centroid",
                                              y_col = "y_centroid",
                                              reduction_name = "spatial",
                                              pt_size = 0.1,
                                              limits = c(400, 1400),
                                              xlim = NULL,
                                              ylim = NULL,
                                              show_axes = FALSE, 
                                              ncol = NULL,                     
                                              colors = c('#6090C1','#ACD2E5','#FEF9B7','#FEE395','#F2724D','#D7312D')) {
  library(Seurat)
  library(ggplot2)
  library(scales)
  library(patchwork) 
  
  # 1. Safely filter samples/groups
  if (!group_col %in% colnames(seurat_obj@meta.data)) {
    stop(paste("Group column not found in metadata:", group_col))
  }
  
  plot_obj <- subset(seurat_obj, cells = rownames(seurat_obj@meta.data[seurat_obj@meta.data[[group_col]] %in% group_id, ]))
  
  # 2. Automatically inject centroid coordinates
  plot_obj <- add_spatial_centroid_reduction(
    seurat_obj     = plot_obj, 
    x_col          = x_col, 
    y_col          = y_col, 
    reduction_name = reduction_name
  )
  
  # 3. Iterate over each feature and build a plot independently
  plot_list <- lapply(features, function(gene) {
    
    if (!gene %in% rownames(plot_obj) && !gene %in% colnames(plot_obj@meta.data)) {
      warning(paste("Feature not found in the object, skipping:", gene))
      return(NULL)
    }
    
    p <- FeaturePlot(
      object = plot_obj, 
      features = gene, 
      reduction = reduction_name,
      pt.size = pt_size
    )
    
    # 4. Overlay the custom gradient colors and clipping
    p <- p + scale_color_gradientn(
      colours = colorRampPalette(colors)(100), 
      limits  = limits,           
      oob     = scales::squish, 
      name    = ""
    )
    
    # 5. Dynamically style each plot theme based on show_axes
    if (show_axes) {
      # 🕵♂️ Debug mode: keep as-is, with a black border and axes to inspect the range
      p <- p + theme_bw() + 
        theme(
          plot.margin = margin(5, 5, 5, 5),
          plot.title = element_text(hjust = 0.5, face = "bold.italic", size = 14, margin = margin(b = 10)),
          legend.position = "right",
          legend.title = element_text(face = "bold", size = 10),
          legend.text = element_text(size = 9),
          panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
          axis.title = element_text(face = "bold", size = 10),
          axis.text = element_text(color = "black", size = 8),
          panel.grid.major = element_line(color = "gray93", size = 0.3),
          panel.grid.minor = element_blank()
        ) +
        labs(x = "Spatial X", y = "Spatial Y")
    } else {
      # 📄 Publication mode: remove the border entirely [modified here]
      p <- p + theme_void() + 
        theme(
          plot.margin = margin(5, 5, 5, 5),
          plot.title = element_text(hjust = 0.5, face = "bold.italic", size = 14, margin = margin(b = 10)),
          legend.position = "right",
          legend.title = element_text(face = "bold", size = 10),
          legend.text = element_text(size = 9),
          # 👈 Key change: replace the original element_rect border with element_blank() to fully remove it
          panel.border = element_blank() 
        )
    }
    
    # 6. Fix the aspect ratio and crop locally
    p <- p + coord_fixed(ratio = 1, xlim = xlim, ylim = ylim)
    
    return(p)
  })
  
  plot_list <- plot_list[!sapply(plot_list, is.null)]
  
  # 7. Combine the plots
  if (is.null(ncol)) ncol <- length(plot_list)
  combined_plot <- wrap_plots(plot_list, ncol = ncol) + 
    plot_layout(guides = "collect") 
  
  return(combined_plot)
}

############################
#######################
#' Bin distance data in a Seurat object into percentiles by group
#' Flexibly adjust percentile distance binning
#' @param seu Seurat object
#' @param dist_col Distance column name
#' @param group_col Group column name
#' @param probs Vector of probability points (e.g. c(0, 0.2, 0.4, 0.6, 0.8, 1))
#' @param new_col_name Name of the newly generated column
add_distance_percentile <- function(seu, 
                                    dist_col = "min_distance", 
                                    group_col = "group", 
                                    probs = seq(0, 1, 0.2),
                                    new_col_name = "distance_percentile") {
  
  meta_data <- seu@meta.data
  
  # 1. Compute the quantiles
  # Use group_by to compute group-specific quantiles
  quantiles <- meta_data %>%
    group_by(.data[[group_col]]) %>%
    summarise(across(all_of(dist_col), ~list(quantile(.x, probs = probs, na.rm = TRUE)))) %>%
    tidyr::unnest_wider(all_of(dist_col), names_sep = "_q")
  
  # 2. Generate label names (e.g. "0-20%", "20-40%"...)
  labels <- sapply(1:(length(probs)-1), function(i) {
    paste0(probs[i]*100, "-", probs[i+1]*100, "%")
  })
  
  # 3. Perform the bin classification
  # Logic: for each row, look up the corresponding quantile boundaries for its group
  meta_data <- meta_data %>%
    left_join(quantiles, by = group_col) %>%
    rowwise() %>%
    mutate(
      !!new_col_name := labels[findInterval(.data[[dist_col]], 
                                            unlist(across(starts_with(paste0(dist_col, "_q")))), 
                                            rightmost.closed = TRUE)]
    ) %>%
    ungroup()
  
  # 4. Write back to the Seurat object
  seu[[new_col_name]] <- meta_data[[new_col_name]]
  
  return(seu)
}



#####################################
##Plot a proportion curve for a continuous variable of a given factor
#' Minimal fast smoothed curve plot (with custom color support)
plot_fast_curve <- function(seu, 
                            x_var, 
                            y_var, 
                            group_var = NULL, 
                            colors = NULL, 
                            max_cells = 10000, 
                            use_relative = FALSE) { # New parameter: relative proportion off by default
  
  # 1. Extract and clean the data
  df <- na.omit(FetchData(seu, vars = c(x_var, y_var, group_var)))
  df[[x_var]] <- as.numeric(df[[x_var]])
  df[[y_var]] <- as.numeric(df[[y_var]])
  
  # 2. Key change: handle the relative proportion of the X-axis
  if (use_relative) {
    min_x <- min(df[[x_var]], na.rm = TRUE)
    max_x <- max(df[[x_var]], na.rm = TRUE)
    
    if (max_x != min_x) {
      # Rescale the X-axis data to between 0 and 1
      df[[x_var]] <- (df[[x_var]] - min_x) / (max_x - min_x)
    } else {
      warning("The maximum and minimum of the X-axis variable are equal, cannot compute relative proportion.")
    }
  }
  
  # 3. Subsample the data (after fixing the relative proportion, to avoid losing extremes)
  if(nrow(df) > max_cells) {
    df <- df[sample(nrow(df), max_cells), ]
  }
  
  # 4. Build the ggplot object
  p <- ggplot(df, aes(x = .data[[x_var]], y = .data[[y_var]]))
  
  if (!is.null(group_var)) {
    df[[group_var]] <- as.factor(df[[group_var]])
    p <- p + aes(color = .data[[group_var]], fill = .data[[group_var]])
  }
  
  # 5. Dynamically adjust the X-axis label text
  x_label <- if (use_relative) paste(x_var, "(Relative Scale 0-1)") else x_var
  
  # 6. Plot and style
  p <- p + 
    geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), alpha = 0.2, size = 1) +
    theme_bw() +
    labs(
      title = paste(y_var, "Dynamics along", x_var), 
      x = x_label, 
      y = y_var
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      panel.grid = element_blank()
    )
  
  # 7. Apply custom colors
  if (!is.null(group_var) && !is.null(colors)) {
    p <- p + scale_color_manual(values = colors) + scale_fill_manual(values = colors)
  }
  
  return(p)
}
