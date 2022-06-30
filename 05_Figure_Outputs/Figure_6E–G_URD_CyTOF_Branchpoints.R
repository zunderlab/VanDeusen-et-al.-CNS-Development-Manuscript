#Corey Williams, University of Virginia
#2019
#Make marker expression vs pseudotime plot for each tip
#Warning: Sometimes timepoint heatmap will be shifted 1-2 pixels out of alignment with expression

print("Start URD_Branchpoints.R")

rm(list = ls(all = TRUE))
.libPaths(c( .libPaths(), "~/R/4.1.1"))

library(URD)
library(viridis)
library(ggplot2)
library(reshape2)
library(cowplot)

# Can't combine more than 2 currently!
URD_DATA_FILES <- c("URD_Final.RData")
COMBINE_ORDER <- list(c(1,2))

# For outputting dendrogram cell plots colored by segment and expression markers
# This should be useful for double-checking things, as a sanity check
OUTPUT_SEG_MARKER_STAGE <- FALSE
SEG_MARKER_STAGE_DIR <- "Segment_Marker_Stage_Plots"

# Current choices for comparison method are "-logPval", "1-pVal", "ksd", "ksp",
# "emd", "diffmean", "diffmedian" ----See code below for how these work
SEG_COMPARE_METHOD <- "-logPval"
STAT_CUTOFF <- 100 # Set to 0 for no cutoff

# Used for trimming marker names, to remove the Metal/Mass info
MARKER_REGEX <- "_[A-Za-z]+[0-9]+Di"

# Plotting order is generated manually.  You can use these to change it.
GLOBAL_MARKER_ORDER <- NULL # Leave NULL unless you want to specify
MARKER_ORDER <- NULL # Leave NULL unless you want to specify
PARENT_ORDER <- NULL # Leave NULL unless you want to specify
CHILD_ORDER <- NULL # Leave NULL unless you want to specify

# For the rectangles on top of the plot that indicate segments
PARENT_RECT_HEIGHT <- 1
PARENT_RECT_TEXT_SIZE <- 4

# Control dot size
DOT_SIZE_MIN = 0
DOT_SIZE_MAX = 7
DOT_SIZE_INF = 1.5 # How big to make the Inf stat dots
DOT_SCALE_MIN = 0 # percentage out of (Max-Min)
DOT_SCALE_MAX = 100 # percentage out of (Max-Min)
DOT_SCALE = 'size' # DOT_SCALE must be either 'size' or 'radius'

# Use this to "ramp up" the color, (make things appear more yellow overall on the viridis scale).
COLOR_SCALE_MAX = 50 # percentage of max expression to set as max color

# Output PNG size
PNG_WIDTH = 6
PNG_HEIGHT= 10

# Set up segment compare function, depending on user choice
if (SEG_COMPARE_METHOD=="-logPval") {
  seg_compare <- function(s1, s2) {
    # Add weight column because t.test() needs this
    s1mat <- cbind(weight=1, location=s1)
    s2mat <- cbind(weight=1, location=s2)
    comparison <- -log(t.test(s1mat,s2mat)$p.value)
    return(comparison)
  }
}

if (SEG_COMPARE_METHOD=="1-pVal") {
  seg_compare <- function(s1, s2) {
    # Add weight column because t.test() needs this
    s1mat <- cbind(weight=1, location=s1)
    s2mat <- cbind(weight=1, location=s2)
    comparison <- 1-t.test(s1mat,s2mat)$p.value
    return(comparison)
  }
}

if (SEG_COMPARE_METHOD=="ksd") {
  seg_compare <- function(s1, s2) {
    # Add weight column because t.test() needs this
    comparison <- unname(ks.test(s1,s2)$statistic)
    return(comparison)
  }
}

if (SEG_COMPARE_METHOD=="ksp") {
  seg_compare <- function(s1, s2) {
    # Add weight column because t.test() needs this
    comparison <- -log(unname(ks.test(s1,s2)$p.value))
    return(comparison)
  }
}

if (SEG_COMPARE_METHOD=="emd") {
  library(emdist)
  seg_compare <- function(s1, s2) {
    # convert to matrices because emd() needs this
    s1mat <- cbind(weight=1, location=s1)
    s2mat <- cbind(weight=1, location=s2)
    comparison <- emd(s1mat,s2mat)
    return(comparison)
  }
}

if (SEG_COMPARE_METHOD=="diffmean") {
  seg_compare <- function(s1, s2) {
    # convert to matrices because emd() needs this
    comparison <- abs(mean(s1) - mean(s2))
    return(comparison)
  }
}

if (SEG_COMPARE_METHOD=="diffmedian") {
  seg_compare <- function(s1, s2) {
    # convert to matrices because emd() needs this
    comparison <- abs(median(s1) - median(s2))
    return(comparison)
  }
}

# Initialize variables that will store info needed for direct comparison between the two URD dendrograms
seg_counter <- 0
g_plot_df <- list()
g_rect_df <- list()
g_re_order_rect <- list()
g_fill <- list()
g_label <- list()
g_div_lines <- list()
g_roots <- list()

for (u in URD_DATA_FILES) {
  
  cat("Processing ", u, "\n")
  
  load(u)
  markers <- URD_Object.tree@var.genes
  
  # Output figures to confirm segment connectivity (and marker expression while we're at it).
  if (OUTPUT_SEG_MARKER_STAGE) {
    dir.create(SEG_MARKER_STAGE_DIR)
    
    png(paste(SEG_MARKER_STAGE_DIR, "URD_Tree_by_stage.png", sep="/"))
    plotTree(URD_Object.tree, "stage", title="Developmental Stage")
    dev.off()
    png(paste(SEG_MARKER_STAGE_DIR, "URD_Tree_by_segment.png", sep="/"))
    plotTree(URD_Object.tree, "segment", title="URD tree segment")
    dev.off()
    for (m in markers) {
      png(paste0(SEG_MARKER_STAGE_DIR, "/", "URD_Tree_by_",m,".png"))
      p <- plotTree(URD_Object.tree, m, title=m)
      print(p)
      dev.off()
    }
  }
  
  URD_expr_data <- URD_Object.tree@count.data
  
  seg_expr_means <- c() # initialize
  bifurcs <- list() # initialize
  for (s in URD_Object.tree@tree$segments) {
    
    # Get expression means for each segment
    seg_expr <- URD_expr_data[,URD_Object.tree@tree$cells.in.segment[s][[1]]]
    seg_expr_means <- cbind(seg_expr_means, apply(seg_expr, 1, mean))
    
    # Get list of all  bifurcations, storing parent and children for each
    c <- segChildren(URD_Object.tree, s)
    if (length(c) > 1) {
      bifurcs[[s]] <- c
    }
  }
  colnames(seg_expr_means) <- URD_Object.tree@tree$segments
  
  # Calculate t-test p-values between child segments for each bifurcation
  # initialize dataframe to store p-values
  df_colnames <- c("parent", "segment1", "segment2", markers)
  diff_markers <- data.frame(matrix(ncol=length(df_colnames), nrow = 0))
  colnames(diff_markers) <- df_colnames
  for (s in names(bifurcs)) { # for each bifurcation
    for (i in 1:(length(bifurcs[[s]])-1)) { # doing this complicated i-j setup for when there are >2 segments at a bifurcation point (trifurcation, etc.)  This will do all pairwise comparisons
      for (j in (i+1):length(bifurcs[[s]])) { 
        marker_scores <- list() # initialize
        for (m in markers) {
          # Get all data from these segments/markers
          s1 <- URD_expr_data[m,URD_Object.tree@tree$cells.in.segment[bifurcs[[s]][i]][[1]]]
          s2 <- URD_expr_data[m,URD_Object.tree@tree$cells.in.segment[bifurcs[[s]][j]][[1]]]
          
          # Do the comparison, based on method choice aboce
          marker_scores[[m]] <- seg_compare(s1,s2)
        }
        # Build parent/child dataframe
        segnames_df <- data.frame(parent=s,
                                  child1=bifurcs[[s]][i],
                                  child2=bifurcs[[s]][j])
        # Build stat compare dataframe for each marker
        stat_compare_df <- as.data.frame(t(unlist(marker_scores)))
        # Combind parent/child and stat compare dataframes, and store
        diff_markers <- rbind(diff_markers, cbind(segnames_df, stat_compare_df))
      }
    }
  }
  
  # Take care of Inf values (from super low p-Values)
  diffs_vector <- unlist(diff_markers[,markers])
  diff_max <- max(diffs_vector[which(is.finite(diffs_vector))])
  
  diff_markers <- do.call(data.frame, lapply(diff_markers, function(x) replace(x, is.infinite(x), diff_max*DOT_SIZE_INF)))
  
  # This is used to make comparisons always go from left to right
  seg_xs <- URD_Object.tree@tree$segment.layout$x
  names(seg_xs) <- URD_Object.tree@tree$segment.layout[,"segment"]
  
  if (is.null(PARENT_ORDER)) {
    # Initialize object to store segment order
    parent_order = c()
    # Define recursive function to walk down the dendrogram and order segments
    order_segments <- function(mat, seg) {
      parent_order <<- c(parent_order, seg) # special equal sign to make this global instead of local
      rec_mat <- mat[which(!as.numeric(mat[,1])==seg),] # cut out the current segment row
      # Reorder so it always goes left to right
      child_segs <- mat[which(as.numeric(mat[,1])==seg),2:3]
      child_segs <- as.character(unique(unlist(child_segs))) # need this for trifurcations
      child_segs <- child_segs[order(seg_xs[as.character(child_segs)])]
      
      # Iterate on all child segments
      for (s in child_segs) {
        # Only recurse if the child segment has its own children
        if (s %in% mat[,1]) {
          order_segments(rec_mat,s)
        }
      }
    }
    
    # Calculate the root index
    root_index <- which(!(diff_markers[,"parent"] %in%
                            unlist(diff_markers[,c("child1", "child2")])))
    # Need to use unique() here because of trifurcation branchpoints
    root_segment <- as.numeric(unique(diff_markers[root_index,"parent"]))
    # Calculate parent order recursively
    order_segments(diff_markers, root_segment)
  } else {
    parent_order <- PARENT_ORDER
  }
  
  if (is.null(MARKER_ORDER)) {
    # Calculate marker order for plotting
    # First get parent order by index in diff_markers
    parent_idx_order <- unlist(lapply(parent_order, function(x) which(diff_markers[,"parent"]==x)))
    
    marker_order <- c() # initialize
    marker_top_ranked <- apply(diff_markers[,markers], 2, which.max) # get max idxs
    for (i in rev(parent_idx_order)) {
      seg_markers <- markers[which(marker_top_ranked==i)]
      marker_order <- c(marker_order,
                        names(sort(diff_markers[i,seg_markers, drop=FALSE])))
    }
  } else {
    marker_order <- MARKER_ORDER
  }
  
  if (is.null(CHILD_ORDER)) {
    stat_compare_table <- c() # initialize
    child_order <- c() # initialize
    for (p in parent_order) {
      # First re-order, to make sure it goes from left to right
      segs <- segChildren(URD_Object.tree, p)
      segs <- segs[order(seg_xs[segs])]
      for (s in segs) {
        # Get index/indices of this parent/child combination (need >1 for trifurcation!)
        p_c_idx <- which(diff_markers[,"parent"]==p &
                           (diff_markers[,"child1"]==s | diff_markers[,"child2"]==s))
        # Get the comparison stats.  Will be one column for bifurcation, 2 columns for trifurcation
        compare_stats <- t(diff_markers[p_c_idx,marker_order])
        # Take max from each row - need this for trifurcations
        compare_stats <- apply(compare_stats, 1, max)
        # Store in table, for dot size
        stat_compare_table <- cbind(stat_compare_table, compare_stats)
        colnames(stat_compare_table)[ncol(stat_compare_table)] <- s
        # Store child order too, for plotting order
        child_order <- c(child_order, s)
      }
    }
  } else {
    child_order <- CHILD_ORDER
  }
  
  # Assemble the dataframe needed for ggplot
  mean_table <- seg_expr_means[marker_order, child_order]
  
  stat_compare_df <- melt(stat_compare_table)
  colnames(stat_compare_df) <- c("Marker", "Segment", "Stat")
  rownames(stat_compare_df) <- paste(stat_compare_df[,"Segment"], stat_compare_df[,"Marker"], sep="-") # Sanity check for merging.  Might be able to catch mismatched here?
  
  mean_df <- melt(mean_table)
  colnames(mean_df) <- c("Marker", "Segment", "Mean")
  rownames(mean_df) <- paste(mean_df[,"Segment"], mean_df[,"Marker"], sep="-")
  
  plot_df <- cbind(stat_compare_df, Mean=mean_df[,"Mean"])
  plot_df$Segment <- factor(plot_df$Segment, levels=child_order)
  
  scale.func <- switch(
    EXPR = DOT_SCALE,
    'size' = scale_size,
    'radius' = scale_radius,
    stop("'DOT_SCALE' must be either 'size' or 'radius'")
  )
  
  # Adjust p-Val range to use for scaling dots 
  dot_scale_min <- min(diff_markers[,markers]) + (max(diff_markers[,markers])-min(diff_markers[,markers]))*DOT_SCALE_MIN/100
  dot_scale_max <- min(diff_markers[,markers]) + max(diff_markers[,markers])*DOT_SCALE_MAX/100
  
  # Set up borders for dividing lines and colored parent rectangles
  dividing_lines <- c()
  parent_label_locs <- c()
  counter = 0.5
  for (b in parent_order) {
    new_loc <- counter + length(bifurcs[[which(names(bifurcs)==b)]])/2
    parent_label_locs <- c(parent_label_locs, new_loc)
    
    counter <- counter + length(bifurcs[[b]])
    dividing_lines <- c(dividing_lines, counter)
  }
  
  # Set up rectangles that will be used to label parent segments
  rect_df <- data.frame(
    x1=c(0.5,dividing_lines[1:(length(dividing_lines)-1)]),
    x2=dividing_lines,
    y1=rep(length(markers)+0.5, length(parent_order)),
    y2=rep(length(markers)+0.5+PARENT_RECT_HEIGHT, length(parent_order)),
    Parent=parent_order
  )
  
  # Set up colors for rectangles and dendrogram labels
  seg_cols <- scales::hue_pal()(length(URD_Object.tree@tree$segments))
  names(seg_cols) <- (1:length(URD_Object.tree@tree$segments))+seg_counter
  
  # Get new dendrogram segment order (based on going down the tree) for segment re-naming
  root_segment <- as.character(root_segment) # Get root segment
  segment_order <- c(root_segment) # Initialize object to store segment order
  # Define recursive function to get the dendrogram ordering we want
  order_segments <- function(s) {
    c <- segChildren(URD_Object.tree, s)
    c <- c[order(seg_xs[c])]
    if (length(c) > 1)
      for (i in c) {
        segment_order <<- c(segment_order, i)
        order_segments(i)
      }
  }
  # Run the recursive function to get the dendrogram ordering we want
  order_segments(root_segment)
  
  # More re-ordering to convert from original to new segment numbering
  starting_seg_order <- URD_Object.tree@tree$segment.layout[,"segment"]
  re_order_tree <- c()
  for (s in starting_seg_order) {
    re_order_tree <- c(re_order_tree, which(segment_order==s))
  }
  re_order_tree <- re_order_tree + seg_counter
  re_order_tree <- as.character(re_order_tree)
  names(re_order_tree) <- (1:length(re_order_tree)) + seg_counter
  
  re_order_rect <- (1:length(re_order_tree)) + seg_counter
  names(re_order_rect) <- segment_order
  re_order_rect[as.character(unique(plot_df$Segment))]
  
  # ggplot output of dot plot
  gg_dot_plot <- ggplot(plot_df, mapping=aes(x=Segment, y=Marker)) +
    geom_point(mapping = aes_string(size='Stat', color ='Mean')) +
    scale_x_discrete(expand = c(0, 0), labels=re_order_rect[as.character(unique(plot_df$Segment))]) +
    scale_color_viridis(limits=c(0, max(plot_df$Mean)*COLOR_SCALE_MAX/100),
                        oob=scales::squish) +
    scale.func(range = c(DOT_SIZE_MIN, DOT_SIZE_MAX),
               limits = c(dot_scale_min, dot_scale_max)) + #optional scaling method
    theme_cowplot() +
    theme(axis.line=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          panel.border=element_rect(colour="black", fill=NA, size=1)) +
    geom_rect(rect_df, inherit.aes = FALSE, mapping=aes(xmin=x1, xmax=x2,
                                                        ymin=y1, ymax=y2),
              fill=seg_cols[as.character(re_order_rect[rect_df[,'Parent']])],
              show.legend=FALSE) +
    geom_text(rect_df, inherit.aes = FALSE, mapping=aes(x=x1+(x2-x1)/2,
                                                        y=y1+(y2-y1)/2,
                                                        label=re_order_rect[as.character(rect_df[,'Parent'])]),
              size=PARENT_RECT_TEXT_SIZE) +
    guides(size = guide_legend(title = SEG_COMPARE_METHOD),
           shape = "none") +  #sets the legend titles, but setting color title changes colorbar to dots. . .
    geom_vline(xintercept = dividing_lines[1:(length(dividing_lines)-1)]) +
    geom_hline(yintercept = length(markers)+0.5)
  dot_png_name <- paste0("Branchpoint_Dots_", SEG_COMPARE_METHOD, "_", u, ".png")
  ggsave(dot_png_name, gg_dot_plot, width=1.7*PNG_WIDTH, height=1.1*PNG_HEIGHT, units="in")
  
  # Store info for global comparison
  g_plot_df[[length(g_plot_df)+1]] <- plot_df
  g_rect_df[[length(g_rect_df)+1]] <- rect_df
  g_re_order_rect[[length(g_re_order_rect)+1]] <- re_order_rect
  g_fill[[length(g_fill)+1]] <- seg_cols[as.character(re_order_rect[rect_df[,'Parent']])]
  g_label[[length(g_label)+1]] <- re_order_rect[as.character(rect_df[,'Parent'])]
  g_div_lines[[length(g_div_lines)+1]] <- dividing_lines
  g_roots[[length(g_roots)+1]] <- URD_expr_data[,URD_Object.tree@tree$cells.in.segment[root_segment][[1]]]
  
  # ggplot output of labeled/colored dendrogram
  # First grab these tree and segment (label) layouts from the object
  segment.layout <- URD_Object.tree@tree$segment.layout
  tree.layout <- URD_Object.tree@tree$tree.layout
  # Now set up segment label coordinates
  segment.labels <- as.data.frame(segment.layout[,c("segment","x")])
  segment.labels$y <- apply(URD_Object.tree@tree$segment.pseudotime.limits, 1, num.mean)[segment.labels$segment]
  # Now plot
  gg_seg_plot <- ggplot() +
    scale_y_reverse(c(1,0), name="Pseudotime", breaks=seq(0, 1, 0.1)) +
    theme_bw() +
    theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), legend.position="none") +
    geom_segment(data=tree.layout, aes(x=x1, y=y1, xend=x2, yend=y2), color='black', alpha=1, size=1, lineend="square") +
    geom_label(data=segment.labels, aes(x=x, y=y, label=re_order_tree, fill=re_order_tree), alpha=1) +
    scale_fill_manual(values=seg_cols)
  tree_png_name <- paste0("Branchpoint_Tree_", SEG_COMPARE_METHOD, "_", u, ".png")
  ggsave(tree_png_name, gg_seg_plot, width=PNG_WIDTH, height=PNG_HEIGHT, units="in")
  
  seg_counter <- seg_counter + length(segment_order)
}

# Set up combined dataframes for dots (mean/p-Val) and rectangles for segment ID
c_plot_df <- c()
c_rect_df <- c()
for (i in 1:length(g_plot_df)) {
  l_plot_df <- g_plot_df[[i]]
  l_rect_df <- g_rect_df[[i]]
  l_re_order <- g_re_order_rect[[i]]
  
  rename_plot_df <- list()
  rename_rect_df <- list()
  for (n in names(l_re_order)) {
    rename_plot_df[[n]] <- which(l_plot_df$Segment==n)
    rename_rect_df[[n]] <- which(l_rect_df$Parent==n)
  }
  
  l_plot_df_segment <- as.character(l_plot_df$Segment)
  for (n in names(rename_plot_df)) {
    l_plot_df_segment[rename_plot_df[[n]]] <- l_re_order[n]
  }
  l_plot_df$Segment <- factor(x=l_plot_df_segment, levels=unique(l_plot_df_segment))
  
  for (n in names(rename_rect_df)) {
    l_rect_df$Parent[rename_rect_df[[n]]] <- l_re_order[n]
  }
  
  c_plot_df <- rbind(c_plot_df, l_plot_df)
  c_rect_df <- rbind(c_rect_df, l_rect_df)
}

# Get matrices back from plotting dataframe
c_stat_mat <- acast(c_plot_df, c_plot_df[,"Marker"]~c_plot_df[,"Segment"], value.var="Stat")
c_mean_mat <- acast(c_plot_df, c_plot_df[,"Marker"]~c_plot_df[,"Segment"], value.var="Mean")

# Get max values for cutoff
max_stats <- apply(c_stat_mat, 1, max)
below_cutoff <- names(max_stats[which(max_stats < STAT_CUTOFF)])

# Add in comparison of root segments! Can't use this with current dataset
# because there are no cells in TrkB/TrkC root segment!
# ------------------This section is not used in plotting for now
r_marker_scores <- list()
for (m in markers) {
  # Get all data from these segments/markers
  s1 <- g_roots[[1]][m,]
  s2 <- g_roots[[2]][m,]
  
  # Do the comparison, based on method choice aboce
  r_marker_scores[[m]] <- seg_compare(s1,s2)
}

# Get global marker order - from the two URD dendrograms combined
if (is.null(GLOBAL_MARKER_ORDER)) {
  c_marker_order <- c() # initialize
  c_marker_top_ranked <- apply(c_stat_mat, 1, which.max) # get max idxs
  for (i in seq(from=2, to=ncol(c_stat_mat), by=2)) {
    tops <- which(c_marker_top_ranked==i | c_marker_top_ranked==i-1)
    c_seg_markers <- names(tops)
    c_marker_order <- c(c_marker_order, rev(names(tops[order(c_stat_mat[tops,i, drop=FALSE])])))
  }
  c_marker_order <- rev(c_marker_order)
} else {
  c_marker_order <- MARKER_ORDER
}

# Remove low p-Value markers below the cutoff STAT_CUTOFF
c_marker_order <- c_marker_order[which(!(c_marker_order %in% below_cutoff))]
c_stat_ordered <- c_stat_mat[c_marker_order,]

# Clean up parameter names
rownames(c_stat_ordered) <- sub(MARKER_REGEX, "", rownames(c_stat_ordered)) 

# Prepare combined dataframe for ggplot
c_stat_df <- melt(c_stat_ordered)
colnames(c_stat_df) <- c("Marker", "Segment", "Stat")
rownames(c_stat_df) <- paste(c_stat_df[,"Segment"], c_stat_df[,"Marker"],
                             sep="-") # Sanity check for merging.  Might be able to catch mismatched here?

c_mean_df <- melt(c_mean_mat[c_marker_order,])
colnames(c_mean_df) <- c("Marker", "Segment", "Mean")
rownames(c_mean_df) <- paste(c_mean_df[,"Segment"], c_mean_df[,"Marker"], sep="-")

c_plot_df <- cbind(c_stat_df, Mean=c_mean_df[,"Mean"])
c_plot_df$Segment <- factor(c_plot_df$Segment, levels=unique(c_plot_df$Segment))

# Set up rectangles and dividing lines for plot
rect_counter <- 0
div_counter <- 0
c_rect_df <- c()
c_div_lines <- c()
for (i in 1:length(g_rect_df)) {
  l_rect_df <- g_rect_df[[i]]
  l_rect_df[,c("x1", "x2")] <- l_rect_df[,c("x1", "x2")] - 0.5 + rect_counter
  c_rect_df <- rbind(c_rect_df, l_rect_df)
  
  l_div_lines <- g_div_lines[[i]] - 0.5 + div_counter
  c_div_lines <- c(c_div_lines, l_div_lines)
  
  rect_counter <- rect_counter + l_rect_df[nrow(l_rect_df),"x2"]
  div_counter <- div_counter + l_div_lines[length(l_div_lines)]
}

c_rect_df[,c("x1", "x2")] <- c_rect_df[,c("x1", "x2")] + 0.5
c_div_lines <- c_div_lines + 0.5
c_rect_df[,c("y1", "y2")] <- c_rect_df[,c("y1", "y2")] - length(below_cutoff)

# Recover color and segment names to use in rectangles
c_fill <- c()
c_label <- c()
for (i in 1:length(g_fill)) {
  c_fill <- c(c_fill, g_fill[[i]])
  c_label <- c(c_label, g_label[[i]])
}

# Adjust p-Val range to use for scaling dots 
dot_scale_min <- min(c_stat_mat) + (max(c_stat_mat)-min(c_stat_mat))*DOT_SCALE_MIN/100
dot_scale_max <- min(c_stat_mat) + max(c_stat_mat)*DOT_SCALE_MAX/100

# ggplot output
p <- ggplot(c_plot_df, mapping=aes(x=Segment, y=Marker)) +
  geom_point(mapping = aes_string(size='Stat', color ='Mean')) +
  scale_x_discrete(expand = c(0, 0), labels=unique(c_plot_df$Segment)) +
  scale_color_viridis(limits=c(0, max(c_plot_df$Mean)*COLOR_SCALE_MAX/100),
                      oob=scales::squish) +
  scale.func(range = c(DOT_SIZE_MIN, DOT_SIZE_MAX),
             limits = c(dot_scale_min, dot_scale_max)) + #optional scaling method
  theme_cowplot() +
  theme(axis.line=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        panel.border=element_rect(colour="black", fill=NA, size=1)) +
  geom_rect(c_rect_df, inherit.aes = FALSE,
            mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2),
            fill=c_fill,
            show.legend=FALSE) +
  geom_text(c_rect_df, inherit.aes = FALSE,
            mapping=aes(x=x1+(x2-x1)/2,
                        y=y1+(y2-y1)/2,
                        label=c_label),
            size=PARENT_RECT_TEXT_SIZE) +
  guides(size = guide_legend(title = SEG_COMPARE_METHOD),
         shape = "none") +  #sets the legend titles, but setting color title changes colorbar to dots. . .
  geom_vline(xintercept = c_div_lines[1:(length(c_div_lines)-1)]) +
  geom_hline(yintercept = length(c_marker_order)+0.5)
dot_png_name <- paste0("Branchpoint_Dots_Combined_", SEG_COMPARE_METHOD, "_", u, ".png")
ggsave(dot_png_name, p, width=1.6*PNG_WIDTH, height=PNG_HEIGHT, units="in")

print("Finish URD_Branchpoints.R")