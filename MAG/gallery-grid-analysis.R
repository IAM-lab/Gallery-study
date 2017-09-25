##################################################################################
# NAME:         gallery-grid-analysis.R
# AUTHOUR:      Alan Davies
# DATE:         15/09/2017
# INSTITUTION:  Interaction Analysis and Modelling Lab (IAM), University of Manchester
# DESCRIPTION:  Version of Markov chain transition analysis using a grid system
#               based on DBSCAN cluster size using Bhattacharyya distance.
##################################################################################

#---------------------------------------------------------------------------------
# FUNCTION:     loadDataFiles()
# INPUT:        void
# OUTPUT:       list
# DESCRIPTION:  Returns a with data frames for the 2 test conditions filtered
#               for unique fixation index events
#---------------------------------------------------------------------------------
loadDataFiles <- function()
{
    data_files <- list()
    cond1_data_path <- paste0(getwd(), "/gallery-study/data/condition1.csv")
    cond2_data_path <- paste0(getwd(), "/gallery-study/data/condition2.csv")
    
    # conditions
    data_file <- read.csv(cond1_data_path, header = TRUE, na.strings = c(" ", "NA", "-"))
    data_files[["condition1"]] <- data_file
    data_file <- read.csv(cond2_data_path, header = TRUE, na.strings = c(" ", "NA", "-"))
    data_files[["condition2"]] <- data_file
    
    return(data_files)
}

#---------------------------------------------------------------------------------
# FUNCTION:     convertToMatrix()
# INPUT:        data.frame
# OUTPUT:       matrix
# DESCRIPTION:  Returns matrix representing probabilty of AOI transitions
#               
#---------------------------------------------------------------------------------
convertToMatrix <- function(data)
{
    tmp0 <- data
    tmp <- tmp0 %>% group_by(participant) %>% mutate(to = lead(AOI))
    tmp2 <- tmp[complete.cases(tmp), ]
    with(tmp2, table(AOI, to))
    out_mat <- as.matrix(with(tmp2, table(AOI, to)))
    diag(out_mat) <- 0
    return(out_mat)
}

#---------------------------------------------------------------------------------
# FUNCTION:     convertToProbMatrix()
# INPUT:        matrix
# OUTPUT:       matrix
# DESCRIPTION:  Converts to probability matrix by dividing by rows sums
#               
#---------------------------------------------------------------------------------
convertToProbMatrix <- function(data)
{
    prob_mat <- data / rowSums(data)
    prob_mat[is.na(prob_mat)] <- 0
    return(prob_mat)
}

#---------------------------------------------------------------------------------
# FUNCTION:     trimMatrixZeros()
# INPUT:        matrix, matrix
# OUTPUT:       list
# DESCRIPTION:  Returns sub matrix removing rows/cols where rows and columns 
#               in both matrices  are all zeros
#---------------------------------------------------------------------------------
trimMatrixZeros <- function(m1, m2)
{
    matrix_list <- list()
    row_col_to_remove <- NULL
    
    for(x in 1:nrow(m1))
    {
        if((rowSums(m1)[x] == 0) && (rowSums(m2)[x] == 0))
        {  
            if((colSums(m1)[x] == 0) && (colSums(m2)[x] == 0))
            {
                row_col_to_remove <- c(row_col_to_remove, x)
            }
        }
    } 
    
    if(length(row_col_to_remove) > 0)
    {
        # removed the selected rows and columns
        m1 <- m1[-row_col_to_remove, ]; m1 <- m1[, -row_col_to_remove] 
        m2 <- m2[-row_col_to_remove, ]; m2 <- m2[, -row_col_to_remove] 
        
        # store in list entries
        matrix_list[[1]] <- m1
        matrix_list[[2]] <- m2
        
        return(matrix_list)
    }
    return(FALSE)
}

#---------------------------------------------------------------------------------
# FUNCTION:     loadPackages(package.args)
# INPUT:        vector
# OUTPUT:       void
# DESCRIPTION:  Loads required packages.
#                
#---------------------------------------------------------------------------------
loadPackages <- function(package.args)
{ 
    for(i in package.args)
    {
        if(!is.element(i, .packages(all.available = TRUE)))
        {
            cat("\nPackage <", i, "> not found, attempting to add it...")
            install.packages(i)
        }
        library(i, character.only = TRUE)
    }
}

#---------------------------------------------------------------------------------
# FUNCTION:     initialize()
# INPUT:        void
# OUTPUT:       void
# DESCRIPTION:  Set up function for adding packages and other source data
#               
#---------------------------------------------------------------------------------
initialize <- function()
{
    # load packages
    package.args <- c("msm", "png", "gplots", "dplyr", "plyr", "matrixStats", "topicmodels", "lattice")
    loadPackages(package.args)
}

#---------------------------------------------------------------------------------
# FUNCTION:     determineImageOffset()
# INPUT:        list, list
# OUTPUT:       vector (horizontal, vertical offsets)
# DESCRIPTION:  Determine the image offsets taking account of the large black borders
#               Tobii applies to the images. Start centre axis and work in until pixel
#               RGB values are not black.
#---------------------------------------------------------------------------------
determineImageOffset <- function(data, stimuli)
{
    # load and plot the stimulus image
    stimuli_image_str <- gsub(".jpg", "", stimuli$file_name)
    stimuli_image_str <- gsub(".png", "", stimuli_image_str)
    bkg_img <- readPNG(paste0(getwd(), "/gallery-study/stimuli/", stimuli_image_str, ".png"), info = TRUE)
    plot.new()
    limits <- par()
    
    # get the dimensions
    dims <- attr(bkg_img, "dim")
    
    # get vertical offset
    for(i in 1:dims[1])
    {
        # get the center point of image and pixel value
        row_offset <- floor(dims[2] / 2)
        px <- bkg_img[i, row_offset, 1:3]
        
        # get offset if pixel colour is not black
        if((px[1] > 0) || (px[2] > 0) || (px[3] > 0))
        {
            vertical_offset <- (i - 1)
            break
        }
    }
    
    # get horizontal offset
    for(i in 1:dims[2])
    {
        # get the center point of image and pixel value
        col_offset <- floor(dims[1] / 2)
        px <- bkg_img[col_offset, i, 1:3]
        
        # get offset if pixel colour is not black
        if((px[1] > 0) || (px[2] > 0) || (px[3] > 0))
        {
            horizontal_offset <- (i - 1)
            break
        }
    }
    return(c(horizontal_offset, vertical_offset))
}

#---------------------------------------------------------------------------------
# FUNCTION:     drawGrid()
# INPUT:        list, list
# OUTPUT:       
# DESCRIPTION:  
#               
#---------------------------------------------------------------------------------
drawGrid <- function(grid_args, stimuli)
{
    # load and plot the stimulus image
    stimuli_image_str <- gsub(".jpg", "", stimuli$file_name)
    stimuli_image_str <- gsub(".png", "", stimuli_image_str)
    bkg_img <- readPNG(paste0(getwd(), "/gallery-study/stimuli/", stimuli_image_str, ".png"), info = TRUE)
    plot.new()
    limits <- par()
    rasterImage(bkg_img, limits$usr[1], limits$usr[3], limits$usr[2], limits$usr[4])
    
    dims <- attr(bkg_img, "dim")
}

#---------------------------------------------------------------------------------
# FUNCTION:     determineGridOffset()
# INPUT:        int, int, int
# OUTPUT:       vector (x, y offset)
# DESCRIPTION:  Determine the x and y offsets to center the grid in the 
#               middle of the image
#---------------------------------------------------------------------------------
determineGridOffset <- function(width, height, grid_len)
{
    x_offset <- (width / grid_len) - floor(width / grid_len)
    y_offset <- (height / grid_len) - floor(height / grid_len)
    return(c(x_offset, y_offset))
}

#---------------------------------------------------------------------------------
# FUNCTION:     removeRepeatFixation()
# INPUT:        data.frame
# OUTPUT:       data.frame
# DESCRIPTION:  Remove any duplicate fixations 
#               
#---------------------------------------------------------------------------------
removeRepeatFixation <- function(data)
{
    filtered_data <- data[FALSE, ]
    data <- data[data$GazeEventType == "Fixation", ]
    participants <- as.vector(data[!duplicated(data$ParticipantName), "ParticipantName"])
    
    for(i in 1:length(participants))
    {
        participant_data <- data[data$ParticipantName == participants[i], ]
        participant_data <- participant_data[!duplicated(participant_data$FixationIndex), ]
        filtered_data <- rbind(filtered_data, participant_data)
    }
    return(filtered_data)
}

#---------------------------------------------------------------------------------
# FUNCTION:     detectCellHits()
# INPUT:        int, int, list
# OUTPUT:       int
# DESCRIPTION:  Detect the number of hits for the current stimulus in the
#               designated grid cell area
#---------------------------------------------------------------------------------
detectCellHits <- function(cell_x, cell_y, cell_len, data)
{
    # extract data
    data <- data[1, c("FixationPointX", "FixationPointY")]
    fixation_x <- data[1, 1]
    fixation_y <- data[1, 2]
    
    # return no hit for any missing data
    if(is.na(fixation_x) || is.na(fixation_y))
    {
        return(FALSE)
    }
    
    # check present inside cell 
    if((fixation_x > cell_x && fixation_x < (cell_x + cell_len)) && 
       (fixation_y > cell_y && fixation_y < (cell_y + cell_len)))
    {
        return(TRUE)
    }
    return(FALSE)
}

#---------------------------------------------------------------------------------
# FUNCTION:     addAOIlabels()
# INPUT:        int
# OUTPUT:       void
# DESCRIPTION:  Lables from alphabet i.e. A-Z then AA etc. 
#               
#---------------------------------------------------------------------------------
addAOIlabels <- function(letter_len) 
{
    a <- rep(LETTERS, length.out = letter_len)
    grp <- cumsum(a == "A")
    vapply(seq_along(a), function(x) paste(rep(a[x], grp[x]), collapse = ""), character(1L))
}

#---------------------------------------------------------------------------------
# FUNCTION:     generateTransitionData()
# INPUT:        data.frame, int, list
# OUTPUT:       data.frame
# DESCRIPTION:  Builds the sequence of AOIs visited by the participants
#               
#---------------------------------------------------------------------------------
generateTransitionData <- function(data, condition_num, grid_args, perm_trans = FALSE)
{
    data_frames <- list()
    mat <- matrix(addAOIlabels(grid_args$num_cells), nrow = grid_args$num_rows, ncol = grid_args$num_cols, byrow = TRUE)
    
    participant_vector <- NULL
    AOI_vector <- NULL
    
    # get specified condition
    if(!perm_trans){
        condition <- data[data$Condition == condition_num, ]
    } else {
        condition <- data
    }
    
    # get the participants for the selected condition
    participants <- as.vector(condition[!duplicated(condition$ParticipantName), "ParticipantName"])
    
    for(i in 1:length(participants))
    {
        # get participant
        participant_data <- data[data$ParticipantName == participants[i], ]
        
        for(j in 1:nrow(participant_data))
        {
            for(m in grid_args$num_rows:1)
            {
                for(n in 1:grid_args$num_cols)
                {
                    # calculate relative position of cell
                    cell_x <- grid_args$images_offsets[1] + grid_args$grid_offsets[1] + (grid_args$dims * (n - 1))
                    cell_y <- grid_args$images_offsets[2] + grid_args$grid_offsets[2] + (grid_args$dims * (m - 1))
                    
                    # detect hit inside cell
                    hit <- detectCellHits(cell_x, cell_y, grid_args$dims, participant_data[j, ])
                    if(hit)
                    {
                        AOI_vector <- c(AOI_vector, mat[m, n])
                    }
                }
            }
        }
        AOIs <- addAOIlabels(grid_args$num_cells)
        participant_vector <- rep(participants[i], length(AOI_vector))
        df <- data.frame(participant = participant_vector, AOI = AOI_vector)
        df$AOI <- factor(df$AOI, levels = AOIs)
        
        # add data frame to list
        data_frames[[i]] <- df
        participant_data <- NULL
    }
    hit_data <- ldply(data_frames, data.frame)
    
    return(hit_data)    
}

#---------------------------------------------------------------------------------
# FUNCTION:     randomPermutations()
# INPUT:        data.frame, int, int, int
# OUTPUT:       vector
# DESCRIPTION:  Generate 2 random sub groups the same size as the original groups
#               and compare distances for the selected amout of permutations.
#---------------------------------------------------------------------------------
randomPermutations <- function(data, group1_size, grid_args, perms = 10)
{
    progress <- 0
    distances <- numeric(perms)
    
    # extract participants
    participants <- as.data.frame(data[!duplicated(data$ParticipantName), "ParticipantName"])
    colnames(participants) <- "participant"
    
    # create progress bar
    progress_bar <- winProgressBar(title = "Computing Distances", min = 0, max = perms, width = 300)
    getWinProgressBar(progress_bar)
    
    for(i in 1:perms)
    {
        # get a random subset of participants and store in group 1 whatever is left put in second group
        group1_participants <- participants[sample(unique(nrow(participants)), group1_size), ]
        group2_participants <- as.vector(participants[!(participants$participant %in% group1_participants), ])
        
        # extract groups
        group1 <- data[which(data$ParticipantName %in% group1_participants), ]
        group2 <- data[which(data$ParticipantName %in% group2_participants), ]
        
        # generate transition data
        group1_trans <- generateTransitionData(group1, NA, grid_args, perm_trans = TRUE)
        group2_trans <- generateTransitionData(group2, NA, grid_args, perm_trans = TRUE)
        
        # convert to matrix format
        m1 <- convertToMatrix(group1_trans)
        m2 <- convertToMatrix(group2_trans)
        
        # trim rows and cols containing 0 in both matrices
        trimed_mats <- trimMatrixZeros(m1, m2)
        m1 <- trimed_mats[[1]]
        m2 <- trimed_mats[[2]]
        
        # convert to probability matrices
        m1 <- convertToProbMatrix(m1)
        m2 <- convertToProbMatrix(m2)
        
        # get distance 
        distance_result <- HellDistance(m1, m2) 
        
        # store computed distance in vector and return
        distances[i] <- distance_result
        
        # display progress bar
        setWinProgressBar(progress_bar, i, title = paste(round(i / perms * 100, 0), "% processed [Computing Distance]"))
    }
    # close the progress bar
    close(progress_bar)
    
    return(distances)
}

#---------------------------------------------------------------------------------
# FUNCTION:     calculatePvalue()
# INPUT:        vector, vector
# OUTPUT:       void
# DESCRIPTION:  Calculate p-value (% of values > correct/incorrect value)
#               
#---------------------------------------------------------------------------------
calculatePvalue <- function(distance, shuffled_distances)
{
    gtr <- length(shuffled_distances[shuffled_distances > distance])
    pvalue <- gtr / length(shuffled_distances) 
    return(pvalue)
}

#---------------------------------------------------------------------------------
# FUNCTION:     generateDensityPlot()
# INPUT:        double, vector, list
# OUTPUT:       void
# DESCRIPTION:  Output density plot 
#               
#---------------------------------------------------------------------------------
generateDensityPlot <- function(distance_result, shuffled_distances, stimuli_data)
{
    density_plot <- density(shuffled_distances)
    plot(density_plot, type = "n", main = stimuli_data$label, xlab = "Hellinger Distance")
    polygon(density_plot, col = "lightgray", border = "grey")
    rug(shuffled_distances, col = ifelse(shuffled_distances == distance_result, 'blue', 'red'))
    print(abline(v = distance_result, col = "purple"))
    print(density_plot)
}

#---------------------------------------------------------------------------------
# FUNCTION:     outputReport()
# INPUT:        list
# OUTPUT:       void
# DESCRIPTION:  Output results of distance and p-value 
#               
#---------------------------------------------------------------------------------
outputReport <- function(report_args)
{
    cat("\n\nData: ", report_args[["dataname"]], "\n")
    df <- data.frame(JSD = report_args[["HEL"]], p.value = report_args[["p-value"]], Permutations = report_args[["permutations"]])
    print(df) 
    cat("\nGroup 1 (n): ", report_args[["groupsize1"]])
    cat("\nGroup 2 (n): ", report_args[["groupsize2"]], "\n")  
}

#---------------------------------------------------------------------------------
# FUNCTION:     HellDistance()
# INPUT:        matrix, matrix
# OUTPUT:       Summed Hellinger Distance
# DESCRIPTION:  Hellinger Distance is defined as: 1/sqrt(2) * sqrt(sum(suqare(sqrt(pi - squrt(qi))))
#---------------------------------------------------------------------------------
HellDistance <- function(m1,m2)
{
    length_of_matrix <- nrow(m1)
    HPQ <- 0; HJ <- 0
    for (i in 1:length_of_matrix)
    {
        for (j in 1:length_of_matrix)
        {
            HJ <- HJ + ((sqrt(m1[i, j]) - sqrt(m2[i, j])) ^ 2)
        }
        HPQ <- HPQ + ((1 / sqrt(2)) * sqrt(HJ))
        HJ <- 0
    }
    return(HPQ / length_of_matrix) 
}

#---------------------------------------------------------------------------------
# FUNCTION:     convertToPrior()
# INPUT:        matrix
# OUTPUT:       matrix
# DESCRIPTION:  Optimised Bayesian prior (Dirichlet) method
#
#---------------------------------------------------------------------------------
convertToPrior <- function(data)
{
    m <- 0
    matrix_len <- nrow(data) 
    for(i in 1:matrix_len)
    {
        m <- sum(data[i, ] * 10)
        for(j in 1:matrix_len)
            data[i, j] <- ((data[i, j] * 10) + 1) / (m + matrix_len)
    }
    return(data)
}

#---------------------------------------------------------------------------------
# FUNCTION:     main()
# INPUT:        void
# OUTPUT:       void
# DESCRIPTION:  Main function. 
#               Makes all subsequent function calls.     
#---------------------------------------------------------------------------------
main <- function()
{
    permutations <- 2
    num_conditions <- 2
    gridsizes <- c(50, 50, 45, 49, 48, 44, 40, 43) * 2 # eps values from DBSCAN
    
    initialize()
    opened_files <- loadDataFiles()
    conditions1 <- opened_files[["condition1"]]
    conditions2 <- opened_files[["condition2"]]
    
    conditions1 <- as.data.frame(append(conditions1, list("Condition" = 1), after = 1))
    conditions2 <- as.data.frame(append(conditions2, list("Condition" = 2), after = 1))
    merged_data <- rbind(conditions1, conditions2)
    merged_data <- removeRepeatFixation(merged_data)
    
    # get number of participant in each conditions
    num_participants_cond1 <- length(as.vector(conditions1[!duplicated(conditions1$ParticipantName), "ParticipantName"]))
    num_participants_cond2 <- length(as.vector(conditions2[!duplicated(conditions2$ParticipantName), "ParticipantName"]))
    
    # load additional source files with common analysis functions
    common_source_data <- paste0(getwd(), "/gallery-study/gallerymetadata.R")
    source(common_source_data) 
    
    # time the execution of the code
    start_time <- proc.time()
    
    # loop over stimuli
    for(i in 1:length(study_stimuli))
    {
        m1 <- NULL
        m2 <- NULL
        
        # set current stimulus
        stimuli_data <- getStimuliMetaData(study_stimuli[i])
        current_stimuli_data <- merged_data[merged_data$MediaName == stimuli_data$file_name, ]
        
        # calculate dimensions / offsets
        grid_args <- list()
        grid_args[["w"]] <- current_stimuli_data[1, "MediaWidth"] 
        grid_args[["h"]] <- current_stimuli_data[1, "MediaHeight"]
        grid_args[["min_dim"]] <- min(grid_args$w, grid_args$h)
        
        grid_args[["dims"]] <- gridsizes[i]
        grid_args[["images_offsets"]] <- determineImageOffset(grid_args$dims, stimuli_data)
        grid_args[["grid_offsets"]] <- determineGridOffset(grid_args$w, grid_args$h, grid_args$dims)
        
        if(grid_args[["h"]] < grid_args[["w"]])
        {
            grid_args[["num_rows"]] <- floor((grid_args$min_dim - ((grid_args$images_offsets[2] * 2) + (grid_args$grid_offsets[2] * 2))) / grid_args$dims)
            grid_args[["num_cols"]] <- grid_args[["num_rows"]]
        } 
        else 
        {
            grid_args[["num_cols"]] <- floor((grid_args$min_dim - ((grid_args$images_offsets[1] * 2) + (grid_args$grid_offsets[1] * 2))) / grid_args$dims)
            grid_args[["num_rows"]] <- grid_args[["num_cols"]]
        }
        grid_args[["num_cells"]] <- grid_args$num_rows * grid_args$num_cols
        
        # plot painting
        drawGrid(grid_args, stimuli_data)
        
        conditions1_participants <- current_stimuli_data[current_stimuli_data$Condition == 1, ]
        conditions1_participants <- as.vector(conditions1_participants[!duplicated(conditions1_participants$ParticipantName), "ParticipantName"])
        
        # get transition data for both condition groups
        cond_one <- generateTransitionData(current_stimuli_data, 1, grid_args)
        cond_two <- generateTransitionData(current_stimuli_data, 2, grid_args)
        
        # convert to matrix
        m1 <- convertToMatrix(cond_one)
        m2 <- convertToMatrix(cond_two)
        
        # trim cols and rows that are both 0 in both matrices
        trimed_mats <- trimMatrixZeros(m1, m2)
        if(trimed_mats)
        {
            m1 <-  trimed_mats[[1]]
            m2 <-  trimed_mats[[2]]
        }
        
        # convert to probability matrices 
        m1 <- convertToProbMatrix(m1)
        m2 <- convertToProbMatrix(m2)
        
        # work out distance between correct and incorrect groups
        distance_result <- HellDistance(m1, m2)
       
        # carry out permutation test
        shuffled_distances <- randomPermutations(current_stimuli_data, length(conditions1_participants), grid_args, perms = permutations)
        
        # output density plot
        generateDensityPlot(distance_result, shuffled_distances, stimuli_data)
        
        # generate report args
        report_args <- list()
        report_args[["dataname"]] <- stimuli_data$label
        report_args[["HEL"]] <- distance_result
        report_args[["p-value"]] <- calculatePvalue(distance_result, shuffled_distances) 
        report_args[["permutations"]] <- permutations
        report_args[["groupsize1"]] <- num_participants_cond1
        report_args[["groupsize2"]] <- num_participants_cond2
        outputReport(report_args)
        
        #generateLevelPlot(cond_one, cond_two, stimuli_data$label)
    }
    
    # stop timing
    print(proc.time() - start_time)
}

# run main
main()




