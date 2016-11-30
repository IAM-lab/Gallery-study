##################################################################################
# NAME:         MCA grid plus.R
# AUTHOUR:      Alan Davies
# DATE:         06/09/2016
# INSTITUTION:  Interaction Analysis and Modelling Lab (IAM), University of Manchester
# DESCRIPTION:  Version of Markov chain transition analysis using a grid system
#               based on DBSCAN cluster size.
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
    cond1_data_path <- paste0(getwd(), "/MAG/data/condition1.csv")
    cond2_data_path <- paste0(getwd(), "/MAG/data/condition2.csv")
    
    # conditions
    data_file <- read.csv(cond1_data_path, header = TRUE, na.strings = c(" ", "NA", "-"))
    data_files[["condition1"]] <- data_file
    data_file <- read.csv(cond2_data_path, header = TRUE, na.strings = c(" ", "NA", "-"))
    data_files[["condition2"]] <- data_file
    
    return(data_files)
}

#---------------------------------------------------------------------------------
# FUNCTION:     buildTransitionMatrix()
# INPUT:        data.frame, int, int, BOOL
# OUTPUT:       matrix
# DESCRIPTION:  Converts transition quantities into probability scores
#               
#---------------------------------------------------------------------------------
buildTransitionMatrix <- function(data, rows, cols, bayesian_prior = FALSE) 
{     
    # create empty matrix of the correct size
    dims <- rows * cols
    mat <- matrix(0, nrow = dims, ncol = dims)
    names <- seq(from = 1, to = nrow(mat), by = 1)
    dimnames(mat) <- list(names, names)
    
    # calculate transition table and convert to data frame
    msm_matrix <- statetable.msm(data$AOI, subject = 1)    
    df_matrix <- as.data.frame(msm_matrix)
    
    # reshape into a standard 2D matrix
    matrix <- as.matrix(xtabs(Freq ~ from + to, df_matrix))
    row_names <- seq(from = 1, to = nrow(matrix), by = 1)    
    rownames(matrix) <- row_names
    
    # convert from frequency to probability values
    prob_matrix <- matrix / rowSums(matrix)
    m <- as.data.frame(prob_matrix)
    
    for(k in 1:length(prob_matrix))
    {
        from <- m[k, 1]
        to <- m[k, 2]
        freq <- m[k, 3]
        
        for(i in 1:dims)
        {
            for(j in 1:dims)
            {
                if((i == from) && (j == to) && (freq != 0))
                {
                    mat[i, j] <- freq   
                }
            }
        }
    }

    # if Bayesian prior, then add 1 to each transition
    if(bayesian_prior)
    {
        mat + 1
    }
    return(mat)
}

#---------------------------------------------------------------------------------
# FUNCTION:     calculateDistance()
# INPUT:        matrix, matrix
# OUTPUT:       eclidean distance
# DESCRIPTION:  Computes Jensen Shannon divergence/distance between two matrices. 
#               
#---------------------------------------------------------------------------------
calculateDistance <- function(m1, m2)
{
    results <- list()
    m <- generateM(m1, m2)
    results[["JSD"]] <- (0.5 * KullbackLeiblerDivergence(m1, m)) + (0.5 * KullbackLeiblerDivergence(m2, m))
    results[["JSDist"]] <- sqrt(results[["JSD"]]) / log(2, base = exp(1))
    
    return(results)
}

#---------------------------------------------------------------------------------
# FUNCTION:     KullbackLeiblerDivergence()
# INPUT:        matrix, matrix
# OUTPUT:       Summed KL divergence
# DESCRIPTION:  Computes the KL divergence between two matrices assumed to
#               have the same dimensions.
#---------------------------------------------------------------------------------
KullbackLeiblerDivergence <- function(m1, m2)
{
    dist <- 0
    summed_state <- 0
    matrix_length <- min(nrow(m1), nrow(m2))
    
    for(i in 1:matrix_length)
    {
        for(j in 1:matrix_length)
        {
            # calculate KLD per row
            if((m1[i, j] != 0) && (m2[i, j] != 0))
                dist <- dist + m1[i, j] * log(m1[i, j] / m2[i, j], base = exp(1))
        }
        # compute summed KLD
        summed_state <- summed_state + dist
        dist <- 0
    }
    return(summed_state / matrix_length)           
}

#---------------------------------------------------------------------------------
# FUNCTION:     generateM()
# INPUT:        matrix
# OUTPUT:       matrix
# DESCRIPTION:  m = 1/2(p + q) Part of the JSD
#               
#---------------------------------------------------------------------------------
generateM <- function(m1, m2)
{
    matrix_length <- nrow(m1)
    m <- matrix(0, nrow = matrix_length, ncol = matrix_length)
    
    for(i in 1:matrix_length)
        for(j in 1:matrix_length)
            m[i, j] <- 0.5 * (m1[i, j] + m2[i, j])
    
    return(m)
}

#---------------------------------------------------------------------------------
# FUNCTION:     calculatePvalue()
# INPUT:        vector, vector
# OUTPUT:       void
# DESCRIPTION:  Calculate p-value (% of values > correct/incorrect value)
#               
#---------------------------------------------------------------------------------
calculatePvalue <- function(correct_and_incorrect, shuffled_distances)
{
    gtr <- length(shuffled_distances[shuffled_distances > correct_and_incorrect])
    pvalue <- gtr / length(shuffled_distances) 
    return(pvalue)
}

#---------------------------------------------------------------------------------
# FUNCTION:     generateReport()
# INPUT:        list
# OUTPUT:       void
# DESCRIPTION:  Outputs report summary of the stimuli and the both
#               distance measures with p-values to screen and file.
#---------------------------------------------------------------------------------
generateReport <- function(report_args)
{
    # output to file
    output_path <- paste0(getwd(), "/output.txt")
    sink(output_path, append = TRUE)
    outputReport(report_args)
    sink()
    
    # output to screen
    outputReport(report_args)
}

#---------------------------------------------------------------------------------
# FUNCTION:     outputReport()
# INPUT:        list
# OUTPUT:       void
# DESCRIPTION:  Generates report summary of the stimuli and the both
#               distance measures with p-values.
#---------------------------------------------------------------------------------
outputReport <- function(report_args)
{
    cat("\n\nData: ", report_args[["dataname"]], "\n")
    JSD <- report_args[["JSD"]]  
    JSD_pval <- report_args[["JSDpvalue"]]
    JSdist <- report_args[["JSDist"]]
    JSdist_pval <- report_args[["JSDistpvalue"]]
    perm <- report_args[["permutations"]]
    
    df <- data.frame(JSD = JSD, p.value.JSD = JSD_pval, JSdist = JSdist, p.value.JSDist = JSdist_pval, xP = perm)
    print(df) 
    cat("\nGroup 1 (n): ", report_args[["groupsize1"]])
    cat("\nGroup 2 (n): ", report_args[["groupsize2"]])  
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
    package.args <- c("msm", "png")
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
    bkg_img <- readPNG(paste0(getwd(), "/MAG/stimuli/", stimuli_image_str, ".png"), info = TRUE)
    plot.new()
    limits <- par()
    rasterImage(bkg_img, limits$usr[1], limits$usr[3], limits$usr[2], limits$usr[4])
    
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
# FUNCTION:     detectCellHits()
# INPUT:        int, int, list
# OUTPUT:       int
# DESCRIPTION:  Detect the number of hits for the current stimulus in the
#               designated grid cell area
#---------------------------------------------------------------------------------
detectCellHits <- function(cell_x, cell_y, cell_len, data)
{
    fixation_points <- NULL
    
    data <- data[1, c("FixationPointX", "FixationPointY")]
    
    # extract data
    fixation_x <- data[1, 1]
    fixation_y <- data[1, 2]
    
    # return no hit for any missing data
    if(is.na(fixation_x) || is.na(fixation_y))
    {
        return(FALSE)
    }
    
    # check present inside cell and increments hits and total FD
    if((fixation_x > cell_x && fixation_x < (cell_x + cell_len)) && 
       (fixation_y > cell_y && fixation_y < (cell_y + cell_len)))
    {
        return(TRUE)
    }
    return(FALSE)
}

#---------------------------------------------------------------------------------
# FUNCTION:     filterFixationData()
# INPUT:        data.frame
# OUTPUT:       data.frame
# DESCRIPTION:  Selects all the gaze event types that are fixations and
#               remove any items that were not related to a stimulus. Get unique
#               fixations per participant and merge into new data.frame object
#---------------------------------------------------------------------------------
filterFixationData <- function(data)
{
    fixation_data <- data[data$GazeEventType == "Fixation", ]
    filtered_data <- fixation_data[!duplicated(fixation_data$FixationIndex), ]
    filtered_data <- na.omit(filtered_data[ ,c("FixationPointX", "FixationPointY", "GazeEventDuration")])
    
    return(filtered_data)
}

#---------------------------------------------------------------------------------
# FUNCTION:     getStudyParticipants()
# INPUT:        data.frame
# OUTPUT:       vector
# DESCRIPTION:  Extract the unique participant names from the data frame
#               
#---------------------------------------------------------------------------------
getStudyParticipants <- function(data)
{
    return(as.vector(data[!duplicated(data$ParticipantName), 1]))
}

#---------------------------------------------------------------------------------
# FUNCTION:     initializeGrid()
# INPUT:        data.frame, int, list
# OUTPUT:       ?
# DESCRIPTION:  Prepares any arguments for the grid generation function and loads
#               all necessary elements into list elements
#---------------------------------------------------------------------------------
initializeGrid <- function(data, grid_size, stimuli)
{
    grid_args <- list()
    grid_args[["grid_size"]] <- grid_size
    grid_args[["grid_data"]] <- data
    grid_args[["stimuli_data"]] <- stimuli
    grid_args[["boundary_detection"]] <- TRUE
    grid_args[["best_fit"]] <- FALSE
    
    return(grid_args)
}

#---------------------------------------------------------------------------------
# FUNCTION:     generateUOIGrid()
# INPUT:        int, int
# OUTPUT:       list
# DESCRIPTION:  Generates an empty grid of the containg id, and location of
#               cells with 0 hits
#---------------------------------------------------------------------------------
generateUOIGrid <- function(num_rows, num_cols, images_offsets, grid_offsets, dims)
{
    UOI <- list()
    id_count <- 0
    
    if(num_cols < 1) num_cols <- 1
    
    for(i in num_rows:1)
    {
        UOI[[i]] <- list()
        
        for(j in 1:num_cols)
        {
            cell_x <- images_offsets[1] + grid_offsets[1] + (dims * (j - 1))
            cell_y <- images_offsets[2] + grid_offsets[2] + (dims * (i - 1))
            id <- id_count
        
            UOI[[i]][[j]] <- c(id, cell_x, cell_y, 0, 0)
            id_count <- id_count + 1
        }
    }
    return(UOI)
}

#---------------------------------------------------------------------------------
# FUNCTION:     randomPermutations()
# INPUT:        data.frame, int, int, int
# OUTPUT:       vector
# DESCRIPTION:  Generate 2 random sub groups the same size as the original groups
#               and compare distances for the selected amout of permutations.
#---------------------------------------------------------------------------------
randomPermutations <- function(group1, group2, group1_size, row_len, col_len, stimuli_index, perms = 10)
{
    progress <- 0
    distance <- NULL
    group1_participants <- list()
    group2_participants <- list()
    shuffled_group1_participants <- list()
    shuffled_group2_participants <- list()
    shuffled_group1 <- NULL 
    shuffled_group2 <- NULL
    
    # create progress bar
    #progress_bar <- winProgressBar(title = "Computing Distances", min = 0, max = perms, width = 300)
    #getWinProgressBar(progress_bar)
    
    for(j in 1:perms)
    {
        # get the current stimuli for each group
        group1_stimuli <- group1[[stimuli_index]]
        group2_stimuli <- group2[[stimuli_index]]
        combined_stimuli <- rbind(as.data.frame(group1_stimuli), as.data.frame(group2_stimuli))
        
        # extract the participants for each group and combine in a single list
        group1_participants <- as.list(group1_stimuli[!duplicated(group1_stimuli$participant), 1])
        group2_participants <- as.list(group2_stimuli[!duplicated(group2_stimuli$participant), 1])
        participants <- append(group1_participants, group2_participants)
        
        for(i in 1:length(group1_participants))
        {
            # loop over group 1 and get a random element from participants and 'pop' it
            rand_ind <- sample(1:length(participants), 1)
            shuffled_group1_participants[[i]] <- participants[[rand_ind]]
            participants[[rand_ind]] <- NULL
        }
        # put whatever is left into the second group
        shuffled_group2_participants <- participants     
        shuffled_group1_participants <- as.vector(unlist(shuffled_group1_participants))
        shuffled_group2_participants <- as.vector(unlist(shuffled_group2_participants))
        
        # use the participant groups to extract their data from the data frames in the new groups
        for(i in 1:length(shuffled_group1_participants))
        {
            df <- combined_stimuli[combined_stimuli$participant == shuffled_group1_participants[i], ]
            shuffled_group1 <- rbind(shuffled_group1, df)
        }
        for(i in 1:length(shuffled_group2_participants))
        {
            df <- combined_stimuli[combined_stimuli$participant == shuffled_group2_participants[i], ]
            shuffled_group2 <- rbind(shuffled_group2, df)
        }
       
        m_1 <- buildTransitionMatrix(shuffled_group1, row_len, col_len, bayesian_prior = TRUE)
        m_2 <- buildTransitionMatrix(shuffled_group2, row_len, col_len, bayesian_prior = TRUE)
        rownames(m_1) <- NULL
        rownames(m_2) <- NULL
        
        # get distance 
        distance_results <- calculateDistance(m_1, m_2)
        
        # store computed distance in vector and return
        distance <- c(distance, distance_results[['JSDist']])
        
        # display progress bar
        #setWinProgressBar(progress_bar, j, title = paste(round(j / perms * 100, 0), "% processed [Computing Distance]"))
    }    
    # close the progress bar
    #close(progress_bar)
    
    return(signif(distance, 3))
}

#---------------------------------------------------------------------------------
# FUNCTION:     runMarkovtransitionAnalysis()
# INPUT:        list, list
# OUTPUT:       void
# DESCRIPTION:   
#                    
#---------------------------------------------------------------------------------
runMarkovtransitionAnalysis <- function(group1, group2, grid1, grid2)
{
    for(i in 1:length(study_stimuli))
    {
        m_1 <- NULL
        m_2 <- NULL
        
        stimuli_data <- getStimuliMetaData(study_stimuli[i])
        current_stimuli <- stimuli_data$ref_name
        combined_data <- rbind(group1, group2)
        
        row_len <- length(grid1[[i]])
        col_len <- length(grid1[[i]][[1]])
        
        # build transition matrices
        m_1 <- buildTransitionMatrix(group1[[i]], length(grid1[[i]]), length(grid1[[i]][[1]]), bayesian_prior = TRUE)
        m_2 <- buildTransitionMatrix(group2[[i]], length(grid2[[i]]), length(grid2[[i]][[1]]), bayesian_prior = TRUE)
        rownames(m_1) <- NULL
        rownames(m_2) <- NULL
        
        distance_results <- calculateDistance(m_1, m_2)
        cond1_and_cond2_results <- distance_results[['JSDist']]
        JSD <- distance_results[['JSD']]
        ylabel <- "Jensen-Shannon Distance"
        
        # repeat and add initial data
        shuffled_distances <- randomPermutations(group1, group2, length(group1), row_len, col_len, i, perms = permutations)
        shuffled_distances <- append(shuffled_distances, cond1_and_cond2_results)
        
        # generate density plot
        density_plot <- density(shuffled_distances)
        plot(density_plot, type = "n", main = paste0(stimuli_data$label, " (",ylabel,")"), xlab = "Jensen-Shannon Distance")
        polygon(density_plot, col = "lightgray", border = "grey")
        rug(shuffled_distances, col = ifelse(shuffled_distances == cond1_and_cond2_results, 'blue', 'red'))
        abline(v = cond1_and_cond2_results, col = "purple")
        
        # write plot to file
        dev.copy(png, paste0(getwd(), "/outputgraphs", stimuli_data$label, " (",ylabel,") density.png"))
        dev.off()
        
        # create report args
        report_args <- list()
        report_args[["dataname"]] <- stimuli_data$label
        report_args[["JSD"]] <- JSD
        report_args[["JSDist"]] <- cond1_and_cond2_results
        report_args[["JSDistpvalue"]] <- calculatePvalue(cond1_and_cond2_results, shuffled_distances) 
        report_args[["JSDpvalue"]] <- calculatePvalue(JSD, shuffled_distances)
        report_args[["permutations"]] <- permutations
        report_args[["groupsize1"]] <- 22
        report_args[["groupsize2"]] <- 22
        
        # generate report summary
        generateReport(report_args)
    }
}

#---------------------------------------------------------------------------------
# FUNCTION:     generateHeatMaps()
# INPUT:        list, list
# OUTPUT:       void
# DESCRIPTION:  Wrapper to call heatmap generation function for all
#               instances
#---------------------------------------------------------------------------------
showHeatMaps <- function(grid1, grid2)
{
    generateHeatMaps(grid1, "duration", "1")
    generateHeatMaps(grid2, "duration", "2")
    #generateHeatMaps(grid1, "count", "1")
    #generateHeatMaps(grid2, "count", "2")
}

#---------------------------------------------------------------------------------
# FUNCTION:     generateHeatMaps()
# INPUT:        list, list
# OUTPUT:       void
# DESCRIPTION:  Generate heat maps using fixation data for either the count or
#               duration. 
#---------------------------------------------------------------------------------
generateHeatMaps <- function(grid, type = "duration", condition)
{
    metric <- ifelse(type == "duration", 4, 5)
    graph_title <- ifelse(type == "duration", " \n(Total fixation duration)\n", " (Fixation count)")
    my_palette <- colorRampPalette(c("white", "orange", "red"))(n = 1000)
    
    for(i in 1:length(grid))
    {
        num_rows <- length(grid[[i]])
        mat <- matrix(0, nrow = num_rows, ncol = length(grid[[i]][[1]]))

        for(m in 1:num_rows)
        {
            num_cols <- length(grid[[i]][[m]])
            for(n in 1:num_cols)
            {
                mat[m, n] <- grid[[i]][[m]][[n]][[metric]]
            }
        }
        # generate heatmap
        title_text <- paste0("Stimuli ", i, " Condtion ", condition, graph_title)
        heatmap.2(x = mat, cellnote = mat, col = my_palette, main = title_text, symm = TRUE, dendrogram = "none", 
                  Rowv = FALSE, trace = "none", density.info = "none", notecol = "black", 
                  sepwidth = c(0.01, 0.01), sepcolor = "black", colsep = 1:ncol(mat), rowsep = 1:nrow(mat), srtCol = 0)
    }
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
    participants <- list()
    gridsizes <- c(250, 250, 250, 250, 250, 250, 250, 250) 
    UOI <- list()
    grids_cond1 <- list()
    grids_cond2 <- list()
    stimuli_sequences_cond1 <- list()
    stimuli_sequences_cond2 <- list()
    id_count <- 0
    permutations <<- 10000
    
    initialize()
    
    opened_files <- loadDataFiles()
    participants[["condition1"]] <- getStudyParticipants(opened_files[["condition1"]])
    participants[["condition2"]] <- getStudyParticipants(opened_files[["condition2"]])
    
    # load additional source files with common analysis functions
    common_source_data <- paste0(getwd(), "/MAG/gallerymetadata.R")
    source(common_source_data) 
    
    # time the execution of the code
    start_time <- proc.time()
    
    # loop over conditions, participants and stimuli
    for(i in 1:length(participants))
    {
        current_data <- opened_files[[i]]
        
        for(j in 1:length(study_stimuli))
        {
            UOI_grid <- NULL
            AOI_vector <- NULL
            durations_vector <- NULL
            participant_vector <- NULL
            sequence_data <- NULL
            
            # set current stimulus
            stimuli_data <- getStimuliMetaData(study_stimuli[j])
            current_stimuli_data <- current_data[current_data$MediaName == stimuli_data$file_name, ]
            current_stimuli_data <- current_stimuli_data[current_stimuli_data$GazeEventType == "Fixation", ]
            
            # create UOI's with grid
            grid_args <- initializeGrid(current_stimuli_data, gridsizes[j], stimuli_data) # gridsizes in loop
            
            # calculate dimensions / offsets
            w <- current_stimuli_data[1, "MediaWidth"] 
            h <- current_stimuli_data[1, "MediaHeight"]
            dims <- gridsizes[j]
            images_offsets <- determineImageOffset(dims, stimuli_data)
            grid_offsets <- determineGridOffset(w, h, dims)
            num_cols <- floor((w - (images_offsets[1] * 2)) / dims)
            num_rows <- floor((h - (images_offsets[2] * 2)) / dims)
            if(num_cols < 1) num_cols <- 1
            
            # generate grid
            UOI_grid <- generateUOIGrid(num_rows, num_cols, images_offsets, grid_offsets, dims)
            
            for(k in 1:length(participants[[i]]))
            {
                # subset for participant and remove duplicate fixation indexes 
                current_participant <- current_stimuli_data[current_stimuli_data$ParticipantName == participants[[i]][k], ]
                current_participant <- current_participant[!duplicated(current_participant$FixationIndex), ]
                
                for(l in 1:nrow(current_participant))
                {
                    for(m in num_rows:1)
                    {
                        for(n in 1:num_cols)
                        {
                            # calculate relative position of cell
                            cell_x <- images_offsets[1] + grid_offsets[1] + (dims * (n - 1))
                            cell_y <- images_offsets[2] + grid_offsets[2] + (dims * (m - 1))
                            
                            # detect hit inside cell
                            hit <- detectCellHits(cell_x, cell_y, dims, current_participant[l, ])
                            if(hit)
                            {
                                # increment hits and append AOI sequence data
                                durations_vector <- c(durations_vector, current_participant[l, "GazeEventDuration"])
                                duration <- current_participant[l, "GazeEventDuration"]
                                participant_vector <- c(participant_vector, participants[[i]][k])
                                AOI_vector <- c(AOI_vector, UOI_grid[[m]][[n]][1])
                                
                                UOI_grid[[m]][[n]][4] <- UOI_grid[[m]][[n]][4] + duration
                                UOI_grid[[m]][[n]][5] <- UOI_grid[[m]][[n]][5] + 1
                            }
                        }
                    }
                }
            }
            # store each stimulus grid / sequence data df in a list
            sequence_data <- data.frame(participant = participant_vector, fixationDuration = durations_vector, AOI = AOI_vector)
            #write.csv(sequence_data, paste0("cond_", i, "_stimuli_", j, "_sequence_data.csv"))
            if(i == 1)
            {
                stimuli_sequences_cond1[[j]] <- sequence_data
                grids_cond1[[j]] <- UOI_grid
            } else {
                stimuli_sequences_cond2[[j]] <- sequence_data
                grids_cond2[[j]] <- UOI_grid
            }
        }
    }
    # run the Markov analysis
    runMarkovtransitionAnalysis(stimuli_sequences_cond1, stimuli_sequences_cond2, grids_cond1, grids_cond2)
    showHeatMaps(grids_cond1, grids_cond2)
    
    # stop timing
    print(proc.time() - start_time)
}

# run main
main()
