##################################################################################
# NAME:         gallery-vis.R
# AUTHOUR:      Alan Davies
# DATE:         09/11/2016
# INSTITUTION:  Interaction Analysis and Modelling Lab (IAM), University of Manchester
# DESCRIPTION:  Visualisation of gallery data, using a combinatorial metric approach
#               and displaying the normalised outputs in a matrix
##################################################################################

#---------------------------------------------------------------------------------
# FUNCTION:     loadDataFiles()
# INPUT:        void
# OUTPUT:       list
# DESCRIPTION:  Returns data for study conditions
#              
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
# FUNCTION:     loadSequenceData()
# INPUT:        void
# OUTPUT:       list
# DESCRIPTION:  Returns all the sequence data files in the directory for
#               both conditions
#---------------------------------------------------------------------------------
loadSequenceData <- function()
{
    data_files <- list()
    sequence_path <- paste0(getwd(), "/MAG/data/sequence-data/")
    
    # set the dir path and access all files
    files_list <- list.files(path = sequence_path, pattern = "csv", all.files = TRUE,
                             full.names = TRUE, recursive = TRUE,
                             ignore.case = TRUE, include.dirs = TRUE)
    
    # store file refrences in list 
    for(i in 1:length(files_list))
        data_files[[basename(files_list[i])]] <- read.csv(files_list[i], header = TRUE, na.strings=c(" ", "NA", "-"))     
    
    return(data_files)
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
# FUNCTION:     getParticipantData()
# INPUT:        data.frame, data.frame
# OUTPUT:       data.frame
# DESCRIPTION:  Extract all records related to the list of participants provided
#               
#---------------------------------------------------------------------------------
getParticipantData <- function(participants, data)
{
    return(data[data$ParticipantName %in% participants, ])
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
    package.args <- c("ggplot2", "lattice", "gplots", "RColorBrewer", "latticeExtra") 
    loadPackages(package.args)
}

#---------------------------------------------------------------------------------
# FUNCTION:     compareSequences()
# INPUT:        data.frame, data.frame
# OUTPUT:       LD
# DESCRIPTION:  Return the Levenshtein distance between the two scanpaths provided
#               after converting them to character sequences
#---------------------------------------------------------------------------------
compareSequences <- function(participant, other)
{
    participant <- toString(lapply(participant$AOI, as.character))
    participant <- stripRedundentChars(participant)
    other <- toString(lapply(other$AOI, as.character))
    other <- stripRedundentChars(other)
    
    return(adist(participant, other))
}

#---------------------------------------------------------------------------------
# FUNCTION:     compareSequences()
# INPUT:        data.frame, data.frame
# OUTPUT:       int
# DESCRIPTION:  The difference between the largest and smallest fixation duration
#               of the compared participants
#---------------------------------------------------------------------------------
compareFixationDurations <- function(participant, other)
{
    participant <- sum(participant$fixationDuration)
    other <- sum(other$fixationDuration)
    
    if(participant != other)
    {
        the_max <- max(participant, other)
        the_min <- min(participant, other)
        diff <- the_max - the_min
    } else {
        diff <- 0
    }
    return(diff)
}

#---------------------------------------------------------------------------------
# FUNCTION:     compareCommonAOI()
# INPUT:        data.frame, data.frame
# OUTPUT:       int
# DESCRIPTION:  Returns the number of shared AOIs that the compared participants
#               have in common
#---------------------------------------------------------------------------------
compareCommonAOI <- function(participant, other) 
{
    # remove participant from other
    participant_name <- unique(participant$participant)
    
    # if they have AOIs return number of shared AOIs
    if(length(other$AOI) > 0 && length(participant$AOI) > 0)
    {
        shared_AOI <- other[other$AOI %in% participant$AOI, 4]
        return(length(unique(as.vector(shared_AOI))))
    } 
    
    return(NA)
}

#---------------------------------------------------------------------------------
# FUNCTION:     stripRedundentChars()
# INPUT:        String
# OUTPUT:       String
# DESCRIPTION:  Remove comma and space from sequence string
#               
#---------------------------------------------------------------------------------
stripRedundentChars <- function(char_string)
{
    char_string <- gsub(",", "", char_string)
    char_string <- gsub(" ", "", char_string)
    return(char_string)    
}

#---------------------------------------------------------------------------------
# FUNCTION:     outputMetricPlots()
# INPUT:        list, data.frame, int, int
# OUTPUT:       void
# DESCRIPTION:  Generate plots for each metric
#               
#---------------------------------------------------------------------------------
outputMetricPlots <-function(args, participants_names, stimuli_num, cond_num)
{
    for(i in 1:length(args))
    {
        data <- args[[i]]
        if(i == 1)
        {
            pal <- "RdBu"
        } else {
            pal <- "OrRd"
        }
       
        title <- paste0("Group ", cond_num, ": ", conditionName(i), "\n", getStimuliMetaData(study_stimuli[stimuli_num])["label"])
        colpal <- colorRampPalette(brewer.pal(11, pal))
        ld = levelplot(data, col.regions = colpal, at = unique(c(seq(from = 0, 
                       to = max(data, na.rm = TRUE), by = 0.1))), scales = list(x = list(rot = 90)), 
                       main = title, xlab = "Participant", ylab = "Participant") 
                        
        print(ld)
    }
}

#---------------------------------------------------------------------------------
# FUNCTION:     conditionName()
# INPUT:        int
# OUTPUT:       void
# DESCRIPTION:  Wrapper function for switch to display appropriate labels for
#               the different metics
#---------------------------------------------------------------------------------
conditionName <- function(i)
{
    switch(i, "1" = "Levenshtein distance", "2" = "Total fixation durations", "3" = "Common AOIs")
}

#---------------------------------------------------------------------------------
# FUNCTION:     outputCombinedMetricPlot()
# INPUT:        list, data.frame, int, int
# OUTPUT:       void
# DESCRIPTION:  Generate plots for combined normalised metrics
#               
#---------------------------------------------------------------------------------
outputCombinedMetricPlot <- function(args, participants, stimuli_num, cond_num)
{
    dims <- length(participants)
    mat_names <- as.vector(participants)
    combined_mat <- matrix(NA, nrow = dims, ncol = dims, dimnames = list(mat_names, mat_names)) 
    
    for(i in 1:length(args))
    {
        data <- args[[i]]
        for(j in 1:nrow(combined_mat))
        {
            for(k in 1:ncol(combined_mat))
            {
                if(i == 1)
                {
                    combined_mat[j, k] <- data[j, k]
                } else {
                    if(!is.na(combined_mat[j, k]))
                    {
                        combined_mat[j, k] <- combined_mat[j, k] + data[j, k]
                    }
                }
            }
        }
    }
    combined_mat <- log10(combined_mat)
    title <- paste0("Group ", cond_num, "\n", getStimuliMetaData(study_stimuli[stimuli_num])["label"])
    colpal <- colorRampPalette(brewer.pal(11, "RdBu"))
    ld = levelplot(combined_mat, col.regions = colpal, at = unique(c(seq(from = 0, 
                   to = max(combined_mat, na.rm = TRUE), by = 0.1))), scales = list(x = list(rot = 90)), 
                   main = title, xlab = "Participant", ylab = "Participant") 
    
    print(ld)
}

#---------------------------------------------------------------------------------
# FUNCTION:     outputSummaryStatsReport()
# INPUT:        list, int, int
# OUTPUT:       void
# DESCRIPTION:  Generate summary statistics report
#               
#---------------------------------------------------------------------------------
outputSummaryStatsReport <- function(args, stimuli_num, cond_num)
{
    stimuli <- toString(getStimuliMetaData(study_stimuli[stimuli_num])["label"])
    cat("\nCondition: ", cond_num, " Stimuli: ", stimuli)
    cat("\nLevenshtein dist: max(", max(args[[1]], na.rm = TRUE), "), min(", min(args[[1]], na.rm = TRUE), 
        "), mean(", mean(args[[1]], na.rm = TRUE), "), sd(", sd(args[[1]], na.rm = TRUE), ")")
    cat("\nTotal fixation duration (difference): max(", max(args[[2]], na.rm = TRUE), "), min(", min(args[[2]], na.rm = TRUE), 
        "), mean(", mean(args[[2]], na.rm = TRUE), "), sd(", sd(args[[2]], na.rm = TRUE), ")")
    cat("\nCommon AOI (number): max(", max(args[[3]], na.rm = TRUE), "), min(", min(args[[3]], na.rm = TRUE), 
        "), mean(", mean(args[[3]], na.rm = TRUE), "), sd(", sd(args[[3]], na.rm = TRUE), ")\n")
}

#---------------------------------------------------------------------------------
# FUNCTION:     combineMetrics()
# INPUT:        list
# OUTPUT:       void
# DESCRIPTION:  Compare one participant against the others based on a number of
#               metrics (levenshtein distance, total fixation difference,
#               and number of common AOIs)
#---------------------------------------------------------------------------------
combineMetrics <- function(args)
{
    metric.args <- list()
    participants <- args[[1]]
    condition <- args[[2]]
    sequence_data <- args[[3]]
    stimuli <- args[[4]]
    stimuli_index <- args[[5]]
    condition_index <- args[[6]]
        
    stimuli_name <- getStimuliMetaData(study_stimuli[stimuli_index])
    stimuli_data <- condition[condition$MediaName == stimuli_name, ]
    stimuli_data <- filterFixationData(stimuli_data)
    
    current_stimuli <- as.data.frame(sequence_data[stimuli])
    dims <- length(participants)
    mat_names <- as.vector(participants)
    
    # create metric matracies for each stimulus 
    ld_mat <- matrix(NA, nrow = dims, ncol = dims, dimnames = list(mat_names, mat_names)) # levenshtein distance
    fd_mat <- matrix(NA, nrow = dims, ncol = dims, dimnames = list(mat_names, mat_names)) # difference in total fixation duration
    ca_mat <- matrix(NA, nrow = dims, ncol = dims, dimnames = list(mat_names, mat_names)) # number of common AOIs

    # compare metrics for one participant against all others and store in matrix cell
    for(i in 1:length(participants))
    {
        names(current_stimuli) <- c("num", "participant", "fixationDuration", "AOI")
        current_participant_sd <- current_stimuli[current_stimuli$participant == participants[i], ]
        if(nrow(current_participant_sd) > 0)
        {
            for(j in 1:length(participants))
            {
                stimuli_data <- stimuli_data[stimuli_data$ParticipantName == participants[j], ]
                other_participant <- current_stimuli[current_stimuli$participant == participants[j], ]
                current_participant <- current_stimuli[current_stimuli$participant == participants[i], ]
                ld_mat[i, j] <- compareSequences(current_participant_sd, other_participant)
                fd_mat[i, j] <- compareFixationDurations(current_participant_sd, other_participant)
                ca_mat[i, j] <- compareCommonAOI(current_participant_sd, other_participant)
            }
        }
    }
    
    # store the matraces in a list
    metric.args[["levenshtein"]] <- ld_mat
    metric.args[["fixationduration"]] <- fd_mat
    metric.args[["commonaoi"]] <- ca_mat
    
    # generate plots and summary statistics
    outputMetricPlots(metric.args, participants, stimuli_index, condition_index)
    outputCombinedMetricPlot(metric.args, participants, stimuli_index, condition_index)
    outputSummaryStatsReport(metric.args, stimuli_index, condition_index)
}

#---------------------------------------------------------------------------------
# FUNCTION:     getTotalFixationDurations()
# INPUT:        data.frame
# OUTPUT:       void
# DESCRIPTION:  Produces data frames showing total fixation durations per
#               painting
#---------------------------------------------------------------------------------
getTotalFixationDurations <- function(condition)
{
    stimuli_id <- NULL
    total_FD <- NULL
    
    for(i in 1:length(study_stimuli))
    {
        stimuli_name <- getStimuliMetaData(study_stimuli[i])
        stimuli_data <- condition[condition$MediaName == stimuli_name, ]
        stimuli_data <- filterFixationData(stimuli_data)
    
        stimuli_id <- c(stimuli_id, stimuli_name[["label"]])
        total_FD <- c(total_FD, sum(stimuli_data$GazeEventDuration))
    }
    df <- data.frame(stimuli = stimuli_id, TFD = total_FD)
    return(df)
}

#---------------------------------------------------------------------------------
# FUNCTION:     getTotalFixationDurations()
# INPUT:        data.frame, data.frame
# OUTPUT:       void
# DESCRIPTION:  Produces plots combining condition 1 and 2 showing total fixation 
#               durations per painting
#---------------------------------------------------------------------------------
plotFixationData <- function(condition1, condition2)
{
    condition1$condition <- c(rep("1", nrow(condition1)))
    condition2$condition <- c(rep("2", nrow(condition2)))
    conditions <- rbind(condition1, condition2)
    
    cond_plot <- ggplot(conditions, aes(factor(stimuli), TFD, fill = condition)) +
        geom_bar(stat = "identity", position = "dodge") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        labs(title = "Total fixation duration", x = "Painting", y = "Total fixation duration (ms)") + 
        scale_fill_discrete(name = "Condition\n") + scale_fill_grey()
    
    print(cond_plot)
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
    combine.args <- list()
    participants <- list()
    initialize()
    
    # get list of participants id's in both conditions
    opened_files <- loadDataFiles()
    participants[["condition1"]] <- getStudyParticipants(opened_files[["condition1"]])
    participants[["condition2"]] <- getStudyParticipants(opened_files[["condition2"]])
  
    # load the sequence data files
    sequence_data <- loadSequenceData()
    
    # load additional source files with common analysis functions
    common_source_data <- paste0(getwd(), "/MAG/gallerymetadata.R")
    source(common_source_data) 
    
    # get the data for both conditions
    condition_one <- getParticipantData(participants[["condition1"]], opened_files[["condition1"]])
    condition_two <- getParticipantData(participants[["condition2"]], opened_files[["condition2"]])
    
    # produce plots for the total fixation duration per painting for each condition
    cond1_FD <- getTotalFixationDurations(opened_files[["condition1"]])
    cond2_FD <- getTotalFixationDurations(opened_files[["condition2"]])
    plotFixationData(cond1_FD, cond2_FD)
    
    # output metric matrices
    for(i in 1:length(participants))
    {
        for(j in 1:length(study_stimuli))
        {
            if(i == 1)
            {
                cond_data <- opened_files[["condition1"]]
                participants_data <-  participants[["condition1"]]
            } else {
                cond_data <- opened_files[["condition2"]]
                participants_data <-  participants[["condition2"]]
            }

            # generate args for metric combination
            combine.args[["participants"]] <- participants_data
            combine.args[["condition"]] <- cond_data
            combine.args[["sequence"]] <- sequence_data
            combine.args[["stimuli"]] <- paste0("cond_", i, "_stimuli_", j, "_sequence_data.csv")
            combine.args[["index"]] <- j
            combine.args[["condition_number"]] <- i
            combineMetrics(combine.args)
        }
    }
}

# run main
main()