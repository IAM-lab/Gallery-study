##################################################################################
# NAME:         DBSCAN.R
# AUTHOUR:      Alan Davies
# DATE:         09/08/2016
# INSTITUTION:  Interaction Analysis and Modelling Lab (IAM), University of Manchester
# DESCRIPTION:  
#               
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
    
    # condition 1
    data_file <- read.csv(cond1_data_path, header = TRUE, na.strings = c(" ", "NA", "-"))
    data_files[["condition1"]] <- data_file
    
    # condition 2
    data_file <- read.csv(cond2_data_path, header = TRUE, na.strings = c(" ", "NA", "-"))
    data_files[["condition2"]] <- data_file
    
    return(data_files)
}

#---------------------------------------------------------------------------------
# FUNCTION:     extractFixationDurationData()
# INPUT:        data.frame
# OUTPUT:       data.frame
# DESCRIPTION:  Selects all the gaze event types that are fixations and
#               remove any items that were not related to a stimulus. Get unique
#               fixations per participant and merge into new data.frame object
#---------------------------------------------------------------------------------
extractFixationDurationData <- function(data)
{
    filtered_data <- data[FALSE, ]
    data <- data[data$GazeEventType == "Fixation", ]
    participants <- getStudyParticipants(data)
    stimuli <- getStimuli(data)
    
    for(i in 1:length(participants))
    {
        participant_data <- data[data$ParticipantName == participants[i], ]
        participant_data <- participant_data[!duplicated(participant_data$FixationIndex), ]
        participant_data <- filterStimuliTransitionFixation(participant_data, stimuli) 
        filtered_data <- rbind(filtered_data, participant_data)
    }
    return(filtered_data[filtered_data$MediaName != "", ])
}

#---------------------------------------------------------------------------------
# FUNCTION:     getStimuli()
# INPUT:        data.frame
# OUTPUT:       vector
# DESCRIPTION:  Returns a vector of unique stimuli names
#               
#---------------------------------------------------------------------------------
getStimuli <- function(data)
{
    data <- data[!duplicated(data$MediaName), ]
    data <- data[data$MediaName != "", ]
    return(as.vector(data$MediaName))
}

#---------------------------------------------------------------------------------
# FUNCTION:     extractStimulus()
# INPUT:        data.frame, String
# OUTPUT:       data.frame
# DESCRIPTION:  Retunrs data for a specific stimulus 
#               
#---------------------------------------------------------------------------------
extractStimulus <- function(data, stimuli_name)
{
    return(data[data$MediaName == stimuli_name, ])
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
    package.args <- c("dbscan", "png")
    loadPackages(package.args)
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
# FUNCTION:     filterStimuliTransitionFixation()
# INPUT:        data.frame
# OUTPUT:       data.frame
# DESCRIPTION:  For each participant remove the first fixation for all subsequent
#               stimului except the first one and recombine. This is done to remove
#               the first fixation from a new stimulus as it is still related to where
#               the participant was looking on the previous stimulus
#---------------------------------------------------------------------------------
filterStimuliTransitionFixation <- function(data, stimuli)
{
    filtered_data <- data[FALSE, ]
    for(i in 1:length(stimuli))
    {
        current_stimuli <- extractStimulus(data, stimuli[i])
        if(i > 1) current_stimuli <- current_stimuli[-1, ]
        filtered_data <- rbind(filtered_data, current_stimuli)
    }
    return(filtered_data)
}

#---------------------------------------------------------------------------------
# FUNCTION:     swapStimuliNames()
# INPUT:        void
# OUTPUT:       vector
# DESCRIPTION:  Swap stimuli label for more meaninful painting name
#               using a substitution vector 
#---------------------------------------------------------------------------------
swapStimuliNames <- function()
{
    translation_vector <- c("1917.170med_resized.jpg.png" = "Rhyl Sands", 
                            "1922.5med_resized.jpg" = "Flask walk, hampstead",
                            "1934.2med.jpg.png" = "Self-portrait",
                            "1934.394med_resized.jpg.png" = "When the west with evening glows",
                            "1964.287med.jpg.png" = "14.6.1964",
                            "1968.173med.jpg.png" = "Woman and Suspended Man",
                            "1976.79med_resized.jpg.png" = "Sir Gregory Page-Turner",
                            "1995.184med.jpg.png" = "Release")
    
    return(translation_vector)
}

#---------------------------------------------------------------------------------
# FUNCTION:     runDBSCAN()
# INPUT:        int, data.frame, String, list
# OUTPUT:       void
# DESCRIPTION:  Runs and displays results of DBSCAN algorithm, overlaying clusters
#               ontop of the stimuli images (paintings)
#---------------------------------------------------------------------------------
runDBSCAN <- function(condition, data, stimuli_name, args)
{
    trans_vec <- swapStimuliNames()
    title_str <- paste0("Condition ", condition, ": ", trans_vec[stimuli_name])
    
    # load background image filtering out filetype string data
    stimuli_image_str <- gsub(".jpg", "", stimuli_name)
    stimuli_image_str <- gsub(".png", "", stimuli_image_str)
    bkg_img <- readPNG(paste0(getwd(), "/MAG/stimuli/", stimuli_image_str, ".png"))

    # run DBSCAN and output results
    matrix_data <- as.matrix(na.omit(data))
    db <- dbscan(matrix_data, eps = args[["eps"]], minPts = args[["minPts"]])
    cat("\n", title_str, "\n")
    print(db)

    # output plots and add bkg image to them
    print(pairs(matrix_data, col = db$cluster + 1L))
    print(plot(matrix_data, col = db$cluster + 1L, main = title_str))  #res$cluster, main = title_str))
    limits <- par()
    rasterImage(bkg_img, limits$usr[1], limits$usr[3], limits$usr[2], limits$usr[4])
    print(grid())
    
    # print raw fixation plot
    print(plot(matrix_data, col = "blue", main = title_str))
    rasterImage(bkg_img, limits$usr[1], limits$usr[3], limits$usr[2], limits$usr[4])
    print(grid())
}

#---------------------------------------------------------------------------------
# FUNCTION:     runOPTICS()
# INPUT:        int, data.frame, String, list
# OUTPUT:       void
# DESCRIPTION:  Runs and displays results of OPTICS algorithm
#           
#---------------------------------------------------------------------------------
runOPTICS <- function(condition, data, stimuli_name, args)
{
    trans_vec <- swapStimuliNames()
    title_str <- paste0("Condition ", condition, ": ", trans_vec[stimuli_name])
    
    matrix_data <- as.matrix(na.omit(data))
    opt <- optics(matrix_data, eps = args[["eps"]], minPts = args[["minPts"]], xi = args[["xi"]])
    cat("\n", title_str, "\n")
    print(opt)
    print(plot(opt))
}

#---------------------------------------------------------------------------------
# FUNCTION:     generatekNNDistPlots()
# INPUT:        data.frame, list
# OUTPUT:       void
# DESCRIPTION:  Produces kNN plot to determine optimal eps value by visualising
#               "knee" in plot curve
#---------------------------------------------------------------------------------
generatekNNDistPlots <- function(data, args)
{
    # plot kNN to determine optimal eps value
    matrix_data <- as.matrix(na.omit(data))
    print(kNNdistplot(matrix_data, k = args[["minPts"]]))
}

#---------------------------------------------------------------------------------
# FUNCTION:     main()
# INPUT:        void
# OUTPUT:       void
# DESCRIPTION:  Main function. 
#               Makes all subsequent function calls
#---------------------------------------------------------------------------------
main <- function()
{
    initialize()
    data <- loadDataFiles()
    stimuli <- getStimuli(data[["condition1"]])
    cond_str <- "condition"
    stimulus_data <- list()
    algorithm_args <- list()
    k <- 1
    
    # set the properties of the algorithms 
    algorithm_args[["eps"]] <- 10               # 15
    algorithm_args[["minPts"]] <- 4             # 4
    algorithm_args[["xi"]] <- 0.05              # 0.05
    
    # loop over conditions 
    for(i in 1:length(data))
    {
        fixation_points_df <- NULL
        fixation_points <- NULL
        participant_data <- list()
        
        # get condition and subset
        cond_data <- data[[paste0(cond_str, i)]]
        
        # subset the fixation duration data
        cond_data <- extractFixationDurationData(cond_data)
        
        # get unique participants and name list elements the same
        participants <- getStudyParticipants(cond_data)
        participant_data <- setNames(vector("list", length(participants)), participants)
        
        # loop over stimuli (paintings)
        for(j in 1:length(stimuli))
        {
            # get the current stimuli
            current_stimuli <- extractStimulus(cond_data, stimuli[j])
            current_stimuli <- current_stimuli[!duplicated(current_stimuli), ]
            
            stimulus_data[[k]] <- current_stimuli
            k <- k + 1
            
            for(l in 1:length(participants))
            {
                # aggregate fixation data accross all stimuli
                participant_FD <- current_stimuli[current_stimuli$ParticipantName == participants[l], ]
                if(length(participant_FD$GazeEventDuration) > 0)
                {
                    if(!is.na(participant_FD$FixationPointX) && !is.na(participant_FD$FixationPointY))
                    {
                        # get all x and y fixation points for each participant
                        x <- participant_FD$FixationPointX
                        y <- participant_FD$FixationPointY
                    
                        fixation_points <- cbind(x, y)
                    }
                }
                # put them together in a single data frame
                fixation_points_df <- rbind(fixation_points_df, fixation_points) 
            }
            # find optimal eps values
            generatekNNDistPlots(fixation_points_df, algorithm_args)
            
            # run DBSCAN and OPTICS algorithms
            runDBSCAN(i, fixation_points_df, stimuli[j], algorithm_args)
            #runOPTICS(i, fixation_points_df, stimuli[j], algorithm_args)
        }
    }
}

# run main
main()
