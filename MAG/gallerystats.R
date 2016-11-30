##################################################################################
# NAME:         gallerystats.R
# AUTHOUR:      Alan Davies
# DATE:         24/06/2016
# INSTITUTION:  Interaction Analysis and Modelling Lab (IAM), University of Manchester
# DESCRIPTION:  Standard statistical analysis of Gallery data
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
# FUNCTION:     extractQuestionResponse()
# INPUT:        data.frame
# OUTPUT:       list
# DESCRIPTION:  Gets the individual participants and then their responses
#               to the question posed as well as gender values 
#---------------------------------------------------------------------------------
extractQuestionResponse <- function(data)
{
    response_vals <- list()
    data <- data[!duplicated(data$ParticipantName), ]
    participants <- data$ParticipantName
    gender_only <- gsub("[0-9]", "", participants)
    gender_only <- substring(gender_only, 2)
    
    # store results in list
    response_vals[["question"]] <- as.vector(data$X.Question.Value)
    response_vals[["gender"]] <- as.vector(gender_only)
    
    return(response_vals)
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
# FUNCTION:     normalityTests()
# INPUT:        data.frame, int, String
# OUTPUT:       void
# DESCRIPTION:  Produces histogram and QQ plots as well as some data description
#               functions to help with assessment of normality of data
#---------------------------------------------------------------------------------
normalityTests <- function(data, condition_num, stimuli_name)
{
    # return if no rows in the data
    if(nrow(data) < 1) return()
    
    # create histogram with normality plot overlaid
    norm_plot <- ggplot(data, aes(GazeEventDuration)) + 
    geom_histogram(aes(y = ..density..), colour = "black", fill = "white") +
    labs(x = "Fixation duration", y = "Density") + 
    ggtitle(paste0("Fixation duration condition ", condition_num, " (", stimuli_name, ")")) + 
    stat_function(fun = dnorm, args = list(mean = mean(data$GazeEventDuration, na.rm = TRUE), 
                                           sd = sd(data$GazeEventDuration, na.rm = TRUE)), 
                                           colour = "red", size = 1)
    print(norm_plot)
    
    # create quantile quantile (QQ plots)
    qq_plot <- qplot(sample = data$GazeEventDuration, stat = "qq", 
                     xlab = "Theoretical Quantiles", ylab = "Sample Quantiles",
                     main = paste0("Fixation duration condition ", condition_num, " (", stimuli_name, ")"))
    print(qq_plot)
    
    # run some tests to look at skew and kurtosis
    #print(describe(data$GazeEventDuration))
    #print(stat.desc(data$GazeEventDuration, basic = FALSE, norm = TRUE))
}

#---------------------------------------------------------------------------------
# FUNCTION:     boxPlotsPerStimulus()
# INPUT:        data.frame, int
# OUTPUT:       void
# DESCRIPTION:  Displays a coloured box plot of all the conditions
#               
#---------------------------------------------------------------------------------
boxPlotsPerStimulus <- function(data, condition_num)
{
    # TODO: Add an if for FD/FC
    
    # swap over the stimuli names for more meaningful painting names
    trans_vec <- swapStimuliNames()
   
    # build and show box plot
    box_plot <- ggplot(data, aes(factor(MediaName), GazeEventDuration, fill = MediaName)) +
        ylim(-200, 2000) + geom_boxplot() +
        labs(title = paste0("Fixation duration (all paintings) condition ", condition_num), x = "Painting", y = "Fixation duration (ms)") +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5), legend.position = "none") + 
        scale_x_discrete(labels = trans_vec[data$MediaName])
    
    print(box_plot)
}

#---------------------------------------------------------------------------------
# FUNCTION:     displayViolinPlots()
# INPUT:        data.frame, int
# OUTPUT:       void
# DESCRIPTION:  Displays a violin plot of all the conditions
#               
#---------------------------------------------------------------------------------
displayViolinPlots <- function(data, condition_num)
{
    # swap over the stimuli names for more meaningful painting names
    trans_vec <- swapStimuliNames()
    
    # build and show box plot
    violin_plot <- ggplot(data, aes(factor(MediaName), GazeEventDuration, fill = MediaName)) + 
        ylim(-200, 2000) + 
        geom_violin() +
        labs(title = paste0("Fixation duration (all paintings) condition ", condition_num), x = "Painting", y = "Fixation duration (ms)") +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5), legend.position = "none") + 
        scale_x_discrete(labels = trans_vec[data$MediaName])
    
    print(violin_plot)
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
    package.args <- c("ggplot2", "psych", "pastecs", "lawstat", "plyr")
    loadPackages(package.args)
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
# FUNCTION:     calculateGroupDifferences()
# INPUT:        list, list, String
# OUTPUT:       void
# DESCRIPTION:  Displays mean and SD of groups and compares with a test
#               (independent 2-group Mann-Whitney U test)
#---------------------------------------------------------------------------------
calculateGroupDifferences <- function(group1, group2, metric_type)
{
    cat("\n\n", rep("-", 40), "\n")
    cat("Metric type: ", metric_type)
    cat("\nGroup 1: mean = ", mean(unlist(group1)), " SD = ", sd(unlist(group1)))
    cat("\nGroup 2: mean = ", mean(unlist(group2)), " SD = ", sd(unlist(group2)), "\n")
    diff <- wilcox.test(unlist(group1), unlist(group2))
    print(diff)
    cat("\n", rep("-", 40), "\n")
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
# FUNCTION:     displayBoxPlotAndDataParticipant()
# INPUT:        data.frame, vector, int
# OUTPUT:       data.frame
# DESCRIPTION:  Produces a boxplot for all participants for a condition.
#               Additionally displays group mean and SD
#---------------------------------------------------------------------------------
displayBoxPlotAndDataParticipant <- function(data, participants, condition)
{
    # convert list into data frame object
    fixation_data <- data.frame(matrix(unlist(data), ncol = length(data), byrow = FALSE), stringsAsFactors = FALSE)
    
    # rename the columns, print mean/SD and make boxplot
    colnames(fixation_data) <- participants
    print(boxplot(fixation_data, main = paste0("Condition ", condition), xlab = "Participant", 
                  ylab = "Fixation duration (ms)", par(las = 2, cex = 0.7)))

    return(fixation_data)
}

#---------------------------------------------------------------------------------
# FUNCTION:     processQuestionResults()
# INPUT:        list
# OUTPUT:       
# DESCRIPTION:  Output a histogram of the responses to the question about narration
#              
#---------------------------------------------------------------------------------
processQuestionResults <- function(responses)
{
    # output histogram of responses
    question_data <- as.vector(responses[["question"]])
    print(question_data)
    print(ggplot(data.frame(question_data), aes(x = question_data)) + geom_histogram(fill = "black"))
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
    group_fixations <- list()
    data_transform <- TRUE
    k <- 1
    
    # loop over conditions 
    for(i in 1:length(data))
    {
        cat("\n\nCondition = ", i, "\n\n")
        participant_data <- list()
        
        # get condition and subset
        cond_data <- data[[paste0(cond_str, i)]]
      
        # if condition 2 extract the responses to question before subsetting
        if(i == 2) question_responses <- extractQuestionResponse(cond_data)
            
        # subset the fixation duration data
        cond_data <- extractFixationDurationData(cond_data)
        
        # get unique participants and name list elements the same
        participants <- getStudyParticipants(cond_data)
        participant_data <- setNames(vector("list", length(participants)), participants)
        
        # generate box plots
        boxPlotsPerStimulus(cond_data, i)
        displayViolinPlots(cond_data, i)
        
        # loop over stimuli (paintings)
        for(j in 1:length(stimuli))
        {
            # get the current stimuli
            current_stimuli <- extractStimulus(cond_data, stimuli[j])
            current_stimuli <- current_stimuli[!duplicated(current_stimuli), ]
            
            stimulus_data[[k]] <- current_stimuli
            k <- k + 1
            
            # if we want to transform the data
            if(data_transform)
            {
                transformData(current_stimuli, "log")
            }
            
            
            # run normaility tests on data
            normalityTests(current_stimuli, i, stimuli[j])
            
            for(l in 1:length(participants))
            {
                # aggregate fixation data accross all stimuli
                participant_FD <- current_stimuli[current_stimuli$ParticipantName == participants[l], ]
                if(length(participant_FD$GazeEventDuration) > 0)
                {
                    participant_data[[participants[l]]] <- c(participant_data[[participants[l]]], 
                                                             as.vector(participant_FD$GazeEventDuration))
                }
            }
        }
        
        # shows FD accross all stimuli
        group_fixations[[paste0("condition_", i)]] <- displayBoxPlotAndDataParticipant(participant_data, participants, i)
    }
    calculateGroupDifferences(group_fixations[["condition_1"]], group_fixations[["condition_2"]], "Fixation duration")
    processQuestionResults(question_responses)
}

# run main
main()
