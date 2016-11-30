##################################################################################
# NAME:         pilot metadata.R
# AUTHOUR:      Alan Davies
# DATE:         08/09/2016
# INSTITUTION:  Interaction Analysis and Modelling Lab (IAM), University of Manchester
# DESCRIPTION:  Stores all the stimuli metadata in a handy list structure negating the
#               need for several functions to detrmine properties.
##################################################################################

# globals
study_stimuli <<- c("rhyl_sands", 
                    "flask_walk", 
                    "self_portrait", 
                    "evening_glows", 
                    "14_6_1964", 
                    "woman_man",
                    "sir_gregory", 
                    "release") 

#---------------------------------------------------------------------------------
# FUNCTION:     getStimuliMetaData()
# INPUT:        string
# OUTPUT:       list
# DESCRIPTION:  Returns a list element for the stimuli required with metadata
#               concerning file name and name conventions.
#---------------------------------------------------------------------------------
getStimuliMetaData <- function(stimulus)
{
    condition <- list()
    
    condition[["rhyl_sands"]] <- list(ref_name = "rhyl_sands", file_name = "1917.170med_resized.jpg.png", label = "Rhyl Sands")
    condition[["flask_walk"]] <- list(ref_name = "flask_walk", file_name = "1922.5med_resized.jpg", label = "Flask walk, hampstead")
    condition[["self_portrait"]] <- list(ref_name = "self_portrait", file_name = "1934.2med.jpg.png", label = "Self-portrait")
    condition[["evening_glows"]] <- list(ref_name = "evening_glows", file_name = "1934.394med_resized.jpg.png", label = "When the west with evening glows")
    condition[["14_6_1964"]] <- list(ref_name = "14_6_1964", file_name = "1964.287med.jpg.png", label = "14.6.1964")
    condition[["woman_man"]] <- list(ref_name = "woman_man", file_name = "1968.173med.jpg.png", label = "Woman and Suspended Man")
    condition[["sir_gregory"]] <- list(ref_name = "sir_gregory", file_name = "1976.79med_resized.jpg.png", label = "Sir Gregory Page-Turner")
    condition[["release"]] <- list(ref_name = "release", file_name = "1995.184med.jpg.png", label = "Release")
                            
    return(condition[[stimulus]])   
}