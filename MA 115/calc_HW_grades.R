#-----------------------------
#
#   Calcualte HW grades
#
#-----------------------------

#set wd 
setwd("~/Documents/Work/github/Courses/MA 115/hw_grades/")

#read in data frames
missing_classes <- read.csv("attendance.csv")[,c(1, 2, 12)]
offline <- read.csv("offline_HW.csv")
online_hw <- read.csv("online.csv")

#clean up 
colnames(missing_classes) <- c("Last", "First", "Missed_Classes")

exams <- offline[, c(1, 2, 6, 11)]
colnames(exams)[1:2] <- c("Last", "First")

offline_hw <- offline[, -c(6, 11)]
colnames(offline_hw)[1:2] <- c("Last", "First")

colnames(online_hw)[1:2] <- c("Last", "First")

#cast all NA's as zeros
for(i in 3:14){
  online_hw[is.na(online_hw[,i]),i] <- 0
  offline_hw[is.na(offline_hw[,i]),i] <- 0
}

exams[is.na(exams[,3]),3] <- 0
exams[is.na(exams[,4]),4] <- 0


#check dimensions
dim(exams)[1]
dim(offline_hw)
dim(online_hw) 

#online HW has 126 students while the exams and online HW has 203 
# need to match online HW names to missing classes sheet

#create First Last column in missing dataset
missing_classes$full_name <- paste(missing_classes$First, missing_classes$Last)
online_hw$full_name <- paste(online_hw$First, online_hw$Last)
offline_hw$full_name <- paste(offline_hw$First, offline_hw$Last)
exams$full_name <- paste(exams$First, exams$Last)


#set up average online hw datastructures
no.students <- length(online_hw$full_name)
hw_ave_online <- numeric(no.students)
acceptions <- c(55,81)

for(i in 1:no.students){
  grades.here <- online_hw[i, 3:14]
  if(missing_classes$Missed_Classes[i] <= 3 || i %in% acceptions){
    #drop 2 lowest
    min_grades <- order(grades.here)[1:2]
    hw_ave_online[i] <- mean(as.numeric(grades.here[-min_grades]))
    
  }else{
    #drop 1 lowest
    min_grade <- which.min(grades.here)
    hw_ave_online[i] <- mean(as.numeric(grades.here[-min_grade]))
  }
}


#need to joing missing classes dataset and offline/exams dataset
ind_matches <- matrix(NA, nrow = no.students, ncol = 2)
for(i in 1:no.students){
  match <- which(missing_classes$full_name[i] == offline_hw$full_name)
  ind_matches[i,] <- c(i, ifelse(length(match)>0, match, NA))
}

#manually fill in NA values
missing_classes$full_name[is.na(ind_matches[,2])]
missing_ind <- which(is.na(ind_matches[,2]))

ind_matches[missing_ind[1],2] <- 4
ind_matches[missing_ind[2],2] <- 18
ind_matches[missing_ind[3],2] <- 27
ind_matches[missing_ind[4],2] <- 28
ind_matches[missing_ind[5],2] <- NA
ind_matches[missing_ind[6],2] <- 38
ind_matches[missing_ind[7],2] <- 61
ind_matches[missing_ind[8],2] <- 84
ind_matches[missing_ind[9],2] <- 89
ind_matches[missing_ind[10],2] <- 90
ind_matches[missing_ind[11],2] <- 97
ind_matches[missing_ind[12],2] <- 123
ind_matches[missing_ind[13],2] <- 127 
ind_matches[missing_ind[14],2] <- 134
ind_matches[missing_ind[15],2] <- 167
ind_matches[missing_ind[16],2] <- 177
ind_matches[missing_ind[17],2] <- 178
ind_matches[missing_ind[18],2] <- 179
ind_matches[missing_ind[19],2] <- 186
ind_matches[missing_ind[20],2] <- 199

#set up average offline hw datastructures
no.students <- length(online_hw$full_name)
acceptions <- c(67)
hw_ave_offline <- numeric(no.students)
for(i in 1:no.students){
  grades.here <- offline_hw[ind_matches[i,2],3:14]
  if(missing_classes$Missed_Classes[i] <= 3 || i %in% acceptions){
    #drop 2 lowest
    min_grades <- order(grades.here)[1:2]
    hw_ave_offline[i] <- mean(as.numeric(grades.here[-min_grades]))
  }else{
    #drop 1 lowest
    min_grade <- order(grades.here)[1]
    hw_ave_offline[i] <- mean(as.numeric(grades.here[-min_grade]))
  }
}


#-----------------------------
#
#   Create  HW grades dataframe
#
#-----------------------------
exams_subset <- exams[ind_matches[,2], 3:4]
rownames(exams_subset) <- 1:nrow(exams_subset)

grades_final <- cbind(missing_classes$Last, missing_classes$First, missing_classes$Missed_Classes,
                      hw_ave_online, hw_ave_offline, )
colnames(grades_final) <- c("Last", "First", "Missed_Classes", 
                            "Online_HW_Average","Offline_HW_Average",
                            "Exam1", "Exam2")




