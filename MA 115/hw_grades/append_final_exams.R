#-----------------------------
#
#   Combine grades
#
#-----------------------------

#read in data
setwd("~/Documents/Work/github/Courses/MA 115/hw_grades/")
hw_grades <- read.csv("hw_grades.csv", stringsAsFactors = FALSE)
exam_grades <- read.csv("final_exam_grades.csv", stringsAsFactors = FALSE)

#append final exam grades
hw_grades$final_exam <- exam_grades[,3]

#drop no final exam grades
drop <- c(22, 28, 36, 61, 75, 101, 102, 108, 109)

#final grade dataframe
finals <- hw_grades[-drop,]

#add in Erin Molloy who did No HWs and took No Exams
finals <- rbind(finals, c("Molloy", "Erin","U23758537", 9, rep(0, 4), NA))
finals <- finals[order(finals$Last), ]

#check dimension
dim(finals)

#write out file
write.csv(finals, "final_grade_sheet.csv", row.names = FALSE)

