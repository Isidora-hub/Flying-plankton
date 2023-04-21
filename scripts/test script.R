# this is a github test! 
# I will be loading in a data file and making a figure
# to test that everyone can work in github..
# I want you to produce a plot or an output and export it to the figures folder
# DO THIS TASK WITHOUT ALTERING CODE ABOVE OR BELOW WHERE YOU ARE WORKING

#### Basic github stuff and tips for collaborating...
## 0. See here: https://bookdown.org/rasmith3/r4phe/collaboration-in-github.html

#### When using git by yourself...
## 1. Pull the repo before you start anything
## 2. Save your progress directly and by doing *commits*
## 3. Push your commits to the repo OFTEN-- this keeps the repo up to date on your coding progress

#### When working in teams...
## 4. when collaborating your team should pick features that don’t have overlapping code 
##    (this leads to conflicts and overlapping sections--bad!)
## 5. To make changes, clone the repo, then fork it (you work on the fork)
##    then submit changes through doing the "commit" and "push" to your fork

######## Now you try! See chapter 8.4 -onward on collaborating and cloning a repository with RStudio

############## Let's try

# load some data

# load packages
if (!require("pacman")) install.packages("pacman") # for rapid install if not in library
pacman::p_load(ggplot2)

# check your working directory
getwd()

## note the R script is associated with the R project so the files are all adjacent to this Project in the directory

# let's use the cars dataset, which is loaded with R already
str(cars)


############# C Wall coding for a simple figure
##### when you make your own figure, partition it out like this so we can all see it

# simple historgram with some labels
fig.CW<-hist(cars$speed, col="coral", main="CWall test fig: Car Speed", xlab="speed (mph)"); fig.CW

# export in directory and save as pdf
dev.print(pdf, "figures/CWall.testfig.pdf", height=4, width=5)
dev.off()

# make an output table
anova<-aov(dist~speed, data=cars)
summary(anova) # this is the right output

#export it to wordfile
capture.output(summary(anova),file="output/anova.test.csv")

############# 

# simple historgram with some labels
fig.CT<-hist(cars$speed, col="blue", main="CTran test fig: Car Speed", xlab="speed (mph)"); fig.CT

# export in directory and save as pdf
dev.print(pdf, "figures/CTran.testfig.pdf", height=4, width=5)
dev.off()

# =======
#Josh Change 
fig.JD <- hist(cars$speed, col="salmon", main= "JD test fig: Car Speed", xlab="speed(mph)"); fig.JD

############# Isidora test

# simple historgram with some labels
fig.IE<-hist(cars$speed, col="blue", main="IEchenique test fig: Car Speed", xlab="speed (mph)"); fig.IE

# export in directory and save as pdf
dev.print(pdf, "figures/IEchenique.testfig.pdf", height=4, width=5)
dev.off()