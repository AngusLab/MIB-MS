#Adding exported data (human.kinome and mouse.kinome)

human.kinome<- read.csv("human.kinome.csv")
mouse.kinome<- read.csv("mouse.kinome.csv")

usethis::use_data(human.kinome)
usethis::use_data(mouse.kinome)


