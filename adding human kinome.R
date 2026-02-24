#Adding exported data (human.kinome and mouse.kinome and proteome)

human.kinome<- read.csv("human.kinome.csv")
mouse.kinome<- read.csv("mouse.kinome.csv")
human.proteome<-read.csv("Proteome.names.csv")
mouse.proteome<- read.csv("Mouse_Proteome.names.csv")

usethis::use_data(human.kinome)
usethis::use_data(mouse.kinome)
usethis::use_data(human.proteome)
usethis::use_data(mouse.proteome)


