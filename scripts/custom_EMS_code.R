## Custom EMS figure script- makes the shade of cell dependent on value of cell (from 0-1)
## K Stark 2018

## make new imagine functions for each of the data frames
#Imagine 3 uses a grey colour palette
imagine3 <- function (comm, col = col, fill = TRUE, 
                      xlab = "", ylab = "", yline = 2, xline = 2, sitenames = rownames(comm), 
                      speciesnames = colnames(comm), binary = TRUE) 
{
  reverse <- nrow(comm):1
  comm <- comm[reverse, ]
  par(mar = c(1, 6, 15, 1))
  
  image(1:dim(comm)[2], 1:dim(comm)[1], t(comm), col = c("#00000000","#52525203","#52525205","#52525208","#5252520A","#5252520D","#5252520F","#52525212","#52525214","#52525217","#5252521A","#5252521C","#5252521F","#52525221","#52525224","#52525226","#52525229","#5252522B","#5252522E","#52525230","#52525233","#52525236","#52525238","#5252523B","#5252523D","#52525240","#52525242","#52525245","#52525247","#5252524A","#5252524D","#5252524F","#52525252","#52525254","#52525257","#52525259","#5252525C","#5252525E","#52525261","#52525263","#52525266","#52525269","#5252526B","#5252526E","#52525270","#52525273","#52525275","#52525278","#5252527A","#5252527D","#52525280","#52525282","#52525285","#52525287","#5252528A","#5252528C","#5252528F","#52525291","#52525294","#52525296","#52525299","#5252529C","#5252529E","#525252A1","#525252A3","#525252A6","#525252A8","#525252AD","#525252B0","#525252B3","#525252B5","#525252B8","#525252BA","#525252BD","#525252BF","#525252C2","#525252C4","#525252C7","#525252C9","#525252CC","#525252CF","#525252D1","#525252D4","#525252D6","#525252D9","#525252DB","#525252DE","#525252E0","#525252E3","#525252E6","#525252E8","#525252EB","#525252ED","#525252F0","#525252F2","#525252F5","#525252F7","#525252FA","#525252FC"), xlab = "", ylab = "", axes = FALSE)
  box()
  if (length(sitenames) > 1) {
    axis(2, at = 1:dim(comm)[1], labels = sitenames, las = 1, 
         cex.axis = 1, lwd.ticks = 0)
  }
  if (length(speciesnames) > 1) {
    axis(3, at = 1:dim(comm)[2], labels = speciesnames, las = 2, 
         cex.axis = 1, lwd.ticks = 0)
  }
  mtext(xlab, 3, cex = 1.5, line = xline)
  mtext(ylab, 2, cex = 1.5, line = yline)
}

## Imagine 4 uses blue colour palette
imagine4 <- function (comm, col = col, fill = TRUE, 
                      xlab = "", ylab = "", yline = 2, xline = 2, sitenames = rownames(comm), 
                      speciesnames = colnames(comm), binary = TRUE, add = TRUE) 
{
  reverse <- nrow(comm):1
  comm <- comm[reverse, ]
  par(mar = c(1, 6, 15, 1))
  
  image(1:dim(comm)[2], 1:dim(comm)[1], t(comm), col = c("#00000000","#21908C03","#21908C05","#21908C08","#21908C0A","#21908C0D","#21908C0F","#21908C12","#21908C14","#21908C17","#21908C1A","#21908C1C","#21908C1F","#21908C21","#21908C24","#21908C26","#21908C29","#21908C2B","#21908C2E","#21908C30","#21908C33","#21908C36","#21908C38","#21908C3B","#21908C3D","#21908C40","#21908C42","#21908C45","#21908C47","#21908C4A","#21908C4D","#21908C4F","#21908C52","#21908C54","#21908C57","#21908C59","#21908C5C","#21908C5E","#21908C61","#21908C63","#21908C66","#21908C69","#21908C6B","#21908C6E","#21908C70","#21908C73","#21908C75","#21908C78","#21908C7A","#21908C7D","#21908C80","#21908C82","#21908C85","#21908C87","#21908C8A","#21908C8C","#21908C8F","#21908C91","#21908C94","#21908C96","#21908C99","#21908C9C","#21908C9E","#21908CA1","#21908CA3","#21908CA6","#21908CA8","#21908CAD","#21908CB0","#21908CB3","#21908CB5","#21908CB8","#21908CBA","#21908CBD","#21908CBF","#21908CC2","#21908CC4","#21908CC7","#21908CC9","#21908CCC","#21908CCF","#21908CD1","#21908CD4","#21908CD6","#21908CD9","#21908CDB","#21908CDE","#21908CE0","#21908CE3","#21908CE6","#21908CE8","#21908CEB","#21908CED","#21908CF0","#21908CF2","#21908CF5","#21908CF7","#21908CFA","#21908CFC"), xlab = "", ylab = "", axes = FALSE, add = TRUE)
  box()
  if (length(sitenames) > 1) {
    axis(2, at = 1:dim(comm)[1], labels = sitenames, las = 1, 
         cex.axis = 1, lwd.ticks = 0)
  }
  if (length(speciesnames) > 1) {
    axis(3, at = 1:dim(comm)[2], labels = speciesnames, las = 2, 
         cex.axis = 1, lwd.ticks = 0)
  }
  mtext(xlab, 3, cex = 1.5, line = xline)
  mtext(ylab, 2, cex = 1.5, line = yline)
}

##Imagine 5 uses purple colour palette
imagine5 <- function (comm, col = col, fill = TRUE, 
                      xlab = "", ylab = "", yline = 2, xline = 2, sitenames = rownames(comm), 
                      speciesnames = colnames(comm), binary = TRUE, add = TRUE) 
{
  reverse <- nrow(comm):1
  comm <- comm[reverse, ]
  par(mar = c(1, 6, 15, 1))
  
  image(1:dim(comm)[2], 1:dim(comm)[1], t(comm), col = c("#00000000","#44015403","#44015405","#44015408","#4401540A","#4401540D","#4401540F","#44015412","#44015414","#44015417","#4401541A","#4401541C","#4401541F","#44015421","#44015424","#44015426","#44015429","#4401542B","#4401542E","#44015430","#44015433","#44015436","#44015438","#4401543B","#4401543D","#44015440","#44015442","#44015445","#44015447","#4401544A","#4401544D","#4401544F","#44015452","#44015454","#44015457","#44015459","#4401545C","#4401545E","#44015461","#44015463","#44015466","#44015469","#4401546B","#4401546E","#44015470","#44015473","#44015475","#44015478","#4401547A","#4401547D","#44015480","44015482","#44015487","#4401548A","#4401548C","#4401548F","#44015491","#44015494","#44015496","#44015499","#4401549C","#4401549E","#440154A1","#440154A3","#440154A6","#440154A8","#440154AD","#440154B0","#440154B3","#440154B5","#440154B8","#440154BA","#440154BD","#440154BF","#440154C2","#440154C4","#440154C7","#440154C9","#440154CC","#440154CF","#440154D1","#440154D4","#440154D6","#440154D9","#440154DB","#440154DE","#440154E0","#440154E3","#440154E6","#440154E8","#440154EB","#440154ED","#440154F0","#440154F2","#440154F5","#440154F7","#440154FA","#440154FC"), xlab = "", ylab = "", axes = FALSE, add = TRUE)
  box()
  if (length(sitenames) > 1) {
    axis(2, at = 1:dim(comm)[1], labels = sitenames, las = 1, 
         cex.axis = 1, lwd.ticks = 0)
  }
  if (length(speciesnames) > 1) {
    axis(3, at = 1:dim(comm)[2], labels = speciesnames, las = 2, 
         cex.axis = 1, lwd.ticks = 0)
  }
  mtext(xlab, 3, cex = 1.5, line = xline)
  mtext(ylab, 2, cex = 1.5, line = yline)
}


### read in data. unfortunately have to read in blues/ greys/ purples separately. 

grays_only3 <- read.csv("grays_only3.csv", header = TRUE, stringsAsFactors = FALSE)
row.names(grays_only3) <- grays_only3[,1]
grays_only3 <- grays_only3[,2:56]


blues_only3 <- read.csv("blues_only3.csv", header = TRUE, stringsAsFactors = FALSE)
row.names(blues_only3) <- blues_only3[,1]
blues_only3 <- blues_only3[,2:56]


purples_only3 <- read.csv("purples_only3.csv", header = TRUE, stringsAsFactors = FALSE)
row.names(purples_only3) <- purples_only3[,1]
purples_only3 <- purples_only3[,2:56]

## Make figure! 

imagine3(grays_only3)
imagine4(blues_only3)
imagine5(purples_only3)