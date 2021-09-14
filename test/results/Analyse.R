library(VennDiagram)
library (grid)
library(RColorBrewer)

pic = read.table('PIC_1d0i_ssh.tsv', header = T)
pypic =read.table('Side_SideChaine_1d0i_table.tsv',header = T)

setpic = as.numeric(pic[,'Dd.a'])
setpypic = as.numeric(pypic[,'d.don.acc.'])

#coloration
myCol <- brewer.pal(3, "Pastel2")

# Chart
venn.diagram(
  x = list(setpic, setpypic, setpypic),
  category.names = c("PyPIC" , "PICserver " , ""),
  filename = '#14_venn_diagramm.png',
  output=TRUE,
  
  # Output 
  imagetype="png" ,
  height = 480 , 
  width = 480 , 
  resolution = 300,
  compression = "lzw",
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol,
  
  # Numbers
  cex = .6,
  fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  rotation = 1
)