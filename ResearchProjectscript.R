#code for the Network Biology Research Project
#MSB1014
#October 2019
#================
# call libraries; give citations 
#===============================
library(RCy3)
citation("RCy3")
citation()
#==============
# Import data from DISEASES Database
#====================================
# Experiments
#=============
diseases.experiments <- read.delim("C:/Users/Caroline/Desktop/Network Biology/Project/diseases-experiments.tsv", header=FALSE)
colnames(diseases.experiments) <- c("ProteinID","GeneName","DOID","Disease","Source","Significance","Score")
#=============
# Knowledge
#=============
diseases.knowledge <- read.delim("C:/Users/Caroline/Desktop/Network Biology/Project/diseases-knowledge.tsv", header=FALSE)
colnames(diseases.knowledge) <- c("ProteinID","GeneName","DOID","Disease","Source","Significance","Score")
# ==============================
#
#
# creating the STRING network
# network is a PPI network with genes from the experiments channel DISEASES database T2DM + PD
#=======================================================================================================
# Experiments
#=============
string.cmd <- paste('string protein query taxonID=9606 limit=0 cutoff=0.5 query="',paste(diseases.experiments$ProteinID, collapse=","),'"',sep="")
commandsGET(string.cmd)
renameNetwork("Diffuse STRING experiment")

loadTableData(diseases.experiments,table.key.column = "query term",data.key.column = "ProteinID")

setVisualStyle("default")

# Set up propagation (Experiments)
#=================================
clearSelection()
disease.score.table <- getTableColumns('node','DOID')
# set seeds from PD gene list; 
PD.nodes <- row.names(disease.score.table)[disease.score.table$DOID == "DOID:14330"]
selectNodes(PD.nodes)
setNodeBorderColorBypass(PD.nodes,"#FF3399")
setNodeBorderWidthBypass(PD.nodes,5)

diffusionBasic()
clearSelection()
setNodeColorMapping("diffusion_output_heat", c(0,0.5), c("#FFFFFF", "#006699"), style="default")
#=======================================================================================================
#
#Knowledge
#=================
# SET UP Propagation on seeds are PD gene list;
# network is  Knowledge-based DISEASES database T2DM + PD
#===========================================================
# Knowledge
#===========
string.cmd.knowledge <- paste('string protein query taxonID=9606 limit=0 cutoff=0.5 query="',paste(diseases.knowledge$ProteinID, collapse=","),'"',sep="")
commandsGET(string.cmd.knowledge)
renameNetwork("STRING knowledge")

loadTableData(diseases.knowledge,table.key.column = "query term",data.key.column = "ProteinID") 

setVisualStyle("default")

# Set up propagation (Knowledge)
#=================================

clearSelection()
disease.score.table.knowledge <- getTableColumns('node','DOID')
PD.nodes.knowledge <- row.names(disease.score.table.knowledge)[disease.score.table.knowledge$DOID == "DOID:14330"]
selectNodes(PD.nodes.knowledge)
setNodeBorderColorBypass(PD.nodes,"#FF3399")
setNodeBorderWidthBypass(PD.nodes,5)

diffusionBasic()
clearSelection()
setNodeColorMapping("diffusion_output_heat", c(0,0.5), c("#FFFFFF", "#006699"), style="default")
#=================================================================================================