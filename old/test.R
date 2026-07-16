### Collaborators (auto-populated):
###
### joshkenney

source("api/dataRequest.R")
source("api/dataFilter.R")

dataRequest("pdi21")

hsc <- dataFilter(pdi21_clean, phenotype = "hsc", site = "Yale")
source("api/legacy/getTaskLegacy.R")
source

source("api/getTask.R")

prl <- getTask("prl")
prl <- getTaskLegacy("prl")
dataRequest("prl")

source("api/dataRequest.R")


dataRequest("sips_n_items")

source("api/getRedcap.R")
sips_n <- getRedcap("sips_n_items")
createSpss(sips_n)
