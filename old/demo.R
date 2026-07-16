### Collaborators (auto-populated):
###
# remediate data from u01 sites against the NDA data dictionary

source("api/ndaRequest.R")
ndaRequest("deldisk")


ndaRequest("form_1")
