



# 9-9-21
# Get a run where I don't use any distal sites
Rscript gatherInputFiles.R --maxNdistalSites 0

# Also run one where I get big peaks at low cutoff, and keep many of them
Rscript gatherInputFiles.R --maxNdistalSites 20 --coaccessCutoff .035 --peakSize 1000



