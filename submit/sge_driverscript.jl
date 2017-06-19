#sge_driverscript.jl

#Word of warning: Qsub will send mail for EVERY task in an array job
#Hence why I've not done so for the norm, model and plots jobs
email = "A.M.Phillips@liverpool.ac.uk"

normChains = 10
normJobs = 1000

modelChains = 100
modelJobs = 1000

plotsJobs = 1000

#########################################
# Import Setup
#########################################
open("bayesprot-import-setup.sh","w") do f
  write(f,"#!/bin/bash\n")
  write(f,"#\$ -M "*email*"\n")
  write(f,"#\$ -m bes\n")
  write(f,"#\$ -V -cwd\n")
  write(f,"julia bayesprot-import-setup.jl")
end

qsubReturn = readall(`qsub bayesprot-import-setup.sh`)
importSetupJobID = parse(Int,match(r"(?<= job )\d+",qsubReturn).match)

#########################################
# Import
#########################################
open("bayesprot-import.sh","w") do f
  write(f,"#!/bin/bash\n")
  write(f,"#\$ -V -cwd\n")
  write(f,"#\$ -o import/out\n")
  write(f,"#\$ -e import/error\n")
  write(f,"#\$ -M "*email*"\n")
  write(f,"#\$ -m bes\n")
  write(f,"#\$ -l h_vmem=8G,h_rt=01:00:00\n")
  write(f,"cd import/results\n")
  write(f,"../../../bin/Rscript ../../import.R HPC")
end

qsubReturn = readall(`qsub -hold_jid $importSetupJobID bayesprot-import.sh`)
importJobID = parse(Int,match(r"(?<= job )\d+",qsubReturn).match)

#########################################
# Norm Setup
#########################################
open("bayesprot-norm-setup.sh","w") do f
  write(f,"#!/bin/bash\n")
  write(f,"#\$ -M "*email*"\n")
  write(f,"#\$ -m bes\n")
  write(f,"#\$ -V -cwd\n")
  write(f,"julia bayesprot-norm-setup.jl $normChains $normJobs")
end

qsubReturn = readall(`qsub -hold_jid $importJobID bayesprot-norm-setup.sh`)
normSetupJobID = parse(Int,match(r"(?<= job )\d+",qsubReturn).match)

#########################################
# Norm
#########################################
open("bayesprot-norm.sh","w") do f
  write(f,"#!/bin/bash\n")
  write(f,"#\$ -V -cwd\n")
  write(f,"#\$ -o norm/out\n")
  write(f,"#\$ -e norm/error\n")
  write(f,"#\$ -l h_vmem=8G,h_rt=01:00:00\n")
  write(f,"sh norm/norm-job\$SGE_TASK_ID.sh")
end

qsubReturn = readall(`qsub -hold_jid $normSetupJobID -t 1-$normJobs bayesprot-norm.sh`)

normJobID = parse(Int,match(r"(?<= job-array )\d+",qsubReturn).match)

#########################################
# Exposures Setup
#########################################
open("bayesprot-exposures-setup.sh","w") do f
  write(f,"#!/bin/bash\n")
  write(f,"#\$ -M "*email*"\n")
  write(f,"#\$ -m bes\n")
  write(f,"#\$ -V -cwd\n")
  write(f,"julia bayesprot-exposures-setup.jl")
end

qsubReturn = readall(`qsub -hold_jid $normJobID bayesprot-exposures-setup.sh`)
exposuresSetupJobID = parse(Int,match(r"(?<= job )\d+",qsubReturn).match)

#########################################
# Exposures
#########################################
open("bayesprot-exposures.sh","w") do f
  write(f,"#!/bin/bash\n")
  write(f,"#\$ -V -cwd\n")
  write(f,"#\$ -o exposures/out\n")
  write(f,"#\$ -e exposures/error\n")
  write(f,"#\$ -M "*email*"\n")
  write(f,"#\$ -m bes\n")
  write(f,"#\$ -l h_vmem=8G,h_rt=01:00:00\n")
  write(f,"cd exposures/results\n")
  write(f,"../../../bin/Rscript ../../exposures.R HPC")
end

qsubReturn  = readall(`qsub -hold_jid $exposuresSetupJobID bayesprot-exposures.sh`)
exposuresJobID = parse(Int,match(r"(?<= job )\d+",qsubReturn).match)

#########################################
# Model Setup
#########################################
open("bayesprot-model-setup.sh","w") do f
  write(f,"#!/bin/bash\n")
  write(f,"#\$ -M "*email*"\n")
  write(f,"#\$ -m bes\n")
  write(f,"#\$ -V -cwd\n")
  write(f,"julia bayesprot-model-setup.jl $modelChains $modelJobs")
end

qsubReturn = readall(`qsub -hold_jid $exposuresJobID bayesprot-model-setup.sh`)
modelSetupJobID = parse(Int,match(r"(?<= job )\d+",qsubReturn).match)

#########################################
# Model
#########################################
open("bayesprot-model.sh","w") do f
  write(f,"#!/bin/bash\n")
  write(f,"#\$ -V -cwd\n")
  write(f,"#\$ -o model/out\n")
  write(f,"#\$ -e model/error\n")
  write(f,"#\$ -l h_vmem=8G,h_rt=01:00:00\n")
  write(f,"sh model/model-job\$SGE_TASK_ID.sh")
end

qsubReturn = readall(`qsub -hold_jid $modelSetupJobID -t 1-$modelJobs bayesprot-model.sh`)
modelJobID = parse(Int,match(r"(?<= job-array )\d+",qsubReturn).match)

#########################################
# Plots Setup
#########################################
open("bayesprot-plots-setup.sh","w") do f
  write(f,"#!/bin/bash\n")
  write(f,"#\$ -M "*email*"\n")
  write(f,"#\$ -m bes\n")
  write(f,"#\$ -V -cwd\n")
  write(f,"julia bayesprot-plots-setup.jl $plotsJobs")
end

qsubReturn = readall(`qsub -hold_jid $modelJobID bayesprot-plots-setup.sh`)
plotsSetupJobID = parse(Int,match(r"(?<= job )\d+",qsubReturn).match)

#########################################
# Plots
#########################################
open("bayesprot-plots.sh","w") do f
  write(f,"#!/bin/bash\n")
  write(f,"#\$ -V -cwd\n")
  write(f,"#\$ -o plots/out\n")
  write(f,"#\$ -e plots/error\n")
  write(f,"#\$ -l h_vmem=8G,h_rt=01:00:00\n")
  write(f,"sh plots-job\$SGE_TASK_ID.sh")
end

qsubReturn = readall(`qsub -hold_jid $plotsSetupJobID -t 1-$plotsJobs bayesprot-plots.sh`)
plotsJobID = parse(Int,match(r"(?<= job-array )\d+",qsubReturn).match)

#########################################
# Output Setup
#########################################
open("bayesprot-output-setup.sh","w") do f
  write(f,"#!/bin/bash\n")
  write(f,"#\$ -M "*email*"\n")
  write(f,"#\$ -m bes\n")
  write(f,"#\$ -V -cwd\n")
  write(f,"julia bayesprot-output-setup.jl")
end

qsubReturn = readall(`qsub -hold_jid $plotsJobID bayesprot-output-setup.sh`)
outputSetupJobID = parse(Int,match(r"(?<= job )\d+",qsubReturn).match)

#########################################
# Output
#########################################
open("bayesprot-output.sh","w") do f
  write(f,"#!/bin/bash\n")
  write(f,"#\$ -V -cwd\n")
  write(f,"#\$ -o output/out\n")
  write(f,"#\$ -e output/error\n")
  write(f,"#\$ -M "*email*"\n")
  write(f,"#\$ -m bes\n")
  write(f,"#\$ -l h_vmem=8G,h_rt=01:00:00\n")
  write(f,"cd output/results\n")
  write(f,"../../../bin/Rscript ../../output.R HPC")
end

run(`qsub -hold_jid $outputSetupJobID bayesprot-output.sh`)
