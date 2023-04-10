# Build/Rebuild Container (Ensure Overwritten Container Not In Use!)
cd /fs1/project/cgm_world/apptainer
apptainer build tensorflow_2.9.2_custom.sif tensorflow_2.9.2_custom.def

# Build From Upstream With No Customization
apptainer build tensorflow_2.9.2_gpu.sif docker://tensorflow/tensorflow:2.10.0-gpu

# Run Interactively (with gpu '--nv')
apptainer shell --nv ./tensorflow_2.9.2_gpu.sif

# Run In Job (with gpu)
apptainer exec --nv ./tensorflow_2.9.2_gpu.sif python --version

# Add additional packages
vim tensorflow_2.9.2_gpu.def
apptainer build tensorflow_2.9.2_custom.sif tensorflow_2.9.2_custom.def

