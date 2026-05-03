# Source me to get the correct configure/build/run environment

# Store tracing and disable (module is *way* too verbose)
{ tracing_=${-//[^x]/}; set +x; } 2>/dev/null

module_load() {
  echo "+ module load $*"
  module load $*
}
module_unload() {
  echo "+ module unload $*"
  module unload $*
}
module_purge() {
  echo "+ module purge"
  module purge
}

# Unload all modules to be certain
[[ ${IFS_RUNTIME_ENV:-unset} == "unset" ]] && module_purge

# Load modules
module_load prgenv/nvidia
module_load nvidia/24.5
if [[ ${IFS_RUNTIME_ENV:-unset} == "unset" ]]; then
  module_load python3/3.11.8-01
  module_load netcdf4/4.9.2
  module_load cmake/3.28.3
  module_load ninja/1.11.1
fi

# Record the RPATH in the executable
export LD_RUN_PATH=$LD_LIBRARY_PATH

export ECBUILD_TOOLCHAIN="./toolchain.cmake"

# Restore tracing to stored setting
{ if [[ -n "$tracing_" ]]; then set -x; else set +x; fi } 2>/dev/null
