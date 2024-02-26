#!/bin/bash

# Run this file with 'source euler_module_set_up.sh'

env2lmod # Switch to new environment for euler to get access to  newer openmpi version
module load gcc/8.2.0
module load openmpi/3.0.1
module load python/3.6.4
module load cmake/3.11.0
