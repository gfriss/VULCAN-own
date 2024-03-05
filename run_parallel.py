import os
import subprocess # to run bash commands like in a terminal
# this is how to use it:
# subprocess.run(['everything', 'you', 'would', 'type', 'word', 'by', 'word'], check = True, text = True)
from mpi4py.MPI import COMM_WORLD as CW # for paralellisation

rank = CW.Get_rank()
size = CW.Get_size()

nsim = 12
sim_per_rank = int(nsim / size) # this is needed to distribure the tasks between tha CPUs

main_folder = '/lustre/astro/gfriss/test_many' # place to store outputs

for i in range(rank*sim_per_rank, (rank+1)*sim_per_rank):   # paralellisation itself, it spreads the task between the CPUs
                                                            # this is the magic, after this just think of it as a normal, sequential loop
    if i < 10:
        new_folder = os.path.join(main_folder, 'sim_0' + str(i)) # my folder naming conventions
    else:
        new_folder = os.path.join(main_folder, 'sim_' + str(i))
    # run the simulation in the given folder
    subprocess.check_call(['python vulcan.py'])


# TO DO:
#   need something to prepare atmosphere and find a good place to store that, chemistry will be the same
#   similar for boundary conditions
#   make a generalised, probably with inputs, script for initial mixing ratios
#   write general vulcan_cfg generator
#   do I need to copy vulcan.py for each function or can I call the same file?
#   after all these can the parallelisation run
#   plus figure out the output structure: how to name and what groupings should I save them in