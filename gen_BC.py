import sys

# this script generates the bottom boundary conditions
# argument order (with values divided by comma inside one argument)
# species, flux, vdep, out_file


species = sys.argv[1].split(',')
flux = sys.argv[2].split(',') 
vdep = sys.argv[3].split(',')
out_file = sys.argv[4]

with open(out_file, 'w') as f:
    f.write('species flux [cm-2 s-1]   v_dep [cm s-1]\n')
    for i in range(len(species)):
        f.write('{}\t\t{}\t\t{}\n'.format(species[i],flux[i],vdep[i]))
