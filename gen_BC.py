import sys

# this script generates the bottom boundary conditions modifying the original file
# argument order (with values divided by comma inside one argument)
# species, flux, vdep, out_file


species = sys.argv[1].split(',')
flux = sys.argv[2].split(',') 
vdep = sys.argv[3].split(',')
out_file = sys.argv[4]

with open(out_file, 'w') as f:
    #f.write('species flux [cm-2 s-1]   v_dep [cm s-1]\n')
    with open('atm/BC_bot_Pearce_B.txt') as og:
        count = 0
        for line in og:
            if count <= len(species) and line[0] != '#' and line.strip() != '':
                found = False
                for i in range(len(species)):
                    if species[i] == line.split()[0]:
                        count += 1
                        f.write('{}\t\t\t{:.1e}\t\t\t\t\t\t{}\n'.format(species[i],flux[i],vdep[i]))
                        found = True
                        break
                if not found:
                        f.write(line)
            else:
                f.write(line)
