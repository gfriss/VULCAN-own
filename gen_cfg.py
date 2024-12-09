import sys

# this script takes the base vulcan_cfg.py file and generates 
# a new one by changing the required varaible (amount not restricted)
# run syntax syntax: 
# python gen_cfg.py new_file variable_to_change_1,new_value_1,type variable_to_change_2,new_value_2,type (etc. if more needed)
cfg_to_change = 'vulcan_cfg_ox.py'
if sys.argv[-2] == 'rerun': # second to last defines if it is a rerun
    cfg_to_change = '/scratch/s2555875/output/cfg_' + sys.argv[-1] + '.txt' # last points to which simulation is being rerun

with open(cfg_to_change, 'r') as f:
    with open(sys.argv[1], 'w') as g:
        for line in f:
            line_unchanged = True
            bits = line.split()
            if len(bits) != 0:
                for input in sys.argv[2:]:
                    var_val = input.split(',')
                    if bits[0] == var_val[0]:
                        if var_val[-1] == 'str':
                            g.write(line.replace(bits[2], '\'' + var_val[1] + '\'')) # something = something so bits[2]
                        if var_val[-1] == 'val': # this option is for both numbers and booleans
                            g.write(line.replace(bits[2], var_val[1])) # something = something so bits[2]
                        line_unchanged = False
                        break
            if line_unchanged:
                g.write(line)
            

