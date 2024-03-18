import sys

# this script takes the base vulcan_cfg.py file and generates 
# a new one by changing the required varaible (amount not restricted)
# run syntax syntax: 
# python gen_cfg.py variable_to_change_1,new_value_1 variable_to_change_2,new_value_2 (etc. if more needed)

with open('vulcan_cfg.py', 'r') as f:
    with open('test_cfg.py', 'w') as g:
        for line in f:
            line_unchanged = True
            bits = line.split()
            if len(bits) != 0:
                for input in sys.argv:
                    var_val = input.split(',')
                    if bits[0] == var_val[0]:
                        g.write(line.replace(bits[-1], '\'' + var_val[1] + '\''))
                        line_unchanged = False
                        break
            if line_unchanged:
                g.write(line)
            

