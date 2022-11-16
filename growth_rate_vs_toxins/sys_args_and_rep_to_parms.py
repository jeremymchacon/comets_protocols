import sys
import math

# when a python script is run on the command line like python myscript.py 2,the arguments live in sys.argv. The first arg is the python file and the rest are the following arguments.

args = sys.argv

rep = int(args[1])
print(f"system argument = {rep}")

# imagine treatments go here
growth_rates = [0.25, 0.5, 1]
space_widths = [0.01, 0.001]

# iterate through these to figure out which to use, and the spatial seed
total_treatments = len(growth_rates) * len(space_widths)
seed = math.floor(rep / total_treatments)
trt_num = 1 + rep % total_treatments
print(f"treatment # = {trt_num}")
counter = 0
for gr in growth_rates:
    growth_rate = gr
    for sw in space_widths:
        counter += 1
        space_width = sw
        if counter == trt_num:
            break
    if counter == trt_num:
        break
print(f"spatial seed = {seed}")
print(f"growth rate = {growth_rate}")
print(f"space width = {space_width}\n")


