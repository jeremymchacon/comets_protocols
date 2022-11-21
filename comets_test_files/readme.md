This set of files can be used to test if comets is working. 

Assuming your COMETS\_HOME and GUROBI\_COMETS\_HOME variables are set, you should be able to go into a directory with these four files and type 

comets\_scr current\_script


Then, COMETS should run on the command line, and the textbook model should grow to 0.1.  

The most common errors occur if the variables above were not set (or not set correctly).

This folder also includes a python script for running comets via python -- "comets_hello_world.py"

To run it, you need python loaded and gurobi loaded.

On MSI, that will look something like:

module load python/3.8.0
module load gurobi/9.0.2
unset PYTHONPATH
unset PYTHONHOME
unset PYTHONSTARTUP

python comets_hello_world.py

This should run without errors, and print a bunch of increasing biomasses. 


