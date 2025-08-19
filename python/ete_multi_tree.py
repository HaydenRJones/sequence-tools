import re
import sys
from ete4 import Tree

# What to root on
root_taxa = 'L251'

# What file to open (default is the first cmd argument)
infile = sys.argv[1]

# Get the name of the file minus the extension
file_name = ''.join(infile.split('.')[:-1])

# Empty strings to store rooted trees and the loci name
new_trees = ''
loci_name = ''

# Open the file and loop over each line in it, reading it as a tree.
file = open(infile, 'r')
for line in file.readlines():
    
    # If we don't already know the loci name strip this out. Previously this was the file name however that is not always best practice
    if not loci_name:
        leaf_name = re.search(r'{}_[A-z0-9\-\.\_]*'.format(root_taxa), line)
        loci_name = leaf_name.group()[len(root_taxa):]
        
    t = Tree(line)
    
    t.set_outgroup(f"{root_taxa}{loci_name}")
    
    new_trees += f'{t.write()}\n'

# Write the new trees to a file called .root
with open(f'{file_name}.root', 'w') as f:
    f.write(new_trees)