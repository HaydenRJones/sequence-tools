import sys
import re

tree = sys.argv[1]

with open(tree, 'r') as file:
    data = file.read().replace('\n', '')
    
# # #

# Extract the gene name from the file input. This will only work if the file is called something like loci.file!
# There are probably other, smarter ways of doing this but it will work for now.

# Once we have the gene name, we can check if it has an underscore in it by trying to split it on that delimiter.
# If the length is more than 1 we know it has one. What we can then do is replace that with a '-', then carry on as usual.

geneName = ''.join(tree.split('/')[-1].split('.')[:-1])

temp = geneName.split('_')
if len(temp) >= 1:
    
    replace = '-'.join(temp)
    
    data = re.sub(f'{geneName}', f'{replace}', data)
    
# # #

p = re.compile('_[_A-Za-z0-9\-]*')

hits = p.finditer(data)

last_end = 0
new_tree = []

for hit in hits:

    string = hit.group().split('_')[1:] # ignore the first slice as it should be blank anyway
    
    string.reverse() # so we have [loci, haplotype]
    
    join_string = ''.join(string)
    
    # Take everything from the last hit end (0 initally) up until the start of our current hit -4 (the length of the sample ID)
    # Append to a list and then add our hit (LOCIHAPLOTYPE_)
    # Then add in the sample ID using the same 4 index
    
    new_tree.append(data[last_end:hit.start() - 4])
    new_tree.append(f'{join_string}_')
    new_tree.append(data[hit.start() - 4:hit.start()])
    
    last_end = hit.end()

new_tree.append(data[last_end:len(data)])
#print(data)
print(''.join(new_tree))