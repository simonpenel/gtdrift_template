import sys
import re

reader = open(sys.argv[1], 'r')
writer = open(sys.argv[2], 'w')
data = reader.readlines()
print("Removing low comlexity from "+sys.argv[1])
for line in data:
	#print("===="+line[0])
	if line[0] != '>':
	    line = re.sub('[a-z]','N',line)
	writer.write(line)
