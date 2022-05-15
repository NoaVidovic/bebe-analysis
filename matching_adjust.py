f = open('matching', 'r')
data = {}

for l in f.readlines()[1:]:
    l = l.split('\t')
    data[int(l[0])] = '\t'.join(l[1:])

s = 'Front\tM1\tM2\tM3\n'
for k in sorted(data.keys()):
    s += f'{k}\t{data[k]}'

f.close()

with open('matching', 'w') as f:
    f.write(s)
