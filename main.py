linesP = [line.rstrip('\n') for line in open('par.o')]
linesS = [line.rstrip('\n') for line in open('seq.o')]

sortP = linesP.sort()
sortS = linesS.sort()

print(sortP == sortS)
