linesP = [line.rstrip('\n') for line in open('seq.eps')]
linesQ = [line.rstrip('\n') for line in open('par.eps')]
#linesR = [line.rstrip('\n') for line in open('test3')]

sortP = linesP.sort()
sortQ = linesQ.sort()
#sortR = linesR.sort()

print(sortP == sortQ)
#print(sortQ == sortR)
#print(sortP == sortR)

