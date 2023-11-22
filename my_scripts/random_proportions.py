import random

def genera(quanti_te_ne_servono):
    numbers = sorted([random.randint(0,1000) for x in range(0, quanti_te_ne_servono-1)])
    numbers.append(1000)
    numbers.insert(0,0)
    intervals = [(numbers[i+1] - numbers[i])/10 for i in range(0,len(numbers)-1)]
    return intervals

#Aureus_S,Aureus_R,Coli,Pseudomonas,Arcobacter
#Arcobacter è quello suscettibile
i=0
while i < 20:
    seq = genera(3)
    seq.insert(2,0)
    seq.insert(2,0)
    if sum(seq) == 100:
        print(seq)
        i += 1
i=0