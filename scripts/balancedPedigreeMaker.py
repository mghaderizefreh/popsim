#!/usr/bin/python

import numpy as np
import numpy.random

base = int(input("number of individuals in base population: "))
thisgennum = base
thisgenid = np.arange(1,base+1)
allgen = np.hstack((thisgenid.reshape(-1,1), np.zeros((base,1)), np.zeros((base,1)))).astype(int)

ngen = int(input("number of generations (excluding base): "))
numOff = np.zeros(ngen).astype(int);
print("  population of generation 0: %d\n" % base)
sum_ = base
for i in range(ngen):
    numOff[i] = int(input("number of offspring per mating in generation %d: " %(i+1)))
    thisgenpop = (np.prod(numOff[:(i+1)].astype(float)/2)*base).astype(int) 
    sum_ += thisgenpop
    print("  population of generation %d: %d\n  Total population: %d" % ((i+1), thisgenpop, sum_ ))

for g in range(ngen):
    noffpermate = numOff[g]
    nextgennum = thisgennum * noffpermate / 2
    nextgenid = np.arange(thisgenid.max()+1, thisgenid.max()+ nextgennum + 1)
    nextgensi = np.random.choice(thisgenid, size = thisgennum/2, replace = False)
    nextgenda = np.array(list(set(thisgenid) - set(nextgensi)))
    np.random.shuffle(nextgenda)

    nextgensir = np.zeros_like(nextgensi)
    nextgendam = np.zeros_like(nextgenda)
    nextgensir[:] = nextgensi[:]
    nextgendam[:] = nextgenda[:]
    
    for i in range(noffpermate-1):
        np.random.shuffle(nextgensi)
        nextgensir = np.hstack((nextgensir, nextgensi))
        np.random.shuffle(nextgenda)
        nextgendam = np.hstack((nextgendam, nextgenda))

    nextgen = np.hstack((nextgenid.reshape(-1,1), nextgensir.reshape(-1,1), nextgendam.reshape(-1,1)))

    allgen = np.vstack((allgen, nextgen))
    thisgennum = nextgennum
    thisgenid = np.zeros_like(nextgenid); thisgenid[:] = nextgenid[:]

np.savetxt('dummypedigree.txt', allgen, fmt='%i')
