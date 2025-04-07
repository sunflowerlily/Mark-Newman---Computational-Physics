from pylab import *

from numpy.random import random

N = 100
tau = 3.053*60
mu = log(2)/tau

z = random(N)

t_dec = -1/mu*log(1-z)
t_dec = sort(t_dec)
decayed = arange(1,N+1)

surrived = -decayed + N

plot(t_dec,surrived)
show()