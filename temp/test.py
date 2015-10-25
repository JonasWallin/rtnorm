
from rtnorm import rtnorm
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats


a = 8.3385004885 
b = 15.0602645553
n = 10000

rtnorm_obj = rtnorm()
X = rtnorm_obj.sample(a, b,size=n)
print(X)
x = np.linspace(a,10,n)
plt.plot(x,stats.truncnorm.pdf(x, a, 10),'r')
plt.plot(x,rtnorm_obj.probabilites(x, a, 10),'g')
plt.hist(X,bins =100, normed=1,alpha=.3)
plt.show()