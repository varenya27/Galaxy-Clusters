from scipy import stats

p=  1 - stats.chi2.cdf(3.84,1)
print(p)