import numpy as np
from scipy.stats import norm, rankdata
from sklearn.mixture import GaussianMixture
from tqdm import tqdm

def empirical_cdf(x, y, X, Y):
    return np.mean((X <= x) & (Y <= y))

def G_uniform(x, y):
        return x * y

def ks_2d_statistic(X, Y, G):

    n = len(X)
    X = np.clip(X,0,1)
    Y = np.clip(Y,0,1)
    
    # maximum of observed points
    D1 = max(
        empirical_cdf(X[i], Y[i], X, Y) - G(X[i], Y[i])
        for i in range(n)
    )

    # maximum distance over all intersection points
    D2 = max(
            empirical_cdf(X[j], Y[i], X, Y) - G(X[j], Y[i])
            for j in range(n) for i in range(n) 
            if (X[j]>X[i] and Y[j]<Y[i])
        )

    # minimum distance over all intersection poitns (with 2/n correction)
    D3 = (2 / n) - min(
                        empirical_cdf(X[j], Y[i], X, Y) - G(X[j], Y[i])
                        for i in range(n) for j in range(n)
                        if (X[j] > X[i] and Y[j] < Y[i])
    )

    # maximum distance among projections of observed points on the right boundary (x = 1)
    D4 = (1/n) - min(
                    empirical_cdf(1, Y[i], X, Y) - G(1, Y[i])
                    for i in range(n)
    )

    # maximum distance among projections of the observed points on top boundary (y = 1)
    D5 = (1/n) - min(
                    empirical_cdf(X[i], 1, X, Y) - G(X[i], 1)
                    for i in range(n)
    )

    # final statistic
    Dn = max(D1, D2, D3, D4, D5)

    return Dn

num_sim = 5000

mean0 = np.array([0,0])
mean1 = np.array([3,3])   
cov = np.array([[1, 0.5], [0.5, 1]])

c_alpha = 0.4141
power_values = []
n = 15
epsilon_list = [0.1, 0.2, 0.4]
for epsilon in epsilon_list:
    Dn_list = []
    for k in tqdm(range(num_sim)):
        # alternative distribution
        gmm = GaussianMixture(n_components=2)
        gmm.weights_= np.array([1-epsilon, epsilon])
        gmm.means_= np.array([mean0,mean1])
        gmm.covariances_= np.array([cov,cov])
        gmm.precisions_cholesky_ = np.linalg.cholesky(np.linalg.inv(cov))[None, :, :].repeat(2, axis=0)
        samples, _ = gmm.sample(n)

        X1 = samples[:,0]
        X2 = samples[:,1]

        Y1 = norm(loc=0, scale=1).cdf(X1)
        Y2 = norm(loc=0.5*X1, scale=np.sqrt(0.75)).cdf(X2)
        
        Dn_list.append(ks_2d_statistic(Y1, Y2, G_uniform))

    power = np.mean(np.array(Dn_list)>c_alpha)
    power_values.append(power)

print("power values:", power_values)