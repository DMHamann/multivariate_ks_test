import numpy as np

def empirical_cdf(x, y, X, Y):
    return np.mean((X <= x) & (Y <= y))

def G_uniform(x, y):
        return x * y

def ks_2d_statistic(X, Y, G):
    """
    The two dimensional Kolmogorov-Smirnov test statistic algorithm.
    """

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