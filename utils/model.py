from scipy.optimize import curve_fit
import numpy as np

### mechanistic model ###
def hill_activation(x, K, n, ymin, ymax):
    
    return ymin + (ymax - ymin) * (x**n / (K**n + x**n))

def hill_activation_and(x, K1, K2, n1, n2, ymin1, ymin2, ymax1, ymax2):
    
    x1, x2 = x
    return hill_activation(x1, K1, n1, ymin1, ymax1) * hill_activation(x2, K2, n2, ymin2, ymax2)

### optimization

#objective function
def compute_sse(t, data, p0, bounds, func):
    
    popt, _ = curve_fit(func, t, data, p0=p0, bounds=bounds)
    sim = func(t, *popt)
    return sum([(act - pred)**2 for act, pred in zip(data, sim)])

#generate random numbers from a uniform distribution for the initial guesses
def random_search(n, ts, data, bounds, func):
    
    init_guess = []
    for k in range(n):
        p0 = [np.random.uniform(low=low, high=high) for low, high in zip(bounds[0], bounds[1])]
        error = compute_sse(ts, data, p0, bounds, func)
        init_guess.append((error, p0))
    return sorted(init_guess)[0][1]