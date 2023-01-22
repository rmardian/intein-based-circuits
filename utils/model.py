from scipy.optimize import curve_fit
import numpy as np

### mechanistic model ###

def hill_activation(x, K, n, eps):

    K_, n_, eps_ = 10**K, 10**n, 10**eps
    return ((x**n_ + (eps_*(K_**n_))) / (K_**n_ + x**n_))

def hill_activation_single(x, ag, K, n, eps):

    return 10**ag * hill_activation(x, K, n, eps)

def hill_activation_and(x, ag, K, n, eps, K1, n1, eps1, K2, n2, eps2):

    x1, x2 = x
    return 10**ag * hill_activation(hill_activation(x1, K1, n1, eps1)*hill_activation(x2, K2, n2, eps2), K, n, eps)

def hill_activation_and_fixed(x, ag, K, n, eps):

    K1, n1, eps1, K2, n2, eps2 = -0.06, 0.15, -2.06, 0.91, 0.29, -3.11
    x1, x2 = x
    return 10**ag * hill_activation(hill_activation(x1, K1, n1, eps1)*hill_activation(x2, K2, n2, eps2), K, n, eps)

def inverse_hill(rpu, ag, K, n, eps):
    
    ag_, K_, n_, eps_ = 10**ag, 10**K, 10**n, 10**eps
    #return (((rpu * (K_**n_)) - (ag_ * eps_ * (K_**n_))) / (ag_ - rpu))**(1/n_)
    return (((rpu - (ag_ * eps_)) * (K_**n_))/(ag_ - rpu))**(1/n_)

### version 0.4

def hill_activation_and_4(x, ag, K, n, eps, an, ac, dn, aan, ddn, dc, K1, n1, eps1, K2, n2, eps2):

    x1, x2 = x
    b = (ac * hill_activation(x2, K2, n2, eps2)) - (an * hill_activation(x1, K1, n1, eps1)) + ddn
    en1 =  -1 * b + (np.sqrt(((np.power(b, 2) - (4 * dn * aan))/(2 * dn))))
    ec1 =  (1*dc) * ((ac * hill_activation(x2, K2, n2, eps2)) - (an * hill_activation(x1, K1, n1, eps1)) + (dn * en1))

    en2 =  -1 * b - (np.sqrt(((np.power(b, 2) - (4 * dn * aan))/(2 * dn))))
    ec2 =  (1*dc) * ((ac * hill_activation(x2, K2, n2, eps2)) - (an * hill_activation(x1, K1, n1, eps1)) + (dn * en2))

    #return 10**ag * hill_activation((en2 * ec2), K, n, eps)
    return 10**ag * hill_activation((en1 * ec1), K, n, eps)

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