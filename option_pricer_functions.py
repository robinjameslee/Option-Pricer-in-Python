import math
from math import sqrt
import pandas as pd
import numpy as np
import datetime
from scipy.stats import norm, qmc

def black_scholes_options(S, K, T, sigma, r, q, OptionType='C'):
    d1 = (np.log(S / K) + (r - q) * T) / (sigma * sqrt(T)) + 0.5 * sigma * sqrt(T)
    d2 = (np.log(S / K) + (r - q) * T) / (sigma * sqrt(T)) - 0.5 * sigma * sqrt(T)
    if OptionType == 'C':
        option_price = S * math.exp(-q * T) * norm.cdf(d1) - K * math.exp(-r * T) * norm.cdf(d2)
    else:
        option_price = K * math.exp(-r * T) * norm.cdf(-d2) - S * math.exp(-q * T) * norm.cdf(-d1)
    return option_price

def black_scholes_vega(S, K, T, sigma, r, q, OptionType='C'):
    d1 = (np.log(S / K) + (r - q) * T) / (sigma * sqrt(T)) + 0.5 * sigma * sqrt(T)
    vega = S * np.exp(-q * T) * norm.pdf(d1) * sqrt(T)
    return vega

def iv_calculator(option_premium, S, K, T, r, q, OptionType='C'):
    sigmahat = sqrt(2 * abs((np.log(S / K) + (r - q) * T) / T))
    tol = 1e-8
    nmax = 100
    sigmadiff = 1
    n = 1
    sigma = sigmahat
    while (sigmadiff >= tol and n < nmax):
        bs_option_price = black_scholes_options(S, K, T, sigma, r, q, OptionType)
        vega = black_scholes_vega(S, K, T, sigma, r, q, OptionType)
        increment = (bs_option_price - option_premium) / vega
        sigma -= increment
        n += 1
        sigmadiff = abs(increment)
    return sigma

def american_options_binomial_tree(S, K, T, sigma, r, OptionType='C', n=1000):
    dt = T / n
    u = math.exp(sigma * sqrt(dt))
    d = 1 / u
    p = (math.exp(r * dt) - d) / (u - d)
    discounting_factor = math.exp(-r * dt)
    stock_price = np.zeros((n + 1, n + 1))
    option_price = np.zeros((n + 1, n + 1))
    
    stock_price[0, 0] = S
    
    # calculate stock prices
    for i in range(1, n + 1):
        stock_price[i, 0] = stock_price[i - 1, 0] * u
        for j in range(1, i + 1):
            stock_price[i, j] = stock_price[i - 1, j - 1] * d
    
    # calculate option prices at maturity
    for j in range(n + 1):
        if OptionType == 'C':
            option_price[n, j] = max(0, stock_price[n, j] - K)
        else:
            option_price[n, j] = max(0, K - stock_price[n, j])
    
    # calculate option prices at each node
    for i in range(n - 1, -1, -1):
        for j in range(i + 1):
            if OptionType == 'C':
                option_price[i, j] = max(stock_price[i, j] - K, discounting_factor * (p * option_price[i + 1, j] + (1 - p) * option_price[i + 1, j + 1]))
            else:
                option_price[i, j] = max(K - stock_price[i, j], discounting_factor * (p * option_price[i + 1, j] + (1 - p) * option_price[i + 1, j + 1]))
    return option_price[0, 0]
    
def geometric_asian_options_closed_form(S, K, T, sigma, r, OptionType='C', n=1000):
    sigsqT = sigma * sigma * T * (n + 1) * (2 * n + 1) / (6 * n * n)
    muT = 0.5 * sigsqT + (r - 0.5 * sigma ** 2) * T * (n + 1) / (2 * n)
    d1 = (np.log(S / K) + (muT + 0.5 * sigsqT)) / sqrt(sigsqT)
    d2 = d1 - sqrt(sigsqT)
    
    if OptionType == 'C':
        option_price = math.exp(-r * T) * (S * math.exp(muT) * norm.cdf(d1) - K * norm.cdf(d2))
    else:
        option_price = math.exp(-r * T) * (K * norm.cdf(-d2) - S * math.exp(muT) * norm.cdf(-d1))
    return option_price

def arithmetic_asian_options(S, K, T, sigma, r, OptionType='C', n=1000, monto_carlo_num_paths=100000, use_control_variate=False):
    # initialize variables
    arithPayoff = np.zeros(monto_carlo_num_paths)
    geoPayoff = np.zeros(monto_carlo_num_paths)
    
    # Set random seed
    seed = 100
    np.random.seed(seed)

    # Pre-calculate the drift and generate the samples
    dt = T / n
    drift = np.exp((r - 0.5 * sigma * sigma) * dt)
    growth_factors = drift * np.exp(sigma * np.sqrt(dt) * np.random.standard_normal((monto_carlo_num_paths, n)))
    
    for i in range(monto_carlo_num_paths):
        Spath = np.zeros(n)
        Spath[0] = S * growth_factors[i, 0];
        
        for j in range(1, n):
            Spath[j] = Spath[j - 1] * growth_factors[i, j]

        arithMean = Spath.mean() #  Arithmetic mean
        geoMean = math.exp((1 / n) * sum(np.log(Spath))) # Geometric mean
        
        if OptionType == 'C':
            arithPayoff[i] = math.exp(-r * T) * max(arithMean - K, 0)
            geoPayoff[i] = math.exp(-r * T) * max(geoMean - K, 0)
        else:
            arithPayoff[i] = math.exp(-r * T) * max(K - arithMean, 0)
            geoPayoff[i] = math.exp(-r * T) * max(K - geoMean, 0)
    
    if use_control_variate:
        # Control variate method
        covXY = np.cov(arithPayoff, geoPayoff)[0][1]
        theta = covXY / np.var(geoPayoff)
        geo_close_form = geometric_asian_options_closed_form(S, K, T, sigma, r, OptionType, n)
        Z = arithPayoff + theta * (geo_close_form - geoPayoff)
        mean = np.mean(Z)
        std = np.std(Z)
    else:
        # Standard Monte Carlo
        mean = np.mean(arithPayoff);
        std = np.std(arithPayoff);
    
    ci = [mean - 1.96 * std / sqrt(monto_carlo_num_paths), mean + 1.96 * std / sqrt(monto_carlo_num_paths)]
    return mean, ci

def geometric_basket_options_closed_form(S_1, S_2, K, T, sigma_1, sigma_2, correlation, r, OptionType='C'):
    b0 = sqrt(S_1 * S_2)
    sigma_basket = sqrt(sigma_1  * sigma_1 + sigma_2 * sigma_2 + 2 * correlation * sigma_1 * sigma_2) / 2
    mu_basket = r - 0.5 * (sigma_1 * sigma_1 + sigma_2 * sigma_2) / 2 + sigma_basket * sigma_basket / 2
    
    d1 = (np.log(b0 / K) + (mu_basket + 0.5 * sigma_basket * sigma_basket) * T) / (sigma_basket * sqrt(T))
    d2 = d1 - sigma_basket * sqrt(T)
    
    if OptionType == 'C':
        option_price = math.exp(-r * T) * (b0 * np.exp(mu_basket * T) * norm.cdf(d1) - K * norm.cdf(d2))
    else:
        option_price = math.exp(-r * T) * (K * norm.cdf(-d2) - b0 * np.exp(mu_basket * T) * norm.cdf(-d1))
    return option_price

def arithmetic_basket_options(S_1, S_2, K, T, sigma_1, sigma_2, correlation, r, OptionType='C', monto_carlo_num_paths=100000, use_control_variate=False):
    # Set random seed
    seed = 100
    np.random.seed(seed)

    Z1 = np.random.standard_normal(monto_carlo_num_paths)
    Z2 = correlation * Z1 + sqrt(1 - correlation * correlation) * np.random.standard_normal(monto_carlo_num_paths)

    S_1_samples = S_1 * np.exp((r - 0.5 * sigma_1 * sigma_1) * T + sigma_1 * sqrt(T) * Z1)
    S_2_samples = S_2 * np.exp((r - 0.5 * sigma_2 * sigma_2) * T + sigma_2 * sqrt(T) * Z2)
    
    arithMean = (S_1_samples + S_2_samples) / 2
    geoMean = np.sqrt(S_1_samples * S_2_samples)
        
    if OptionType == 'C':
        arithPayoff = math.exp(-r * T) * np.maximum(arithMean - K, 0)
        geoPayoff = math.exp(-r * T) * np.maximum(geoMean - K, 0)
    else:
        arithPayoff = math.exp(-r * T) * np.maximum(K - arithMean, 0)
        geoPayoff = math.exp(-r * T) * np.maximum(K - geoMean, 0)
    
    if use_control_variate:
        # Control variate method
        covXY = np.cov(arithPayoff, geoPayoff)[0][1]
        theta = covXY / np.var(geoPayoff)
        geo_close_form = geometric_basket_options_closed_form(S_1, S_2, K, T, sigma_1, sigma_2, correlation, r, OptionType)
        Z = arithPayoff + theta * (geo_close_form - geoPayoff)
        mean = np.mean(Z)
        std = np.std(Z)
    else:
        # Standard Monte Carlo
        mean = np.mean(arithPayoff);
        std = np.std(arithPayoff);
    
    ci = [mean - 1.96 * std / sqrt(monto_carlo_num_paths), mean + 1.96 * std / sqrt(monto_carlo_num_paths)]
    return mean, ci

def kiko_options(S, K, T, sigma, r, barrier_lower, barrier_upper, rebate, n):
    # set the random seed
    seed = 100
    np.random.seed(seed)
    
    # generate the paths of stock prices
    n = int(n)
    M = int(100000)
    dt = T / n
    deltaS = 0.02
    values = {'S': np.zeros(M), 'up': np.zeros(M), 'down': np.zeros(M)}
    value_local = {}
    sequencer = qmc.Sobol(d=n, seed=seed)
    X = np.array(sequencer.random(n=M)) # uniform samples
    Z = norm.ppf(X) # standard normal samples

    # scaled samples    
    samples = (r - 0.5 * sigma * sigma) * dt + sigma * sqrt(dt) * Z
    all_samples = {
        'S': S * np.exp(np.cumsum(samples, axis=1)),
        'up': (S + deltaS) * np.exp(np.cumsum(samples, axis=1)),
        'down': (S - deltaS) * np.exp(np.cumsum(samples, axis=1)),
    }

    for name, samples in all_samples.items():
        for i, path in enumerate(samples):
            if path.max() >= barrier_upper: # (1) knock-out happened
                ko_time = np.argmax(samples[i] >= barrier_upper)
                values[name][i] = rebate * np.exp(-ko_time * r * dt)
            elif path.min() <= barrier_lower: # (2) knock-in happend
                values[name][i] = np.exp(- r * T) * max(K - path[-1], 0)
        value_local[name] = np.mean(values[name])

    mean = np.mean(values['S'])
    std = np.std(values['S'])
    ci = [mean - 1.96 * std / sqrt(M), mean + 1.96 * std / sqrt(M)]
    delta = (value_local['up'] - value_local['down']) / (2 * deltaS)
    return mean, delta