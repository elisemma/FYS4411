import numpy as np

from time import time


#  def block(filename, verbose=False):
def block(x, verbose=False):
    #  x = np.loadtxt(filename)

    n = len(x)
    d = np.log2(n)

    assert n == 2**d, "n must be a power of 2"
    d = int(d)

    s, gamma = np.zeros(d), np.zeros(d)
    mu = np.mean(x)
    t0 = time()

    # estimate the auto-covariance and variances
    # for each blocking transformation
    for i in np.arange(d):
        n = len(x)
        # estimate autocovariance of x
        gamma[i] = (n) ** (-1) * sum((x[0 : (n - 1)] - mu) * (x[1:n] - mu))
        # estimate variance of x
        s[i] = np.var(x)
        # perform blocking transformation
        x = 0.5 * (x[0::2] + x[1::2])

    #  print(s)

    # generate the test observator M_k from the theorem
    M = (np.cumsum(((gamma / s) ** 2 * 2 ** np.arange(1, d + 1)[::-1])[::-1]))[::-1]

    # we need a list of magic numbers
    q = np.array(
        [
            6.634897,
            9.210340,
            11.344867,
            13.276704,
            15.086272,
            16.811894,
            18.475307,
            20.090235,
            21.665994,
            23.209251,
            24.724970,
            26.216967,
            27.688250,
            29.141238,
            30.577914,
            31.999927,
            33.408664,
            34.805306,
            36.190869,
            37.566235,
            38.932173,
            40.289360,
            41.638398,
            42.979820,
            44.314105,
            45.641683,
            46.962942,
            48.278236,
            49.587884,
            50.892181,
        ]
    )

    k = -1
    # use magic to determine when we should have stopped blocking
    for k in np.arange(0, d):
        if M[k] < q[k]:
            break
    if k >= d - 1:
        print("Warning: Use more data")
    #  variance = s[k] / 2 ** (d - k)
    variance = s[k] / (2 ** (d - k))
    stderr = variance**0.5
    if verbose:
        #  print(f"Runtime: {time() - t0} sec")
        time_str = f"{time() - t0:.8f}s"
        print("Blocking Statistics :")
        #  print(filename)
        print(
            f"{'iter': <5}| {'time': <{len(time_str)}} | {'average': <19}| std. error"
        )
        print(f"{k: <5}| {time_str} | {mu: <19}| {stderr}")
    return stderr


if __name__ == "__main__":
    ans = block(np.loadtxt("output/output.dat"))

    print(ans)
