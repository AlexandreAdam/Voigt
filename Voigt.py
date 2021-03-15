import cmath
import numpy as np
import sys
from scipy.special import wofz # Faddeeva function
from scipy.integrate import quad

def voigt(a:float, u):
    t = [.314240376, .947788391, 1.59768264, 2.27950708, 3.02063703, 3.8897249]
    c = [1.01172805, -.75197147, 1.2557727e-2, 1.00220082e-2, -2.42068135e-4, 5.00848061e-7]
    s = [1.393237, .231152406, -.155351466, 6.21836624e-3, 9.19082986e-5, -6.27525958e-7]
    y1 = a + 1.5
    y2 = y1**2
    y3 = a + 3
    region2 = (a >= 0.85)|(np.abs(u) <= a*18.1 + 1.65) 
    regionExp = np.abs(u) >= 12

    voigt = np.zeros(u.size)
    voigt[regionExp] += np.exp(-u[regionExp]**2)

    for i in range(6):
        r  = u - t[i]
        r2 = r**2
        d  = 1 / (r2 + y2)
        d1 = y1 * d
        d2 = r * d
        voigt[~region2] += a*(c[i]*(r[~region2]*d2[~region2] - 1.5*d1[~region2])  \
                             + s[i]*y3*d2[~region2])/ (r2[~region2] + 2.25)

        r  = u + t[i]
        r2 = r**2
        d  = 1./(r2 + y2)
        d3 = y1 * d
        d4 = r * d
        voigt[~region2] += a*(c[i]*(r[~region2]*d4[~region2] - 1.5*d3[~region2])  \
                            - s[i]*y3*d4[~region2])/(r2[~region2] + 2.25)

        r  = u - t[i]
        d  = 1./(r**2 + y2)
        d1 = y1 * d
        d2 = r * d
        r  = u + t[i]
        d  = 1./(r**2 + y2)
        d3 = y1 * d
        d4 = r * d
        voigt[region2] += c[i]*(d1[region2] + d3[region2]) - s[i]*(d2[region2] - d4[region2])

    return voigt


def voigt_wofz(a:float, u):
    u = np.atleast_1d(u)
    t = 1 / 4 / a**2
    z = (1 - 1j * u/a) / 2 / np.sqrt(t)
    profile = wofz(1j * z).real 
    return profile


def voigt_true(a:float, u):
    # Ground truth, but slow
    u = np.atleast_1d(u)
    @np.vectorize
    def integral(_u):
        integrand = lambda x: np.exp(-x**2) / ((_u - x)**2 + a**2 )
        out = a / np.pi * quad(integrand, -np.inf, np.inf)[0]
        return out
    out = integral(u)
    return out



def test_voigt():
    import matplotlib.pyplot as plt
    from scipy.integrate import trapz, romb
    import time
    fig = plt.figure()
    u, h = np.linspace(-20, 20, 2049, retstep=True)
    a = [0.002778, 2, 5, 10]
    # a = [2, 5, 10]
    fig.add_axes((.1, .3, .8, .6))
    plt.title("Voigt Profile with Faddeeva function")
    colors = ["k", "r", "b", "g"]
    start = time.time()
    for i, _a in enumerate(a):
        plt.plot(u, voigt_wofz(_a, u), color=colors[i], ls="-", label=f"a={_a}")
        # plt.plot(u, voigt(_a, u), color=colors[i], ls="-.", label=f"Pseudo Voigt a={_a}")
        # plt.plot(u, voigt_true(_a, u), color=colors[i], label=f"True profile a={_a}")
    end = time.time()
    print(f"Took {end - start:.3f} s for wofz_func")
    start = time.time()
    for i, _a in enumerate(a):
        voigt_true(_a, u)
    end = time.time()
    print(f"Took {end - start:.3f} s for true_func")
    plt.legend()
    plt.ylabel("H(a, u)")
    fig.add_axes((.1, .1, .8, .2))
    for i, _a in enumerate(a):
        res1 = voigt_wofz(_a, u)
        print(romb(res1, h))
        res2 = voigt_true(_a, u )
        print(romb(res2, h))
        plt.plot(u, (res2 - res1)/res1 * 100, color=colors[i])
    plt.xlabel("u")
    plt.ylim(-1e-3, 1e-3)
    plt.ylabel(r"Relative error ($\%$)")
    plt.savefig("test_voigt.png", bbox_inches="tight")
    plt.show()

if __name__ == '__main__':
    test_voigt()
