from scipy.optimize import fsolve

def alpha(alpha_ab, alpha_star, r, x_a):
    """
    Iterative solve of implicit function for seperation factor, alpha for membrane seperation
    
    Parameters:
    - alpha_current: current estimate of alpha
    - alpha_star: initial or guessed value of alpha* (constant in the equation)
    - x_A: mole fraction of component A
    - r: relative flow rate or any other constant in the denominator

    Returns:
    - Next value of alpha

    """
    numerator = x_a * (alpha_ab - 1) + 1 - r * alpha_ab
    denominator = x_a * (alpha_ab - 1) + 1 - r
    return alpha_ab - alpha_star * (numerator / denominator)


def main():
    # example usage
    alpha_star = 3.47
    x_a = 0.5
    r = 0.5 / 1
    alpha_ab_guess = 0.5

    print(fsolve(alpha, alpha_ab_guess, args=(alpha_star, x_a, r)))

if "__name__" == "__main__":
    main()

    