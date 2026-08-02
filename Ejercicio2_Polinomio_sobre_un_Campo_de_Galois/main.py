from pathlib import Path
import sys

sys.path.append(str(Path(__file__).resolve().parent.parent))

from Ejercicio1_Campos_de_Galois.galois_field import GaloisField
from gf_poly import GFPoly


def main():
    
    field = GaloisField(4, 0b1100)

    zero = field[0]
    one = field[field.one_value]

    alpha_2 = field[field.exponent_table[2]]
    alpha_3 = field[field.exponent_table[3]]
    alpha_5 = field[field.exponent_table[5]]

    # P(x) = alpha^2*x^2 + alpha^5*x + 1
    p = GFPoly(
        field,
        [alpha_2, alpha_5, one]
    )

    # Q(x) = x + alpha^3
    q = GFPoly(
        field,
        [one, alpha_3]
    )

    print("P(x) =", p)
    print("Q(x) =", q)

    print("\nOperaciones:")
    print("P + Q =", p + q)
    print("P * Q =", p * q)

    quotient = p // q
    remainder = p % q

    print("P // Q =", quotient)
    print("P % Q  =", remainder)

    print("\nVerificación de la división:")
    print("P == Q*(P//Q) + P%Q")
    print(p == q * quotient + remainder)

    print("\nEscalado:")
    print("alpha^2 * P =", p.scale(alpha_2))

    print("\nEvaluación:")
    print("P(alpha^3) =", p(alpha_3))

    print("\nConstrucción a partir de raíces:")

    roots = [alpha_2, alpha_3, alpha_5]

    roots_polynomial = GFPoly.from_roots(
        field,
        roots
    )

    print("R(x) =", roots_polynomial)

    print("\nVerificación de las raíces:")

    for root in roots:
        print(
            f"R({root}) =",
            roots_polynomial(root)
        )


if __name__ == "__main__":
    main()