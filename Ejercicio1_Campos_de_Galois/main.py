from galois_field import GaloisField


def main():

    field = GaloisField(4, 0b1100)

    print("Elementos de GF(2^4):")
    print()

    for exponent, value in enumerate(field.exponent_table):
        element = field[value]
        print(f"alpha^{exponent:2d} = {element}")

    print("0         =", field[0])

    print("\nOperaciones de prueba:")
    print()

    a = field[0b0110]  # alpha^5
    b = field[0b1101]  # alpha^7

    print("a =", a, "= alpha^5")
    print("b =", b, "= alpha^7")
    print()

    print("a + b =", a + b)
    print("a * b =", a * b)
    print("a / b =", a / b)
    print("a^3   =", a ** 3)
    print("a^-1  =", a.inverse())

    print("\nVerificaciones:")
    print()

    print("a * a^-1 =", a * a.inverse())
    print("a / a     =", a / a)
    print("a + a     =", a + a)
    print("0 * a     =", field[0] * a)


if __name__ == "__main__":
    main()