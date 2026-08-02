from pathlib import Path
import sys

sys.path.append(str(Path(__file__).resolve().parent.parent))

from Ejercicio1_Campos_de_Galois.galois_field import GaloisField


class GFPoly:

    def __init__(self, field, coefficients):
        if not isinstance(field, GaloisField):
            raise TypeError(
                "field debe ser una instancia de GaloisField."
            )

        if not isinstance(coefficients, (list, tuple)):
            raise TypeError(
                "Los coeficientes deben entregarse en una lista o tupla."
            )

        if len(coefficients) == 0:
            raise ValueError(
                "La lista de coeficientes no puede estar vacía."
            )

        self.field = field

        for coefficient in coefficients:
            self._check_coefficient(coefficient)

        self.coefficients = list(coefficients)

        self._remove_leading_zeros()

    def _check_coefficient(self, coefficient):

        if not isinstance(coefficient, GaloisField):
            raise TypeError(
                "Todos los coeficientes deben ser elementos "
                "de GaloisField."
            )

        if (
            coefficient.m != self.field.m
            or coefficient.primitive_polynomial
            != self.field.primitive_polynomial
        ):
            raise ValueError(
                "Todos los coeficientes deben pertenecer "
                "al mismo campo."
            )

    def _check_polynomial(self, other):

        if not isinstance(other, GFPoly):
            raise TypeError(
                "La operación requiere otro polinomio GFPoly."
            )

        if (
            self.field.m != other.field.m
            or self.field.primitive_polynomial
            != other.field.primitive_polynomial
        ):
            raise ValueError(
                "Los polinomios pertenecen a campos diferentes."
            )

    def _remove_leading_zeros(self):

        while (
            len(self.coefficients) > 1
            and self.coefficients[0].value == 0
        ):
            self.coefficients.pop(0)

    @property
    def degree(self):

        if self.is_zero():
            return -1

        return len(self.coefficients) - 1

    def is_zero(self):

        return (
            len(self.coefficients) == 1
            and self.coefficients[0].value == 0
        )

    def __add__(self, other):

        self._check_polynomial(other)

        length = max(
            len(self.coefficients),
            len(other.coefficients)
        )

        zero = self.field[0]

        coefficients_a = (
            [zero] * (length - len(self.coefficients))
            + self.coefficients
        )

        coefficients_b = (
            [zero] * (length - len(other.coefficients))
            + other.coefficients
        )

        result = [
            a + b
            for a, b in zip(coefficients_a, coefficients_b)
        ]

        return GFPoly(self.field, result)

    def __sub__(self, other):

        return self + other

    def __mul__(self, other):

        self._check_polynomial(other)

        if self.is_zero() or other.is_zero():
            return GFPoly(self.field, [self.field[0]])

        result_length = (
            len(self.coefficients)
            + len(other.coefficients)
            - 1
        )

        result = [
            self.field[0]
            for _ in range(result_length)
        ]

        for i, coefficient_a in enumerate(self.coefficients):
            for j, coefficient_b in enumerate(other.coefficients):
                result[i + j] = (
                    result[i + j]
                    + coefficient_a * coefficient_b
                )

        return GFPoly(self.field, result)

    def divmod(self, divisor):

        self._check_polynomial(divisor)

        if divisor.is_zero():
            raise ZeroDivisionError(
                "No se puede dividir por el polinomio nulo."
            )

        if self.degree < divisor.degree:
            quotient = GFPoly(
                self.field,
                [self.field[0]]
            )

            remainder = GFPoly(
                self.field,
                self.coefficients
            )

            return quotient, remainder

        quotient_length = self.degree - divisor.degree + 1

        quotient_coefficients = [
            self.field[0]
            for _ in range(quotient_length)
        ]

        remainder = GFPoly(
            self.field,
            self.coefficients
        )

        while (
            not remainder.is_zero()
            and remainder.degree >= divisor.degree
        ):
            leading_coefficient = (
                remainder.coefficients[0]
                / divisor.coefficients[0]
            )

            degree_difference = (
                remainder.degree - divisor.degree
            )

            quotient_position = (
                quotient_length - degree_difference - 1
            )

            quotient_coefficients[quotient_position] = (
                leading_coefficient
            )

            term_coefficients = (
                [leading_coefficient]
                + [self.field[0]] * degree_difference
            )

            term = GFPoly(
                self.field,
                term_coefficients
            )

            remainder = remainder + term * divisor

        quotient = GFPoly(
            self.field,
            quotient_coefficients
        )

        return quotient, remainder

    def __floordiv__(self, other):

        quotient, _ = self.divmod(other)

        return quotient

    def __mod__(self, other):

        _, remainder = self.divmod(other)

        return remainder

    def scale(self, scalar):

        self._check_coefficient(scalar)

        result = [
            coefficient * scalar
            for coefficient in self.coefficients
        ]

        return GFPoly(self.field, result)

    def evaluate(self, x):

        self._check_coefficient(x)

        result = self.field[0]

        for coefficient in self.coefficients:
            result = result * x + coefficient

        return result

    def __call__(self, x):

        return self.evaluate(x)

    @classmethod
    def from_roots(cls, field, roots):

        if not isinstance(roots, (list, tuple)):
            raise TypeError(
                "Las raíces deben entregarse en una lista o tupla."
            )

        polynomial = cls(
            field,
            [field[field.one_value]]
        )

        one = field[field.one_value]

        for root in roots:
            if not isinstance(root, GaloisField):
                raise TypeError(
                    "Todas las raíces deben ser elementos "
                    "de GaloisField."
                )

            if (
                root.m != field.m
                or root.primitive_polynomial
                != field.primitive_polynomial
            ):
                raise ValueError(
                    "Todas las raíces deben pertenecer "
                    "al mismo campo."
                )

            factor = cls(
                field,
                [one, root]
            )

            polynomial = polynomial * factor

        return polynomial

    def __eq__(self, other):

        if not isinstance(other, GFPoly):
            return False

        same_field = (
            self.field.m == other.field.m
            and self.field.primitive_polynomial
            == other.field.primitive_polynomial
        )

        return (
            same_field
            and self.coefficients == other.coefficients
        )

    def __repr__(self):

        coefficients = [
            format(coefficient.value, f"0{self.field.m}b")
            for coefficient in self.coefficients
        ]

        return f"GFPoly({coefficients})"