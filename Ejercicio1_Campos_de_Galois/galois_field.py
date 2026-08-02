class GaloisField:

    def __init__(
        self,
        m,
        primitive_polynomial,
        value=0,
        exponent_table=None,
        logarithm_table=None
    ):
        if not isinstance(m, int) or m <= 0:
            raise ValueError("m debe ser un entero positivo.")

        if not isinstance(primitive_polynomial, int):
            raise TypeError(
                "El polinomio primitivo debe representarse como un entero."
            )

        self.m = m
        self.primitive_polynomial = primitive_polynomial
        self.order = 2 ** m
        self.value = value

        if not 0 <= primitive_polynomial < self.order:
            raise ValueError(
                f"El polinomio primitivo debe estar entre "
                f"0 y {self.order - 1}."
            )

        if not isinstance(value, int):
            raise TypeError(
                "El valor del elemento debe ser un entero."
            )

        if not 0 <= value < self.order:
            raise ValueError(
                f"El elemento debe estar entre 0 y {self.order - 1}."
            )

        self.one_value = 1 << (self.m - 1)

        if exponent_table is None or logarithm_table is None:
            self.exponent_table, self.logarithm_table = (
                self._build_tables()
            )

            self.elements = [
                GaloisField(
                    m=m,
                    primitive_polynomial=primitive_polynomial,
                    value=i,
                    exponent_table=self.exponent_table,
                    logarithm_table=self.logarithm_table
                )
                for i in range(self.order)
            ]

            for element in self.elements:
                element.elements = self.elements

        else:
            self.exponent_table = exponent_table
            self.logarithm_table = logarithm_table

    def _build_tables(self):

        exponent_table = [0] * (self.order - 1)
        logarithm_table = [None] * self.order

        current = self.one_value

        for exponent in range(self.order - 1):

            if logarithm_table[current] is not None:
                raise ValueError(
                    "El polinomio ingresado no es primitivo."
                )

            exponent_table[exponent] = current
            logarithm_table[current] = exponent

            current = self._multiply_by_alpha(current)

        if current != self.one_value:
            raise ValueError(
                "El polinomio ingresado no es primitivo."
            )

        return exponent_table, logarithm_table

    def _multiply_by_alpha(self, value):

        overflow = value & 1

        value >>= 1

        if overflow:

            value ^= self.primitive_polynomial

        return value

    def _check_other(self, other):

        if not isinstance(other, GaloisField):
            raise TypeError(
                "La operación requiere otro elemento de Galois."
            )

        if (
            self.m != other.m
            or self.primitive_polynomial
            != other.primitive_polynomial
        ):
            raise ValueError(
                "Los elementos pertenecen a campos diferentes."
            )

    def __getitem__(self, value):

        if not isinstance(value, int):
            raise TypeError(
                "El valor del elemento debe ser un entero."
            )

        if not 0 <= value < self.order:
            raise ValueError(
                f"El elemento debe estar entre 0 y {self.order - 1}."
            )

        return self.elements[value]

    def __add__(self, other):

        self._check_other(other)

        result = self.value ^ other.value

        return self.elements[result]

    def __sub__(self, other):

        return self + other

    def __mul__(self, other):

        self._check_other(other)

        if self.value == 0 or other.value == 0:
            return self.elements[0]

        exponent_a = self.logarithm_table[self.value]
        exponent_b = self.logarithm_table[other.value]

        result_exponent = (
            exponent_a + exponent_b
        ) % (self.order - 1)

        result_value = self.exponent_table[result_exponent]

        return self.elements[result_value]

    def inverse(self):

        if self.value == 0:
            raise ZeroDivisionError(
                "El elemento 0 no tiene inverso multiplicativo."
            )

        exponent = self.logarithm_table[self.value]

        inverse_exponent = (
            -exponent
        ) % (self.order - 1)

        result_value = self.exponent_table[inverse_exponent]

        return self.elements[result_value]

    def __floordiv__(self, other):

        self._check_other(other)

        if other.value == 0:
            raise ZeroDivisionError(
                "No se puede dividir por cero."
            )

        if self.value == 0:
            return self.elements[0]

        exponent_a = self.logarithm_table[self.value]
        exponent_b = self.logarithm_table[other.value]

        result_exponent = (
            exponent_a - exponent_b
        ) % (self.order - 1)

        result_value = self.exponent_table[result_exponent]

        return self.elements[result_value]

    def __truediv__(self, other):

        return self // other

    def __pow__(self, exponent):

        if not isinstance(exponent, int):
            raise TypeError(
                "El exponente debe ser un entero."
            )

        if exponent < 0:
            raise ValueError(
                "El exponente debe ser mayor o igual que cero."
            )

        # Por convención, a^0 = 1.
        if exponent == 0:
            return self.elements[self.one_value]

        # Para n > 0, 0^n = 0.
        if self.value == 0:
            return self.elements[0]

        base_exponent = self.logarithm_table[self.value]

        result_exponent = (
            base_exponent * exponent
        ) % (self.order - 1)

        result_value = self.exponent_table[result_exponent]

        return self.elements[result_value]

    def __eq__(self, other):

        if not isinstance(other, GaloisField):
            return False

        return (
            self.m == other.m
            and self.primitive_polynomial
            == other.primitive_polynomial
            and self.value == other.value
        )

    def __int__(self):

        return self.value

    def __repr__(self):

        binary = format(self.value, f"0{self.m}b")

        return f"GF({self.value}, 0b{binary})"