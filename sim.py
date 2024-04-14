from typing import Callable, Tuple, Dict, Generic, TypeVar
from dataclasses import dataclass
from enum import Enum
import math
from functools import cache, wraps


class PlusMinus(Enum):
    Plus = 1
    Minus = -1


@dataclass
class Operator:
    a: Callable[[float], float]
    rho: Callable[[float], float]
    b: Callable[[float], float]


@dataclass
class Grid:
    get_adjacent: Callable[[float], Tuple[float, float]]


@dataclass
class CellData:
    time: Tuple[float, float]
    prob: Tuple[float, float]


def integrate(
    function: Callable[[float], float], left: float, right: float, increment=0.001
) -> float:
    # Riemann sums
    x = left
    integral: float = 0
    while x < right:
        integral += function(x) * increment
        x += increment
    return integral


def compute_psi(
    operator: Operator, x0: float, with_cache: bool = False
) -> Callable[[float], float]:
    def out(x: float) -> float:
        return 2 * integrate(
            function=lambda y: operator.b(y) / (operator.a(y) * operator.rho(y)),
            left=x0,
            right=x,
        )

    if with_cache:
        out = wraps(out)(cache(out))
    return out


def compute_v0(
    left: float,
    right: float,
    operator: Operator,
    pm: PlusMinus,
    psi: Callable[[float], float],
    with_cache: bool = False,
) -> Callable[[float], float]:
    # if psi is None:
    #     psi = compute_psi(operator=operator, x0=left, with_cache=True)

    def integrand(y: float) -> float:
        return math.exp(psi(y)) / operator.a(y)

    def integrate_from_to(xy: Tuple[float, float]) -> float:
        return integrate(function=integrand, left=xy[0], right=xy[1])

    denom = integrate_from_to((left, right))
    if pm == PlusMinus.Plus:

        def out(x: float) -> float:
            return integrate_from_to((left, x)) / denom

    elif pm == PlusMinus.Minus:

        def out(x: float) -> float:
            return integrate_from_to((x, right)) / denom

    if with_cache:
        out = wraps(out)(cache(out))
    return out


def compute_v1(
    left: float,
    right: float,
    operator: Operator,
    # pm: PlusMinus,
    psi: Callable[[float], float],
    v0: Callable[[float], float],
    v0plus: Callable[[float], float],
    with_cache: bool = False,
) -> Callable[[float], float]:
    # if psi is None:
    #     psi = compute_psi(operator=operator, x0=left, with_cache=True)
    # if v0 is None:
    #     v0 = compute_v0(
    #         left=left, right=right, operator=operator, pm=pm, psi=psi, with_cache=True
    #     )
    def G(x: float, y: float) -> float:
        if x <= y:
            return (
                2
                * (v0plus(x) - v0plus(left))
                * (v0plus(right) - v0plus(y))
                / (v0plus(right) - v0plus(left))
            )
        else:
            return (
                2
                * (v0plus(y) - v0plus(left))
                * (v0plus(right) - v0plus(x))
                / (v0plus(right) - v0plus(left))
            )

    def integrand(x: float, y: float) -> float:
        return G(x, y) * v0(y) * math.exp(psi(y)) / operator.rho(y)

    def out(x: float) -> float:
        # print("diagnostics for v1 integrand")
        # print(integrand(x, x + (right - x) / 2))
        # print(integrand(x, x - (x - left) / 2))
        return integrate(function=lambda y: integrand(x, y), left=left, right=right)

    if with_cache:
        out = wraps(out)(cache(out))
    return out


def compute_celldata(x: float, G: Grid, L: Operator) -> CellData:
    l, r = G.get_adjacent(x)
    psi = compute_psi(operator=L, x0=l, with_cache=True)
    v0plus = compute_v0(
        left=l, right=r, operator=L, psi=psi, with_cache=True, pm=PlusMinus.Plus
    )
    v0minus = compute_v0(
        left=l, right=r, operator=L, psi=psi, with_cache=True, pm=PlusMinus.Minus
    )
    v1plus = compute_v1(
        left=l, right=r, operator=L, psi=psi, v0=v0plus, v0plus=v0plus, with_cache=True
    )
    v1minus = compute_v1(
        left=l, right=r, operator=L, psi=psi, v0=v0minus, v0plus=v0plus, with_cache=True
    )
    # print(v1plus(x))
    # print(v1minus(x))
    return CellData(
        time=(v1minus(x) / v0minus(x), v1plus(x) / v0plus(x)),
        prob=(v0minus(x), v0plus(x)),
    )


class Simulator:
    pass


if __name__ == "__main__":
    L = Operator(a=lambda _: 1, rho=lambda _: 1, b=lambda _: 0)
    G = Grid(get_adjacent=lambda x: (x - 0.3, x + 0.5))
    print(compute_celldata(x=0, G=G, L=L))
