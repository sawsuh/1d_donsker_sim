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
    v0minus: float
    v0plus: float
    v1minus: float
    v1plus: float


I = TypeVar("I")
O = TypeVar("O")


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


def compute_psi(operator: Operator, x0: float) -> Callable[[float], float]:
    return lambda x: 2 * integrate(
        function=lambda y: operator.b(y) / (operator.a(y) * operator.rho(y)),
        left=x0,
        right=x,
    )


def compute_v0(
    left: float,
    right: float,
    operator: Operator,
    pm: PlusMinus,
    psi: None | Callable[[float], float],
) -> Callable[[float], float]:
    if psi is None:
        psi_helper = compute_psi(operator=operator, x0=left)
        psi = wraps(psi_helper)(cache(psi_helper))

    def integrand(y: float) -> float:
        return math.exp(psi(y)) / operator.a(y)

    def integrate_from_to(xy: Tuple[float, float]) -> float:
        return integrate(function=integrand, left=xy[0], right=xy[1])

    denom = integrate_from_to((left, right))
    if pm == PlusMinus.Plus:
        return lambda x: integrate_from_to((left, x)) / denom
    elif pm == PlusMinus.Minus:
        return lambda x: integrate_from_to((x, right)) / denom


def compute_v1(
    left: float,
    right: float,
    operator: Operator,
    pm: PlusMinus,
    psi: None | Callable[[float], float],
    v0: None | Callable[[float], float],
) -> Callable[[float], float]:
    # use 1 memoization for psi and 1 memoization for v0
    # initialise both here but allow for non memoized
    if psi is None:
        psi_helper = compute_psi(operator=operator, x0=left)
        psi = wraps(psi_helper)(cache(psi_helper))

    if v0 is None:
        v0_helper = compute_v0(
            left=left, right=right, operator=operator, pm=pm, psi=psi
        )
        v0 = wraps(v0_helper)(cache(v0_helper))
    return lambda _: 0
