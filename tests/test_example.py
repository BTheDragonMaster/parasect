# -*- coding: utf-8 -*-

"""Text example module for PARASECT."""

import pytest


def add(a: int, b: int) -> int:
    """Add two numbers.

    :param a: First number.
    :type a: int
    :param b: Second number.
    :type b: int
    :return: Sum of the two numbers.
    :rtype: int
    """
    return a + b


@pytest.fixture
def numbers() -> tuple[int, int]:
    """Return a tuple of two numbers.

    :return: A tuple of two numbers.
    :rtype: tuple[int, int]
    """
    return 3, 5


def test_add(numbers: tuple[int, int]) -> None:
    """Test the add function.

    :param numbers: A tuple of two numbers.
    :type numbers: tuple[int, int]
    :raises AssertionError: If the test fails.
    """
    assert add(*numbers) == 8
