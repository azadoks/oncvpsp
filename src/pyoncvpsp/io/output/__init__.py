"""ONCVPSP output I/O package."""

from . import text
from .text import OncvpspOutputError, OncvpspTextParser

__all__: list[str] = [
    "text",
    "OncvpspTextParser",
    "OncvpspOutputError",
]
