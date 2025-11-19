"""ONCVPSP stdout/stderr parsing."""

# pylint: disable=too-many-lines

from functools import cached_property
import logging
from math import copysign
from os import PathLike
from pathlib import Path
import re
from typing import Any, Literal, Sequence, TextIO

from pyoncvpsp.io.input import OncvpspInput

from .._utils import fort_float

__all__: list[str] = [
    "ONCVPSP_WARNINGS",
    "ONCVPSP_ERRORS",
    "OncvpspOutputError",
    "OncvpspOutputLine",
    "OncvpspOutputReader",
    "OncvpspTextParser",
]

LOGGER: logging.Logger = logging.getLogger(__name__)

# These are printed without stopping the program
ONCVPSP_WARNINGS: dict[str, dict[str, str]] = {
    "ICMOD1_NOT_CONVERGED": {
        "match_string": "WARNING - modcore not converged",
        "subroutine": "modcore",
        "description": "Optimization of polynomial model core charge did not converge",
    },
    "ICMOD4_NM_NOT_CONVERGED": {
        "match_string": "WARNING: not fully converged in 100 steps",
        "subroutine": "modcore3",
        "description": ("Nelder-Mead optimization of model core charge Teter parameters did not converge"),
    },
    "SRATOM_NOT_CONVERGED": {
        "match_string": "sratom: WARNING failed to converge",
        "subroutine": "sratom",
        "description": "Scalar-/non-relativistic all-electron calculation did not converge",
    },
    "RELATOM_NOT_CONVERGED": {
        "match_string": "relatom: WARNING failed to converge",
        "subroutine": "relatom",
        "description": "Relativistic all-electron calculation did not converge",
    },
    "LDIRACFB_NOT_CONVERGED": {
        "match_string": "runconfig: WARNING ldiracfb convergence ERROR n,l,kap,iter=",
        "subroutine": "run_config_r",
        "description": "",
    },
    "FR_NO_AE_SOLUTION": {
        "match_string": "run_config_r: WARNING  for AE atom,",
        "subroutine": "run_config_r",
        "description": "",
    },
    "FR_FULLY_NON_LOCAL_PS_ATOM": {
        "match_string": "run_config_r: WARNING for fully non-local PS atom",
        "subroutine": "run_config_r",
        "description": "",
    },
    "NEGATIVE_ENERGY_GHOST": {
        "match_string": "WARNING - GHOST(-)",
        "subroutine": "run_ghosts",
        "description": "A bound ghost state was detected",
    },
    "POSITIVE_ENERGY_GHOST": {
        "match_string": "WARNING - GHOST(+)",
        "subroutine": "run_ghosts",
        "description": "A highly-localized unbound ghost state was detected",
    },
    "LOCALIZED_BOUND_STATE_FOR_PROJECTOR": {
        "match_string": "WARNING wellstate: localized bound state found for n=",
        "subroutine": "wellstate[_r]",
    },
    "MODERATELY_LOCALIZED_BOUND_STATE_FOR_PROJECTOR": {
        "match_string": "WARNING wellstate: moderately localized bound state found for n=",
        "subroutine": "wellstate[_r]",
    },
    "SCATTERING_STATE_FOR_PROJECTOR": {
        "match_string": "WARNING wellstate: negative energy specified for n=",
        "subroutine": "wellstate[_r]",
    },
    "RUN_DIAG_LSCHVKBBE_NOT_CONVERGED": {
        "match_string": "run_diag: lschvkbbe ERROR",
        "subroutine": "run_diag",
        "description": "Convergence error in finding a bound state solution in a VKB pseudopotential",
    },
}

# These lead to the program immediately terminating
ONCVPSP_ERRORS: dict[str, dict[str, str]] = {
    "FR_LSCHVKBB_NOT_CONVERGED": {
        "match_string": "psatom_r: WARNING lschvkbb convergence error n,l,kap,iter=",
        "subroutine": "psatom_r",
        "description": ("Convergence error in finding a bound state solution in a VKB pseudopotential"),
    },
    "SR_LSCHVKBB_NOT_CONVERGED": {
        "match_string": "psatom: WARNING lschvkbb convergence error n,l,iter=",
        "subroutine": "psatom_nr",
        "description": (
            "Convergence error in finding a scalar-/non-relativistic " "bound state solution in a VKB pseudopotential"
        ),
    },
    "FR_WELL_STATE_NOT_CONVERGED": {
        "match_string": "ERROR wellstate_r: well potential iteration failed to converge, n=",
        "subroutine": "wellstate_r",
        "description": "Well potential solution did not converge",
    },
    "BAD_INPUT": {
        "match_string": "ERROR: test_data found",
        "subroutine": "check_data",
        "description": "Invalid input data",
    },
    "CONST_BASIS_SVD_FAILURE": {
        "match_string": "const_basis: ERTOR const_basis: dgesvd info = ",
        "subroutine": "const_basis",
        "description": "SVD failure in construction of basis for residual energy minimization",
    },
    "DPNINT_TOO_FEW_POINTS": {
        "match_string": "dpnint: interpolation ERROR, n=",
        "subroutine": "dpnint",
        "description": "Too few points for interpolation (n < polynomial order + 1)",
    },
    "DPNINT_EXTRAPOLATION": {
        "match_string": "dpnint: interpolation ERROR - out of range",
        "subroutine": "dpnint",
        "description": "Attempt to extrapolate beyond the range of data",
    },
    "LIBXC_NOT_AVAILABLE": {
        "match_string": "exc_libxc_stub: ERROR iexc = ",
        "subroutine": "exc_libxc_stub",
        "description": "Libxc was requested but is not available",
    },
    "METAGGA_NOT_IMPLEMENTED": {
        "match_string": "functrionals: ERROOR : meta GGA not coded yet ",
        "subroutine": "xc_functl_get_vxc",
        "description": "Meta-GGA functionals are not implemented",
    },
    "READ_ERROR": {
        "match_string": "Read ERROR, input data file line",
        "subroutine": "program",
        "description": "Error reading input data file",
    },
    "NO_CLASSSICAL_TURNING_POINT": {
        "match_string": "ERROR no classical turning point",
        "subroutine": "lsch*",
        "description": "No classical turning point found",
    },
    "SINGULAR_POLYM_MATRIX": {
        "match_string": "modcore:ERROR stop - singular polym matrix",
        "subroutine": "modcore",
        "description": "Singular polynomial matrix encountered in model core charge optimization",
    },
    "DGESV_INPUT_ERROR": {
        "match_string": "modcore:ERROR stop - dgesv input error",
        "subroutine": "modcore",
        "description": "Input error in LAPACK dgesv routine called during model core charge optimization",
    },
    "IRCC_NOT_FOUND": {
        "match_string": "ERROR ircc (core-valence charge crossover",
        "subroutine": "modcore*",
        "description": "Core-valence charge crossover radius not found",
    },
    "REFERENCE_ATOM_NOT_CONVERGED": {
        "match_string": "ERROR all-electron reference atom not converged",
        "subroutine": "program",
        "description": "All-electron reference atom calculation did not converge",
    },
    "FIRST_PSWF_HAS_NODE": {
        "match_string": "ERROR pspot:  first pseudo wave function has node,",
        "subroutine": "pspot",
        "description": "First pseudo wave function has a node, indicating poor psp parameters",
    },
}


def check_warning_in_line(line: str) -> str:
    """Check whether a line contains a known ONCVPSP warning.

    Args:
        line (str): ONCVPSP stdout/stderr line.

    Returns:
        str | None: Error code (key to ONCVPSP_WARNINGS) or None.
    """
    for warning_code, warning_info in ONCVPSP_WARNINGS.items():
        if warning_info["match_string"] in line:
            return warning_code
    return ""


def check_error_in_line(line: str) -> str:
    """Check whether a line contains a known ONCVPSP error.

    Args:
        line (str): ONCVPSP stdout/stderr line.

    Returns:
        str | None: Error code (key to ONCVPSP_ERRORS) or None.
    """
    for error_code, error_info in ONCVPSP_ERRORS.items():
        if error_info["match_string"] in line:
            return error_code
    return ""


def _kappa(l: int, s: float) -> int:
    if s == 0:
        return 0
    j: float = l + s
    return int((l - j) * (2 * j + 1))


class OncvpspOutputError(Exception):
    """Exception raised when an ONCVPSP error is encountered in the output file."""

    def __init__(self, line: "OncvpspOutputLine"):
        self.line: OncvpspOutputLine = line
        if line.error is None:
            raise ValueError("Line contains no error code.")
        error_info: dict[str, str] = ONCVPSP_ERRORS[line.error]
        message: str = (
            f"ONCVPSP Error encountered at line {line.number}:\n"
            f"  Subroutine: {error_info.get('subroutine', 'Unknown')}\n"
            f"  Description: {error_info.get('description', 'No description available')}\n"
            f"  '{line.text.strip()}'\n"
        )
        super().__init__(message)


class OncvpspOutputLine:
    """Text file line with associated line number and warning and error codes (if any)."""

    text: str
    number: int
    warning: str | None
    error: str | None

    def __init__(
        self,
        text: str,
        number: int,
        warning: str | None = None,
        error: str | None = None,
    ) -> None:
        self.text = text
        self.number = number
        self.warning = check_warning_in_line(text) if warning is None else warning
        self.error = check_error_in_line(text) if error is None else error

    @property
    def has_problem(self) -> bool:
        """Whether a warning or error is present on this line"""
        return self.warning is not None or self.error is not None

    def __str__(self) -> str:
        return self.text

    def __int__(self) -> int:
        return int(self.text)

    def __float__(self) -> float:
        return float(self.text)

    def __contains__(self, substring: str) -> bool:
        return substring in self.text

    def __getitem__(self, index: slice) -> str:
        return self.text[index]

    def __len__(self) -> int:
        return len(self.text)

    def __repr__(self) -> str:
        w: Literal["W"] | Literal[" "] = "W" if self.warning else " "
        e: Literal["E"] | Literal[" "] = "E" if self.error else " "
        return f"{self.number:>6d}|{w}{e}|{self.text.rstrip()}"

    def __eq__(self, other: Any) -> bool:
        if isinstance(other, str):
            return self.text == other
        if isinstance(other, OncvpspOutputLine):
            return self.text == other.text and self.warning == other.warning and self.error == other.error
        if isinstance(other, Sequence) and len(other) == 3:
            return self.text == other[0] and self.warning == other[1] and self.error == other[2]
        raise NotImplementedError(f"Cannot compare ONCVPSPOutputLine with {type(other)}")

    def startswith(self, *args, **kwargs) -> bool:
        """Calls str.startswith on the line string."""
        return self.text.startswith(*args, **kwargs)

    def split(self, *args, **kwargs) -> list[str]:
        """Calls str.split on the line string."""
        return self.text.split(*args, **kwargs)

    def strip(self, *args, **kwargs) -> "OncvpspOutputLine":
        """Calls str.strip on the line string and returns a new OncvpspOutputLine containing the stripped line."""
        return OncvpspOutputLine(
            self.text.strip(*args, **kwargs),
            self.number,
            self.warning,
            self.error,
        )

    def rstrip(self, *args, **kwargs) -> "OncvpspOutputLine":
        """Calls str.rstrip on the line string and returns a new OncvpspOutputLine containing the stripped line."""
        return OncvpspOutputLine(
            self.text.rstrip(*args, **kwargs),
            self.number,
            self.warning,
            self.error,
        )


class OncvpspOutputReader:
    """Streams an ONCVPSP output file line by line. At each line:
    - keeps track of line numbers
    - checks and stores warnings and errors with their line numbers
    - stores the lines themselves
    - yields a ONCVPSPOutputLine object
    """

    def __init__(self, file: str | PathLike | TextIO):
        if isinstance(file, (str, PathLike)):
            file_path = Path(file)
            if not file_path.is_file():
                raise FileNotFoundError(f"File not found: {file}")
            self.io = open(file_path, "r", encoding="latin=1")  # pylint: disable=consider-using-with
        else:
            self.io = file
        self.lines: list[OncvpspOutputLine] = []
        self.line_number: int = 0

    def __repr__(self) -> str:
        line_strs: list[str] = []
        for line in self.lines:
            if line.number != self.line_number:
                line_strs.append(" " + line.__repr__())
            else:
                line_strs.append(">" + line.__repr__())
        return "\n".join(line_strs)

    def __next__(self) -> OncvpspOutputLine:
        return self.readline()

    def __getitem__(self, index: int) -> OncvpspOutputLine:
        if index >= len(self.lines):
            while len(self.lines) <= index:
                _: OncvpspOutputLine = next(self)
        line: OncvpspOutputLine = self.lines[index]
        return line

    def readline(self, check_error=True, check_warning=True) -> OncvpspOutputLine:
        """Read the next line from the file or retrieve it from the lines cache.

        Args:
            check_error (bool, optional): Whether to check for errors in the line. Defaults to True.
            check_warning (bool, optional): Whether to check for warnings in the line. Defaults to True.

        Raises:
            StopIteration: End of file reached.
            OncvpspOutputError: An ONCVPSP error was encountered in the line.

        Returns:
            OncvpspOutputLine: The next line.
        """
        if self.line_number < len(self.lines):
            line: OncvpspOutputLine = self.lines[self.line_number]
            self.line_number += 1
            return line

        text: str = self.io.readline()
        self.line_number += 1

        if text == "":
            self.io.close()
            raise StopIteration

        line = OncvpspOutputLine(
            text,
            self.line_number,
            warning="" if not check_warning else None,
            error="" if not check_error else None,
        )
        self.lines.append(line)
        if line.error:
            raise OncvpspOutputError(line)

        return line

    def goto_line_number(self, number: int, **kwargs) -> None:
        """Move the reader to the specified line number.

        Args:
            number (int): Target line number (1-based).

        Raises:
            StopIteration: End of file reached before reaching the target line number.

        Returns:
            None
        """
        start_line_number: int = self.line_number
        if number <= len(self.lines):
            self.line_number = number
            return
        while self.line_number < number:
            try:
                _: OncvpspOutputLine = self.readline(**kwargs)
            except StopIteration as e:
                self.line_number = start_line_number
                raise StopIteration(f"Reached end of file before reaching line number {number}") from e
        return

    def lookahead(self) -> OncvpspOutputLine:
        """Peek at the next line without advancing the reader.

        Returns:
            OncvpspOutputLine: The next line.
        """
        position: int = self.io.tell()
        text: str = self.io.readline()
        line = OncvpspOutputLine(text, self.line_number + 1)
        self.io.seek(position)
        return line

    def goto_line_equals(self, target: str, lookback: bool = True, **kwargs) -> OncvpspOutputLine:
        """Find a line equal to `target`. If `lookback=True`, start the search
        from the first cached line. If not, start with the next unread line.

        Args:
            target (str): Target string.
            lookback (bool, optional): Whether to search the cached lines. Defaults to True.

        Raises:
            StopIteration: End of file reached before finding `target`.

        Returns:
            OncvpspOutputLine: Line equal to `target`.
        """
        start_line_number: int = self.line_number
        if not target.endswith("\n"):
            target += "\n"
        if lookback:
            for line in self.lines:
                if line == target:
                    self.line_number = line.number
                    return line
        try:
            while (line := self.readline(**kwargs)) != target:
                pass
        except StopIteration as e:
            self.line_number = start_line_number
            raise StopIteration(f'Reached end of file before finding line equal to: "{target.strip()}"') from e
        return line

    def goto_line_startswith(self, start: str, lookback: bool = True, **kwargs) -> OncvpspOutputLine:
        """Find a line starting with `start`. If `lookback=True`, start the search
        from the first cached line. If not, start with the next unread line.

        Args:
            start (str): Target string.
            lookback (bool, optional): Whether to search the cached lines. Defaults to True.

        Raises:
            StopIteration: End of file reached before finding line staring with `start`.

        Returns:
            OncvpspOutputLine: Line starting with `start`.
        """
        start_line_number: int = self.line_number
        if lookback:
            for line in self.lines:
                if line.startswith(start):
                    self.line_number = line.number
                    return line
        try:
            while not (line := self.readline(**kwargs)).startswith(start):
                pass
        except StopIteration as e:
            self.line_number = start_line_number
            raise ValueError(f'Reached end of file before finding line starting with: "{start}"') from e
        return line

    def goto_line_contains(self, substring: str, lookback: bool = True, **kwargs) -> OncvpspOutputLine:
        """Find a line contining `substring`. If `lookback=True`, start the search
        from the first cached line. If not, start with the next unread line.

        Args:
            substring (str): Target string.
            lookback (bool, optional): Whether to search the cached lines. Defaults to True.

        Raises:
            StopIteration: End of file reached before finding line containing `substring`.

        Returns:
            OncvpspOutputLine: Line containing `substring`.
        """
        start_line_number: int = self.line_number
        if lookback:
            for line in self.lines:
                if substring in line:
                    self.line_number = line.number
                    return line
        try:
            while not substring in (line := self.readline(**kwargs)):
                pass
        except StopIteration as e:
            self.line_number = start_line_number
            raise ValueError(f'Reached end of file before finding line containing: "{substring}"') from e
        return line

    def close(self) -> None:
        """Close the underlying file object."""
        self.io.close()


def restore_reader_state(func):
    """Decorator which ensures that `self.reader: OncvpspOutputReader`'s
    line number is reset to its original state after the function body.

    Args:
        func (Callable[OncvpspOutputParser, Any]): Function whose first argument is a class instance
        which has a class instance property `reader` with property `line_number: int` and bound method
        `goto_line_number(self, line_number: int)`
    """

    def wrapper(*args, **kwargs):
        self: OncvpspTextParser = args[0]
        original_line_number: int = self.reader.line_number
        result: Any = func(*args, **kwargs)
        self.reader.goto_line_number(original_line_number)
        return result

    return wrapper


class OncvpspTextParser:
    """Parses an ONCVPSP output file given an OncvpspOutputReader object."""

    def __init__(self, file: str | PathLike | TextIO | OncvpspOutputReader):
        if isinstance(file, OncvpspOutputReader):
            self.reader: OncvpspOutputReader = file
        else:
            self.reader = OncvpspOutputReader(file)

    def close(self) -> None:
        """Close the underlying OncvpspOutputReader object."""
        self.reader.close()

    @cached_property
    @restore_reader_state
    def _program_information(self) -> dict[str, str | int | float]:
        program: str = self.reader[0].strip().split()[0]
        if program == "ONCVPSP":
            relativity, _, version, date = self.reader[1].strip().split()
        else:  # program == 'METAPSP'
            relativity = "non-relativistic"
            line: OncvpspOutputLine = self.reader[2].strip()
            version_match: re.Match[str] | None = re.search(r"(\d+\.\d+\.\d+)", line.text)
            version: str = version_match.group(1) if version_match else "unknown"
            date_match: re.Match[str] | None = re.search(r"(\d{2}/\d{2}/\d{4})", line.text)
            date: str = date_match.group(1) if date_match else "01/01/1970"
        major_version, minor_version, patch_version = [int(i) for i in version.split(".")]
        day, month, year = [int(i) for i in date.split("/")]
        return {
            "program": program,
            "relativity": relativity,
            "version": version,
            "date": date,
            "major_version": major_version,
            "minor_version": minor_version,
            "patch_version": patch_version,
            "day": day,
            "month": month,
            "year": year,
        }

    @property
    def program(self) -> str:
        """Program name."""
        return self._program_information["program"]  # type: ignore

    @property
    def major_version(self) -> int:
        """Major version number."""
        return self._program_information["major_version"]  # type: ignore

    @property
    def relativistic(self) -> bool:
        """Whether the pseudopotential is fully-relativistic."""
        return self._program_information["relativity"] in ("relativistic", "fully-relativistic")  # type: ignore

    @cached_property
    @restore_reader_state
    def _atom_definition_input(self) -> dict[str, str | int | float]:
        """Example:
        # atsym  z   nc   nv     iexc    psfile
          Te 52.00    8    3       3      upf
        """
        _: OncvpspOutputLine = self.reader.goto_line_equals("# atsym  z   nc   nv     iexc    psfile\n")
        words: list[str] = self.reader.readline().strip().split()
        return {
            "atsym": words[0],
            "z": float(words[1]),
            "nc": int(words[2]),
            "nv": int(words[3]),
            "iexc": int(words[4]),
            "psfile": words[5],
        }

    @property
    def nc(self) -> int:
        """Number of core states."""
        return self._atom_definition_input["nc"]  # type: ignore

    @property
    def nv(self) -> int:
        """Number of valence states."""
        return self._atom_definition_input["nv"]  # type: ignore

    def _parse_nonrel_reference_configuration(self):
        words: list[str] = self.reader.readline().strip().split()
        n: int = int(words[0])
        l: int = int(words[1])
        f: float = float(words[2])
        eig: float = fort_float(words[3])
        return {"n": n, "l": l, "f": f, "eig": eig}

    def _parse_rel_reference_configuration(self):
        words: list[str] = self.reader.readline().strip().split()
        n: int = int(words[0])
        l: int = int(words[1])
        f: float = float(words[2])
        spins: list[float] = [+1 / 2] if l == 0 else [+1 / 2, -1 / 2]
        eigs: dict[int, float] = {}
        for idx, s in enumerate(spins):
            j: float = l + s
            kappa = int((l - j) * (2 * j + 1))
            eig: float = fort_float(words[3 + idx])
            eigs[kappa] = eig
        return {"n": n, "l": l, "f": f, "eig": eigs}

    def _parse_metapsp_reference_configuration(self):
        words: list[str] = self.reader.readline().strip().split()
        n: int = int(words[0])
        l: int = int(words[1])
        f: float = float(words[2])
        eigs: dict[str, float] = {}
        eigs["mgga"] = float(words[3])
        eigs["pbe"] = float(words[4])
        return {"n": n, "l": l, "f": f, "eigs": eigs}

    @cached_property
    @restore_reader_state
    def _reference_configuration(
        self,
    ) -> list[dict[str, int | float | dict[int, float]]]:
        if self.program == "ONCVPSP":
            if self.major_version <= 3:
                self.reader.goto_line_equals("#   n    l    f        l+1/2             l-1/2\n")
            else:
                self.reader.goto_line_equals("#   n    l    f        energy (Ha)\n")
        elif self.program == "METAPSP":
            self.reader.goto_line_equals("#   n    l    f      MGGA eval (Ha)         PBE         delta")
        else:
            raise ValueError(f"Unknown program: {self.program}")

        reference_configuration = []
        for _ in range(self.nc + self.nv):
            if not self.relativistic and self.program == "ONCVPSP":
                reference_configuration.append(self._parse_nonrel_reference_configuration())
            elif self.relativistic and self.program == "ONCVPSP":
                reference_configuration.append(self._parse_rel_reference_configuration())
            else:  # self.program == "METAPSP":
                reference_configuration.append(self._parse_metapsp_reference_configuration())
        return reference_configuration

    @cached_property
    @restore_reader_state
    def _pseudopotential_and_optimization_input(
        self,
    ):
        self.reader.goto_line_equals("# lmax\n")
        lmax = int(self.reader.readline().strip())
        self.reader.goto_line_equals("#   l,   rc,      ep,       ncon, nbas, qcut\n")
        pseudopotentials: list[dict[str, int | float]] = []
        for _ in range(lmax + 1):
            words: list[str] = self.reader.readline().strip().split()
            pseudopotentials.append(
                {
                    "l": int(words[0]),
                    "rc": float(words[1]),
                    "ep": float(words[2]),
                    "ncon": int(words[3]),
                    "nbas": int(words[4]),
                    "qcut": float(words[5]),
                }
            )
        return {"lmax": lmax, "pseudopotentials": pseudopotentials}

    @property
    def lmax(self) -> int:
        """Maximum angular momentum quantum number."""
        return self._pseudopotential_and_optimization_input["lmax"]  # type: ignore

    @property
    def rc(self) -> list[float]:
        """Core radii for each angular momentum channel."""
        return [pp["rc"] for pp in self._pseudopotential_and_optimization_input["pseudopotentials"]]  # type: ignore

    @property
    def ncon(self) -> list[int]:
        """Number of constraints for each angular momentum channel."""
        assert all(
            pp["l"] == i for i, pp in enumerate(self._pseudopotential_and_optimization_input["pseudopotentials"])
        )
        return [pp["ncon"] for pp in self._pseudopotential_and_optimization_input["pseudopotentials"]]  # type: ignore

    @cached_property
    @restore_reader_state
    def _local_potential_input(self) -> dict[str, int | float]:
        self.reader.goto_line_equals("# lloc, lpopt,  rc(5),   dvloc0\n")
        words = self.reader.readline().strip().split()
        return {
            "lloc": int(words[0]),
            "lpopt": int(words[1]),
            "rc(5)": float(words[2]),
            "dvloc0": float(words[3]),
        }

    @property
    def lloc(self) -> int:
        """Angular momentum channel of the local potential."""
        return self._local_potential_input["lloc"]  # type: ignore

    @cached_property
    @restore_reader_state
    def _vkb_projectors_input(self) -> list[dict[str, int | float]]:
        self.reader.goto_line_equals("# l, nproj, debl\n")
        vkb_projectors = []
        for _ in range(self.lmax + 1):
            words = self.reader.readline().strip().split()
            vkb_projectors.append(
                {
                    "l": int(words[0]),
                    "nproj": int(words[1]),
                    "debl": float(words[2]),
                }
            )
        return vkb_projectors

    @property
    def nproj(self) -> list[int]:
        """Number of projectors for each angular momentum channel."""
        assert all(vkb["l"] == i for i, vkb in enumerate(self._vkb_projectors_input))
        return [vkb["nproj"] for vkb in self._vkb_projectors_input]  # type: ignore

    @cached_property
    @restore_reader_state
    def _model_core_charge_input(self) -> dict[str, int | float]:
        self.reader.goto_line_equals("# icmod, fcfact, rcfact\n")
        words = self.reader.readline().strip().split()
        return {
            "icmod": int(words[0]),
            "fcfact": float(words[1]),
            "rcfact": float(words[2]),
        }

    @property
    def icmod(self) -> int:
        """Model core charge construction method."""
        return self._model_core_charge_input["icmod"]  # type: ignore

    @cached_property
    @restore_reader_state
    def _log_derivative_analysis_input(self) -> dict[str, float | None]:
        self.reader.goto_line_startswith("# epsh1, epsh2, depsh")
        words: list[str] = self.reader.readline().strip().split()
        return {
            "epsh1": float(words[0]),
            "epsh2": float(words[1]),
            "depsh": float(words[2]),
            "rxpsh": float(words[3]) if len(words) > 3 else None,
        }

    @cached_property
    @restore_reader_state
    def _output_grid_input(self) -> dict[str, float]:
        self.reader.goto_line_equals("# rlmax, drl\n")
        words: list[str] = self.reader.readline().strip().split()
        return {
            "rlmax": float(words[0]),
            "drl": float(words[1]),
        }

    @cached_property
    @restore_reader_state
    def _test_configurations_input(self) -> list[list[dict[str, int | float]]]:
        self.reader.goto_line_equals("# ncnf\n")
        ncnf = int(self.reader.readline().strip())
        self.reader.goto_line_equals("#   n    l    f\n")
        test_configurations: list[list[dict[str, int | float]]] = []
        for _ in range(ncnf):
            nvcnf = int(self.reader.readline().strip())
            test_configuration: list[dict[str, int | float]] = [{} for _ in range(nvcnf)]
            for i in range(nvcnf):
                words: list[str] = self.reader.readline().strip().split()
                test_configuration[i]["n"] = int(words[0])
                test_configuration[i]["l"] = int(words[1])
                test_configuration[i]["f"] = float(words[2])
            test_configurations.append(test_configuration)
            self.reader.readline()  # '#\n'
        return test_configurations

    @property
    def ncnf(self) -> int:
        """Number of test configurations."""
        return len(self._test_configurations_input)

    @property
    def nvcnf(self) -> list[int]:
        """Number of valence states in each test configuration."""
        return [len(cfg) for cfg in self._test_configurations_input]  # type: ignore

    @cached_property
    @restore_reader_state
    def _reference_configuration_results(self) -> dict[str, int | float]:
        self.reader.goto_line_equals("Reference configufation results\n")
        return {
            "iterations": int(self.reader.readline().strip().split()[1]),
            "all_electron_total_energy_ha": fort_float(self.reader.readline().strip().split()[4]),
        }

    def _parse_wellstate_information(self, l: int) -> dict[str, int | float]:
        line: OncvpspOutputLine = self.reader.goto_line_contains(f"Wellstate for l ={l:>3d}", lookback=False)
        words: list[str] = line.strip().split()
        l = int(words[4])
        n = int(words[7])
        if len(words) == 9:
            kappa = int(words[9])
        else:
            kappa = 0
        line = self.reader.readline()
        eigenvalue = float(line.strip().split()[-1])
        line = self.reader.readline()
        asymptotic_potential = float(line.strip().split()[-1])
        line = self.reader.readline()
        half_point_radius = float(line.strip().split()[-1])
        return {
            "l": l,
            "n": n,
            "kappa": kappa,
            "eigenvalue": eigenvalue,
            "asymptotic_potential": asymptotic_potential,
            "half_point_radius": half_point_radius,
        }

    def _parse_constraint_consistency_test(self, ncon: int, i: int) -> dict[str, list[float]]:
        self.reader.goto_line_equals("      pswf0 val/der aewf val/der     diff\n", lookback=False)
        pswf0_val_der: list[float] = [0.0] * (ncon + i)
        aewf_val_der: list[float] = [0.0] * (ncon + i)
        for j in range((ncon + i)):
            words: list[str] = self.reader.readline().strip().split()
            pswf0_val_der[j] = float(words[0])
            aewf_val_der[j] = float(words[1])
        return {
            "pswf0_val_der": pswf0_val_der,
            "aewf_val_der": aewf_val_der,
        }

    def _parse_pswf_optimization(self) -> dict[str, float]:
        line: OncvpspOutputLine = self.reader.goto_line_contains("Fraction of norm inside rc", lookback=False)
        fraction_of_norm_inside_rc = float(line.strip().split()[-1])
        line = self.reader.readline()
        # words = line.strip().split()
        # qcut = float(words[4])
        # ecut = float(words[7])
        line = self.reader.readline()
        words: list[str] = line.strip().split()
        ke_infinity = float(words[4])
        e_infinity = float(words[6])
        line = self.reader.readline()
        residual_kinetic_energy_error = fort_float(line.strip().split()[-2])
        return {
            "fraction_of_norm_inside_rc": fraction_of_norm_inside_rc,
            # "qcut": qcut,
            # "ecut": ecut,
            "ke_infinity": ke_infinity,
            "e_infinity": e_infinity,
            "residual_kinetic_energy_error": residual_kinetic_energy_error,
        }

    def _parse_kinetic_energy_consistency_test(self) -> dict[str, float]:
        self.reader.goto_line_equals("    Total kinetic energy consistency test\n", lookback=False)
        self.reader.readline()  # "Fourier integration compared to d^2/dr^2 integral"
        line: OncvpspOutputLine = self.reader.readline()
        words: list[str] = line.strip().split()
        fourier = float(words[1])
        r_space = float(words[3])
        return {
            "fourier": fourier,
            "r_space": r_space,
        }

    def _parse_potential_consistency_at_rc_test(self) -> dict[str, float]:
        self.reader.goto_line_equals("    Potential consistency test at r_c\n", lookback=False)
        line: OncvpspOutputLine = self.reader.readline()
        words: list[str] = line.strip().split()
        vpsp = float(words[1])
        vae = float(words[3])
        return {
            "vpsp": vpsp,
            "vae": vae,
        }

    def _parse_convergence_profile(self) -> dict[str, list[float]]:
        self.reader.goto_line_equals("    Energy error per electron        Cutoff\n", lookback=False)
        self.reader.readline()  # "Ha ... eV ... Ha"
        convergence_profile: dict[str, list[float]] = {
            "energy_error_per_electron_ha": [0.0] * 4,
            "energy_error_per_electron_ev": [0.0] * 4,
            "cutoff_ha": [0.0] * 4,
        }
        for i in range(4):
            words: list[str] = self.reader.readline().strip().split()
            convergence_profile["energy_error_per_electron_ha"][i] = float(words[0])
            convergence_profile["energy_error_per_electron_ev"][i] = float(words[1])
            convergence_profile["cutoff_ha"][i] = float(words[2])
        return convergence_profile

    @cached_property
    @restore_reader_state
    def _optimize_results(self):
        self.reader.goto_line_equals("Begin loop to  construct optimized pseudo wave functions\n")
        self.reader.readline()  # 'and semi-local pseudopoentials for all angular momenta'
        self.reader.readline()  # Blank line
        optimize_results = []
        for l in range(self.lmax + 1):
            wellstates = []
            for _ in range(self.nproj[l]):
                if "Wellstate" in self.reader.lookahead():
                    wellstate: dict[str, int | float] = self._parse_wellstate_information(l)
                    wellstates.append(wellstate)
                    line: OncvpspOutputLine = self.reader.readline()  # Blank line
                    assert len(line.strip()) == 0, repr(line)
                else:
                    wellstates.append({})
            if self.relativistic:
                spins: list[float] = [+1 / 2] if l == 0 else [+1 / 2, -1 / 2]
            else:
                spins = [0.0]
            for s in spins:
                k: int = _kappa(l, s)
                for i in range(self.nproj[l]):
                    line = self.reader.goto_line_startswith("Calculating optimized projector #", lookback=False)
                    assert (
                        int(line.strip().split()[-1]) == i + 1
                    ), f"Expected projector {i + 1}, got {int(line.strip().split()[-1])} from {line}"
                    line = self.reader.goto_line_startswith(" for l=   ", lookback=False)
                    assert (
                        int(line.strip().split()[-1]) == l
                    ), f"Expected l={l}, got l={int(line.strip().split()[-1])} from {line}"
                    constraint_consistency_test = self._parse_constraint_consistency_test(self.ncon[l], i)
                    pswf_optimization: dict[str, float] = self._parse_pswf_optimization()
                    kinetic_energy_consistency_test: dict[str, float] = self._parse_kinetic_energy_consistency_test()
                    potential_consistency_at_rc_test: dict[str, float] = self._parse_potential_consistency_at_rc_test()
                    convergence_profile: dict[str, list[float]] = self._parse_convergence_profile()
                    self.reader.readline()  # Blank line, set up for next "Wellstate" lookahead
                    optimize_results.append(
                        {
                            "l": l,
                            "kappa": k,
                            "i": i,
                            "wellstate_info": wellstates[i],
                            "constraint_consistency_test": constraint_consistency_test,
                            "pswf_optimization": pswf_optimization,
                            "kinetic_energy_consistency_test": kinetic_energy_consistency_test,
                            "potential_consistency_at_rc_test": potential_consistency_at_rc_test,
                            "convergence_profile": convergence_profile,
                        }
                    )
        if not self.relativistic:
            assert len(optimize_results) == sum(self.nproj)
        else:
            assert len(optimize_results) == sum(
                1 * self.nproj[l] if l == 0 else 2 * self.nproj[l] for l in range(self.lmax + 1)
            )
        return optimize_results

    def _parse_nonrel_vkb_projector_construction(self):
        construct_vkb = []
        for l in range(self.lmax + 1):
            if self.nproj[l] < 2:
                continue
            vkb_info = {}
            self.reader.readline()  # Blank line
            line = self.reader.readline()  # B matrix Hermiticity error, ll=   (l)
            assert line.strip().startswith("B matrix Hermiticity error")
            assert int(line.strip().split()[5]) == l
            vkb_info["l"] = l
            vkb_info["kappa"] = 0
            for _ in range((self.nproj[l] * (self.nproj[l] + 1) // 2) - self.nproj[l]):
                line = self.reader.readline()  # i j error
                vkb_info["i"] = int(line.strip().split()[0])
                vkb_info["j"] = int(line.strip().split()[1])
                vkb_info["error"] = fort_float(line.strip().split()[2])
            self.reader.readline()  # Blank line
            line = self.reader.readline()  # Orthonormal projector coefficients
            assert line.strip().startswith("Orthonormal projector coefficients")
            words = line.strip().split()
            vkb_info["coefficients"] = [float(word) for word in words[3:]]
            assert len(vkb_info["coefficients"]) == self.nproj[l]
            construct_vkb.append(vkb_info)
        if self.program == "METAPSP":
            self.reader.goto_line_equals("B matrix rms Hermiticity error, Npairs\n")
            line = self.reader.readline()
            words = line.strip().split()
            _ = float(words[1])  # Hermiticity RMSE
            _ = int(words[2])  # Number of pairs
        return construct_vkb

    def _parse_rel_vkb_projector_construction(self):
        construct_vkb = []
        for l in range(self.lmax + 1):
            if self.nproj[l] < 2:
                continue
            for _ in range(1 if l == 0 else 2):  # kappa
                vkb_info = {}
                self.reader.readline()  # Blank line
                line = next(self.reader)  # B matrix Hermiticity error, ll={l>:3d}  kap={kappa:>3d}
                assert line.strip().startswith("B matrix Hermiticity error")
                assert int(line.strip().split()[5]) == l
                vkb_info["l"] = l
                vkb_info["kappa"] = int(line.strip().split()[-1])
                for _ in range((self.nproj[l] * (self.nproj[l] + 1) // 2) - self.nproj[l]):
                    line = self.reader.readline()  # i j error
                    vkb_info["i"] = int(line.strip().split()[0])
                    vkb_info["j"] = int(line.strip().split()[1])
                    vkb_info["error"] = fort_float(line.strip().split()[2])
                self.reader.readline()  # Blank line
                line = self.reader.readline()  # Orthonormal projector coefficients
                assert line.strip().startswith("Orthonormal projector coefficients")
                words = line.strip().split()
                vkb_info["coefficients"] = [float(word) for word in words[3:]]
                assert len(vkb_info["coefficients"]) == self.nproj[l]
                construct_vkb.append(vkb_info)

        for l in range(1, self.lmax + 1):
            self.reader.readline()  # Blank line
            line = next(self.reader)  # Orthonormal scalar projector coefficients, l ={l>:3d}
            assert line.strip().startswith("Orthonormal scalar projector coefficients")
            assert int(line.strip().split()[-1]) == l
            line = self.reader.readline()
            construct_vkb[l]["scalar_coefficients"] = [float(word) for word in line.strip().split()]
            assert len(construct_vkb[l]["scalar_coefficients"]) == 2 * self.nproj[l]
            line = self.reader.readline()  # Blank line
            line = next(self.reader)  # Orthonormal spin-orbit projector coefficients, l ={l>:3d}
            assert line.strip().startswith("Orthonormal spin-orbit projector coefficients")
            assert int(line.strip().split()[-1]) == l
            line = self.reader.readline()
            construct_vkb[l]["spin_orbit_coefficients"] = [float(word) for word in line.strip().split()]
            assert len(construct_vkb[l]["spin_orbit_coefficients"]) == 2 * self.nproj[l]
        return construct_vkb

    @cached_property
    @restore_reader_state
    def _vkb_projector_construction(self):
        self.reader.goto_line_equals("Construct Vanderbilt / Kleinmman-Bylander projectors\n")
        if self.relativistic:
            return self._parse_rel_vkb_projector_construction()
        return self._parse_nonrel_vkb_projector_construction()

    @cached_property
    @restore_reader_state
    def _pseudo_wavefunction_solutions_mgga(self):
        if self.program != "METAPSP":
            return []
        self.reader.goto_line_equals("Psueodwavefunction solutions for first PS total energy terms\n")
        self.reader.readline()  # '  n  l       MGGA eval      <|KE,VLOC,VKB|>     delta\n'
        self.reader.readline()  # Blank line
        solutions = []
        for l in range(self.lmax + 1):
            line = self.reader.readline()
            words = line.strip().split()
            assert int(words[1]) == l
            solution = {
                "n": int(words[0]),
                "l": l,
                "mgga_eval": float(words[2]),
                "ke_vloc_vkb": float(words[3]),
            }
            solutions.append(solution)
        return solutions

    @cached_property
    @restore_reader_state
    def _model_core_correction_construction_and_analysis(self):  # pylint: disable=too-many-statements
        if self.icmod == 0 or self.icmod > 4:
            return {}
        data = {}
        self.reader.goto_line_equals("Model core correction analysis\n")
        self.reader.readline()  # '  based on d2Exc/dn_idn_j\n'
        self.reader.goto_line_equals("d2excae - all-electron derivatives\n")
        self.reader.readline()  # Blank line
        data["d2exc_ae"] = [[fort_float(i) for i in self.reader.readline().strip().split()] for _ in range(self.nv)]
        self.reader.goto_line_equals("d2excps - pseudofunction derivatives with no core correction\n")
        self.reader.readline()
        data["d2exc_ps_no_modcore"] = [
            [fort_float(i) for i in self.reader.readline().strip().split()] for _ in range(self.nv)
        ]
        line = self.reader.goto_line_startswith("rms 2nd derivative error")
        data["d2exc_rmse_no_modcore"] = fort_float(line.strip().split()[-1])
        if self.icmod == 1:
            self.reader.goto_line_equals("Polynomial model core charge\n")
        elif self.icmod == 2:
            line = self.reader.goto_line_startswith("amplitude prefactor, scale prefactor")
            words = self.reader.readline().strip().split()
            data["amp_prefactor"] = float(words[0])
            data["scale_prefactor"] = float(words[1])
        elif self.icmod == 3:
            line = self.reader.goto_line_startswith("rmatch, rhocmatch")
            words = line.strip().split()
            data["rmatch"] = float(words[2])
            data["rhocmatch"] = float(words[3])
            self.reader.goto_line_equals("Teter function model core charge with specified parameters\n")
        else:  # self.icmod == 4
            line = self.reader.goto_line_startswith("rmatch, rhocmatch")
            words = line.strip().split()
            data["rmatch"] = float(words[2])
            data["rhocmatch"] = float(words[3])
            self.reader.goto_line_equals("Coarse scan for minimum error\n")
            self.reader.readline()  # '  matrix elements: rms 2nd-derivative errors (mHa)\n'
            self.reader.readline()  # '  column index : amplitude prefactor to rhocmatch\n'
            self.reader.readline()  # '  row index : scale prefactor to rmatch\n'
            self.reader.readline()  # Blank line
            line = self.reader.readline()
            data["amp_prefactors_to_rhocmatch"] = [float(word) for word in line.strip().split()]
            data["scale_prefactors_to_rmatch"] = []
            data["d2exc_rmse_grid"] = []
            self.reader.readline()  # Blank line
            while (line := self.reader.readline()).rstrip():
                words = line.strip().split()
                data["scale_prefactors_to_rmatch"].append(float(words[0]))
                data["d2exc_rmse_grid"].append([float(word) for word in words[1:]])
            self.reader.goto_line_equals("Nelder-Mead iteration\n")
            data["nm_num_iters"] = int(self.reader.readline().strip().split()[-2])
            self.reader.goto_line_equals("amplitude prefactor, scale prefactor\n")
            line = self.reader.readline()
            words = line.strip().split()
            data["nm_opt_amp_prefactor"] = float(words[0])
            data["nm_opt_scale_prefactor"] = float(words[1])

        self.reader.goto_line_equals("d2excps - pseudofunction derivatives with core correction\n")
        self.reader.readline()
        data["d2exc_ps_with_modcore"] = [
            [fort_float(i) for i in self.reader.readline().strip().split()] for _ in range(self.nv)
        ]
        line = self.reader.goto_line_startswith("rms 2nd derivative error")
        data["d2exc_rmse_with_modcore"] = float(line.strip().split()[-1])

        return data

    @cached_property
    @restore_reader_state
    def _pseudoatom_total_energy(self) -> float:
        line = self.reader.goto_line_startswith("Pseudoatom total energy")
        return float(line[23:].strip())

    def _parse_nonrel_vkb_diagnostic_test(self) -> list[dict[str, int | float]]:
        self.reader.goto_line_equals("   l    rcore       rmatch      e in        delta e   norm test  slope test\n")
        diagnostics = []
        for l in range(self.lmax + 1):
            self.reader.readline()  # Blank line
            for i in range(self.nproj[l]):
                line = self.reader.readline()
                if line.warning:  # Likely "run_diag: lschvkbbe ERROR ..."
                    line = self.reader.readline()
                words = line.strip().split()
                assert int(words[0]) == l
                # assert abs(float(words[1]) - self.rc[l]) < 1e-5  # This isn't true -- rc is locked to the grid
                vkb_diag_info = {
                    "l": int(words[0]),
                    "i": i,
                    "kappa": 0,
                    "rcore": float(words[1]),
                    "rmatch": float(words[2]),
                    "e_in": float(words[3]),
                    "delta_e": float(words[4]),
                    "norm_test": float(words[5]),
                    "slope_test": float(words[6]),
                }
                diagnostics.append(vkb_diag_info)
        return diagnostics

    def _parse_rel_vkb_diagnostic_test(self) -> list[dict[str, int | float]]:
        self.reader.goto_line_equals(
            "   l  kap   rcore       rmatch      e in        delta e    norm test   slope test\n"
        )
        diagnostics = []
        for l in range(self.lmax + 1):
            self.reader.readline()  # Blank line
            for i in range(self.nproj[l]):
                for _ in range(1 if l == 0 else 2):
                    line = self.reader.readline()
                    words = line.strip().split()
                    assert int(words[0]) == l
                    # This won't be rc, but a nearby point on the mesh
                    # assert abs(float(words[2]) - self.rc[l]) < 1e-5
                    vkb_diag_info = {
                        "l": int(words[0]),
                        "i": i,
                        "kappa": int(words[1]),
                        "rcore": float(words[2]),
                        "rmatch": float(words[3]),
                        "e_in": float(words[4]),
                        "delta_e": float(words[5]),
                        "norm_test": float(words[6]),
                        "slope_test": float(words[7]),
                    }
                    diagnostics.append(vkb_diag_info)
        return diagnostics

    def _parse_rel_vkb_diagnostic_tests(self):
        self.reader.goto_line_equals("  relativistic with J = L +/- 1/2 projectors\n")
        j_diagnostics = self._parse_rel_vkb_diagnostic_test()
        self.reader.goto_line_equals("  relativistic with scalar and spin-orbit non-local projectors\n")
        scalar_so_diagnostics = self._parse_rel_vkb_diagnostic_test()
        return {
            "j_diagnostics": j_diagnostics,
            "scalar_so_diagnostics": scalar_so_diagnostics,
        }

    @cached_property
    @restore_reader_state
    def _vkb_diagnostic_tests(
        self,
    ) -> list[dict[str, int | float]] | dict[str, list[dict[str, int | float]]]:
        if self.relativistic:
            diagnostics = self._parse_rel_vkb_diagnostic_tests()
        else:
            diagnostics = self._parse_nonrel_vkb_diagnostic_test()
        return diagnostics

    def _parse_ghosts_test(self) -> dict[str, list[dict[str, int | float]]]:
        self.reader.goto_line_equals(" Testing for bound ghosts\n", lookback=False)
        self.reader.readline()  # '   n   l    E NL Schr. Eq
        bound_ghosts = []
        while not (line := self.reader.readline()).startswith(" Testing"):
            if len(line.rstrip()) == 0:
                continue
            bound_ghosts.append(
                {
                    "n": int(line[0:4]) if line[0:4].strip() else None,
                    "l": int(line[4:8]),
                    "e_nl_schrod_eq": (float(line[8:24]) if line[8:24].strip() else None),
                    "e_basis_diag": float(line[24:40]),
                    "e_cutoff": float(line[40:50]),
                    "warning": "WARNING" in line,
                }
            )
        self.reader.goto_line_equals(" Testing for highly-localized positive-energy ghosts\n")
        self.reader.readline()  # '       l    <radius>/rc     E Basis Diag.   E Cutoff'
        unbound_ghosts = []
        while (line := self.reader.readline()) not in (
            "Test configuration 0\n",
            "Ghost test for J = L - 1/2\n",
        ):
            words = line.strip().split()
            if len(words) == 0:
                continue
            unbound_ghosts.append(
                {
                    "l": int(words[0]),
                    "radius_over_rc": float(words[1]),
                    "e_basis_diag": float(words[2]),
                    "e_cutoff": float(words[3]),
                    "warning": "WARNING" in line,
                }
            )
        self.reader.goto_line_number(line.number - 1)
        return {
            "bound_ghosts": bound_ghosts,
            "unbound_ghosts": unbound_ghosts,
        }

    @cached_property
    @restore_reader_state
    def _ghosts_tests(
        self,
    ) -> dict[str, list[dict[str, int | float]]] | dict[str, dict[str, list[dict[str, int | float]]]]:
        if self.relativistic:
            self.reader.goto_line_equals("Ghost test for J = L + 1/2\n")
            ghosts_s_up = self._parse_ghosts_test()
            self.reader.goto_line_equals("Ghost test for J = L - 1/2\n")
            ghosts_s_dw = self._parse_ghosts_test()
            return {
                "ghosts_s_up": ghosts_s_up,
                "ghosts_s_dw": ghosts_s_dw,
            }
        return self._parse_ghosts_test()

    def _parse_rel_test_configuration_result(self, test_index):  # pylint: disable=too-many-locals
        self.reader.goto_line_equals(f"Test configuration {test_index}\n", lookback=False)
        self.reader.readline()  # Blank line
        self.reader.readline()  # '   n   l     f        eae           eps        diff'
        nvcnf_i = self.nv if test_index == 0 else self.nvcnf[test_index - 1]
        states = []
        for i_state in range(self.nc + nvcnf_i):
            if i_state < self.nc or test_index == 0:
                n = self._reference_configuration[i_state]["n"]
                l = self._reference_configuration[i_state]["l"]
                # f = self._reference_configuration[i_state]["f"]
            else:
                n = self._test_configurations_input[test_index - 1][i_state - self.nc]["n"]
                l = self._test_configurations_input[test_index - 1][i_state - self.nc]["l"]
                # f = self._test_configurations_input[test_index][i_state - self.nc]["f"]
            n_kappa = 1 if l == 0 else 2
            for _ in range(n_kappa):
                line = self.reader.readline().rstrip()
                state = {
                    "n": int(line[0:4]),
                    "l": int(line[4:8]),
                    "kappa": int(line[8:12]),
                    "f": float(line[12:20]),
                    "e_ae": float(line[20:34]),
                    "e_ps": None,
                }
                assert state["n"] == n, f"Test {test_index}, state {i_state} expected n={n}, got n={state['n']}"
                assert state["l"] == l, f"Test {test_index}, state {i_state} Expected l={l}, got l={state['l']}"
                # assert abs(state["f"] - f) < 1e-12
            if i_state >= self.nc:
                state["e_ps"] = float(line[34:47])
            states.append(state)
        self.reader.readline()  # Blank line
        self.reader.readline()  # '     Total energies and differences'
        words = self.reader.readline().strip().split()
        ae_ref, ae_tst = fort_float(words[1]), fort_float(words[3])
        words = self.reader.readline().strip().split()
        ps_ref, ps_tst = fort_float(words[1]), fort_float(words[3])
        psp_excitation_error = fort_float(self.reader.readline().strip().split()[-1])
        return {
            "states": states,
            "energies": {
                "ae_ref": ae_ref,
                "ae_tst": ae_tst,
                "ps_ref": ps_ref,
                "ps_tst": ps_tst,
                "psp_excitation_error": psp_excitation_error,
            },
        }

    def _parse_nonrel_test_configuration_result(self, test_index):  # pylint: disable=too-many-locals
        self.reader.goto_line_equals(f"Test configuration {test_index}\n", lookback=False)
        self.reader.readline()  # Blank line
        self.reader.readline()  # '   n   l     f        eae           eps        diff'
        nvcnf_i = self.nv if test_index == 0 else self.nvcnf[test_index - 1]
        states = []
        for i_state in range(self.nc + nvcnf_i):
            if i_state < self.nc or test_index == 0:
                n = self._reference_configuration[i_state]["n"]
                l = self._reference_configuration[i_state]["l"]
                f = self._reference_configuration[i_state]["f"]
            else:
                n = self._test_configurations_input[test_index - 1][i_state - self.nc]["n"]
                l = self._test_configurations_input[test_index - 1][i_state - self.nc]["l"]
                f = self._test_configurations_input[test_index - 1][i_state - self.nc]["f"]
            line = self.reader.readline().rstrip()
            state = {
                "n": int(line[0:4]),
                "l": int(line[4:8]),
                "kappa": 0,
                "f": float(line[8:16]),
                "e_ae": float(line[16:30]),
                "e_ps": None,
            }
            assert state["n"] == n, f"Test {test_index}, state {i_state} expected n={n}, got n={state['n']}"
            assert state["l"] == l, f"Test {test_index}, state {i_state} Expected l={l}, got l={state['l']}"
            assert abs(state["f"] - f) < 1e-12, f"Test {test_index}, state {i_state} expected f={f}, got f={state['f']}"
            if i_state >= self.nc:
                state["e_ps"] = float(line[30:44])
            states.append(state)
        self.reader.readline()  # Blank line
        self.reader.readline()  # '     Total energies and differences'
        words = self.reader.readline().strip().split()
        ae_ref, ae_tst = fort_float(words[1]), fort_float(words[3])
        words = self.reader.readline().strip().split()
        ps_ref, ps_tst = fort_float(words[1]), fort_float(words[3])
        psp_excitation_error = fort_float(self.reader.readline().strip().split()[-1])
        return {
            "states": states,
            "energies": {
                "ae_ref": ae_ref,
                "ae_tst": ae_tst,
                "ps_ref": ps_ref,
                "ps_tst": ps_tst,
                "psp_excitation_error": psp_excitation_error,
            },
        }

    @cached_property
    @restore_reader_state
    def _test_configuration_results(self):
        if self.relativistic:
            return [self._parse_rel_test_configuration_result(i) for i in range(self.ncnf)]
        return [self._parse_nonrel_test_configuration_result(i) for i in range(self.ncnf)]

    @cached_property
    @restore_reader_state
    def generation_data(self):
        """Data related to the pseudopotential generation process, including:
        - Program information
        - Atom definition
        - Reference configuration input parameters and eigenvalues
        - Pseudopotential and optimization input parameters
        - Local potential input parmaetrs
        - VKB projector input parameters
        - Model core charge input parameters
        - Log-derivative analysis input parameters
        - Output grid input parameters
        - Test configurations input parameters
        - Reference configuration results (number of iterations, total energy, etc.)
        - Optimize results (wellstate information, optimization tests, convergence profiles)
        - VKB projector construction data (Hermiticity errors, coefficients)
        - Model core correction construction and analysis data (2nd derivative matrices, RMSEs)
        - Pseudoatom total energy
        - VKB diagnostic tests (rcore, rmatch, energy in, delta e, norm test, slope test)
        - Ghosts tests (bound and unbound ghosts)
        - Test configuration results (eigenvalues, total energies, excitation errors)

        Returns:
            dict: Generation data
        """
        return {
            "program_information": self._program_information,
            "atom_definition": self._atom_definition_input,
            "reference_configuration": self._reference_configuration,
            "pseudopotential_and_optimization": self._pseudopotential_and_optimization_input,
            "local_potential": self._local_potential_input,
            "vkb_projectors": self._vkb_projectors_input,
            "model_core_charge": self._model_core_charge_input,
            "log_derivative_analysis": self._log_derivative_analysis_input,
            "output_grid": self._output_grid_input,
            "test_configurations": self._test_configurations_input,
            "reference_configuration_results": self._reference_configuration_results,
            "optimize_results": self._optimize_results,
            "vkb_projector_construction": self._vkb_projector_construction,
            "model_core_correction_construction_and_analysis": self._model_core_correction_construction_and_analysis,
            "pseudoatom_total_energy": self._pseudoatom_total_energy,
            "vkb_diagnostic_tests": self._vkb_diagnostic_tests,
            "ghosts_tests": self._ghosts_tests,
            "test_configuration_results": self._test_configuration_results,
        }

    @cached_property
    @restore_reader_state
    def _semilocal_potentials_plot_data(self):
        self.reader.goto_line_equals(
            " radii, charge, pseudopotentials (ll=0, 1, lmax)\n",
            check_warning=True,
            check_error=True,
        )
        self.reader.readline()  # Blank line
        data = {
            "semilocal_potentials": [{"l": l, "f": []} for l in range(self.lmax + 1)],
            "rho_val": {"f": []},
            "local_potential": {},
        }
        r_list = []
        while (line := self.reader.readline()).startswith("!p"):
            words = line.strip().split()
            r_list.append(float(words[1]))
            data["rho_val"]["f"].append(float(words[2]))
            for l in range(self.lmax + 1):
                data["semilocal_potentials"][l]["f"].append(float(words[3 + l]))
        data["rho_val"]["r"] = r_list
        for l in range(self.lmax + 1):
            data["semilocal_potentials"][l]["r"] = r_list
        if line.startswith(" !L"):  # lloc == 4
            assert self.lloc == 4
            self.reader.goto_line_number(
                line.number - 1,
            )
            data["local_potential"]["f"] = []
            data["local_potential"]["r"] = []
            while (line := self.reader.readline()).startswith(" !L"):
                words = line.strip().split()
                data["local_potential"]["r"].append(float(words[1]))
                data["local_potential"]["f"].append(float(words[2]))
        return data

    @cached_property
    @restore_reader_state
    def _charge_density_plot_data(self):
        self.reader.goto_line_equals(
            " radii, charge, core charge, model core charge\n",
            check_warning=True,
            check_error=True,
        )
        self.reader.readline()  # Blank line
        data = {
            "rho_val": {"f": []},
            "rho_core": {"f": []},
            "rho_model_core": {"f": []},
        }
        r_list = []
        while (line := self.reader.readline()).startswith("!r"):
            words = line.strip().split()
            r_list.append(float(words[1]))
            data["rho_val"]["f"].append(float(words[2]))
            data["rho_core"]["f"].append(float(words[3]))
            data["rho_model_core"]["f"].append(float(words[4]))
        for value in data.values():
            value["r"] = r_list
        return data

    @cached_property
    @restore_reader_state
    def _kinetic_energy_density_plot_data(self):
        if self.program != "METAPSP":
            return {}
        self.reader.goto_line_equals(
            " Metagga taups and taumodps \n",
        )
        self.reader.readline()  # Blank line
        data = {
            "tau": {"tau_ps": {"f": []}, "tau_mod_ps": {"f": []}},
            "vtau": {"vtau_ps": {"f": []}, "vtau_mod_ps": {"f": []}},
        }
        r_list = []
        while (line := self.reader.readline()).startswith("!t"):
            words = line.strip().split()
            r_list.append(float(words[1]))
            data["tau"]["tau_ps"]["f"].append(float(words[2]))
            data["tau"]["tau_mod_ps"]["f"].append(float(words[3]))
        for value in data["tau"].values():
            value["r"] = r_list

        line = self.reader.goto_line_startswith(
            "!vt",
            lookback=False,
        )
        self.reader.goto_line_number(
            line.number - 1,
        )
        r_list = []
        while (line := self.reader.readline()).startswith("!vt"):
            words = line.strip().split()
            r_list.append(float(words[1]))
            data["vtau"]["vtau_ps"]["f"].append(float(words[2]))
            data["vtau"]["vtau_mod_ps"]["f"].append(float(words[3]))
        for value in data["vtau"].values():
            value["r"] = r_list
        return data

    def _parse_nonrel_wavefunctions_and_vkb_plot_data(self):
        wavefunctions = []
        vkb_projectors = []
        for l in range(self.lmax + 1):
            # Parse wavefunctions
            for i in range(self.nproj[l]):
                line = self.reader.goto_line_contains(
                    "all-electron wave function, pseudo w-f",
                    lookback=False,
                )  # There can be warnings here about well states
                words = line.strip().split()
                is_bound = "scattering" not in line
                if is_bound:
                    n = int(words[1].strip(","))
                    assert int(words[3].strip(",")) == l
                else:
                    n = 0  # For scattering states, n is not defined
                    assert int(words[2].strip(",")) == i + 1
                    assert int(words[4].strip(",")) == l
                ae_wf = {
                    "n": n,
                    "l": l,
                    "kappa": 0,
                    "i": i,
                    "is_bound": is_bound,
                    "ae_ps": "ae",
                    "f": [],
                }
                ps_wf = {
                    "n": n,
                    "l": l,
                    "i": i,
                    "kappa": 0,
                    "is_bound": is_bound,
                    "ae_ps": "ps",
                    "f": [],
                }
                self.reader.readline()  # Blank line
                r_list = []
                while (line := self.reader.readline()).startswith("&"):
                    words = line.strip().split()
                    r_list.append(float(words[2]))
                    ae_wf["f"].append(float(words[3]))
                    ps_wf["f"].append(float(words[4]))
                ae_wf["r"] = r_list
                wavefunctions.append(ae_wf)
                ps_wf["r"] = r_list
                wavefunctions.append(ps_wf)
            # Parse VKB projectors
            self.reader.readline()  # Blank line
            projs_l = [{"l": l, "i": i, "kappa": 0, "f": []} for i in range(self.nproj[l])]
            r_list = []
            while (line := self.reader.readline()).startswith("!J"):
                words = line.strip().split()
                r_list.append(float(words[2]))
                for i in range(self.nproj[l]):
                    assert int(words[1]) == l
                    projs_l[i]["f"].append(float(words[3 + i]))
            for proj in projs_l:
                proj["r"] = r_list
            vkb_projectors.extend(projs_l)
        return {
            "wavefunctions": wavefunctions,
            "vkb_projectors": vkb_projectors,
        }

    def _parse_rel_wavefunctions_and_vkb_plot_data(self):  # pylint: disable=too-many-locals
        wavefunctions = []
        vkb_projectors = []
        for l in range(self.lmax + 1):
            for s in (+1 / 2,) if l == 0 else (+1 / 2, -1 / 2):
                kappa = _kappa(l=l, s=s)
                for i in range(self.nproj[l]):
                    line = self.reader.goto_line_contains(
                        "all-electron wave function, pseudo w-f", lookback=False
                    )  # There can be warnings here about well states
                    words = line.strip().split()
                    is_bound = "scattering" not in line
                    if is_bound:
                        n = int(words[1].strip(","))
                        assert int(words[3].strip(",")) == l
                        assert int(line[17:19]) == kappa
                    else:
                        n = 0
                        assert int(words[2].strip(",")) == i + 1
                        assert int(words[4].strip(",")) == l
                        assert int(line[32:34]) == kappa
                    ae_wf = {
                        "n": n,
                        "l": l,
                        "kappa": kappa,
                        "i": i,
                        "is_bound": is_bound,
                        "ae_ps": "ae",
                        "f": [],
                    }
                    ps_wf = {
                        "n": n,
                        "l": l,
                        "kappa": kappa,
                        "i": i,
                        "is_bound": is_bound,
                        "ae_ps": "ps",
                        "f": [],
                    }
                    self.reader.readline()  # Blank line
                    r_list = []
                    while (line := self.reader.readline()).startswith("&"):
                        words = line.strip().split()
                        r_list.append(float(words[2]))
                        ae_wf["f"].append(float(words[3]))
                        ps_wf["f"].append(float(words[4]))
                    ae_wf["r"] = r_list
                    wavefunctions.append(ae_wf)
                    ps_wf["r"] = r_list
                    wavefunctions.append(ps_wf)
                # Parse VKB projectors
                self.reader.readline()  # Blank line
                projs = [{"l": l, "i": i, "kappa": kappa, "r": [], "f": []} for i in range(self.nproj[l])]
                r_list = []
                while (line := self.reader.readline()).startswith("!J"):
                    words = line.strip().split()
                    r_list.append(float(words[2]))
                    for i in range(self.nproj[l]):
                        assert int(words[1]) == copysign(l, kappa)
                        projs[i]["f"].append(float(words[3 + i]))
                for proj in projs:
                    proj["r"] = r_list
                vkb_projectors.extend(projs)
        return {
            "wavefunctions": wavefunctions,
            "vkb_projectors": vkb_projectors,
        }

    @cached_property
    @restore_reader_state
    def _wavefunction_and_vkb_projector_plot_data(self):
        if self.relativistic:
            return self._parse_rel_wavefunctions_and_vkb_plot_data()
        return self._parse_nonrel_wavefunctions_and_vkb_plot_data()

    @cached_property
    @restore_reader_state
    def _convergence_profiles_plot_data(self):
        self.reader.goto_line_startswith(
            "convergence profiles, (ll=0,lmax)",
        )
        self.reader.readline()  # lmax line
        profiles = [{"l": l, "ecut_ha": [], "error_ha": []} for l in range(self.lmax + 1)]
        while (line := self.reader.readline()).startswith("!C"):
            words = line.strip().split()
            l = int(words[1])
            assert profiles[l]["l"] == l
            profiles[l]["ecut_ha"].append(float(words[2]))
            profiles[l]["error_ha"].append(float(words[3]))
        return profiles

    def _parse_nonrel_log_derivative_analysis_plot_data(self):
        log_ders = []
        for l in range(self.lmax + 1):
            self.reader.goto_line_equals(
                f"log derivativve data for plotting, l={l:>2d}\n",
                lookback=False,
                check_warning=True,
                check_error=True,
            )
            r = float(self.reader.readline().strip().split()[-1])
            self.reader.readline()  # 'l, energy, all-electron, pseudopotential
            self.reader.readline()  # Blank line
            log_der_ae = {"l": l, "kappa": 0, "r": r, "ae_ps": "ae", "f": []}
            log_der_ps = {"l": l, "kappa": 0, "r": r, "ae_ps": "ps", "f": []}
            e_list = []
            while (line := self.reader.readline()).startswith("!"):
                words = line.strip().split()
                assert l == int(words[1])
                e_list.append(float(words[2]))
                log_der_ae["f"].append(float(words[3]))
                log_der_ps["f"].append(float(words[4]))
            log_der_ae["e"] = e_list
            log_ders.append(log_der_ae)
            log_der_ps["e"] = e_list
            log_ders.append(log_der_ps)
        return log_ders

    def _parse_rel_log_derivative_analysis_plot_data(self):
        log_ders = []
        for l in range(self.lmax + 1):
            for s in (+1 / 2,) if l == 0 else (+1 / 2, -1 / 2):
                kappa = _kappa(l=l, s=s)
                self.reader.goto_line_equals(
                    f"log derivativve data for plotting, l={l:>2d}\n",
                    lookback=False,
                    check_warning=True,
                    check_error=True,
                )
                line = self.reader.readline()  # 'atan(r * ((d psi(r)/dr)/psi(r))), r={r:6.2f}\n'
                r = float(line.strip().split()[-1])
                self.reader.readline()  # 'l, energy, all-electron, pseudopotential\n'
                self.reader.readline()  # Blank line
                log_der_ae = {"l": l, "kappa": 0, "r": r, "ae_ps": "ae", "f": []}
                log_der_ps = {"l": l, "kappa": 0, "r": r, "ae_ps": "ps", "f": []}
                e_list = []
                while (line := self.reader.readline()).startswith("!"):
                    words = line.strip().split()
                    assert copysign(l, kappa) == int(words[1])
                    e_list.append(float(words[2]))
                    log_der_ae["f"].append(float(words[3]))
                    log_der_ps["f"].append(float(words[4]))
                log_der_ae["e"] = e_list
                log_ders.append(log_der_ae)
                log_der_ps["e"] = e_list
                log_ders.append(log_der_ps)
        return log_ders

    @cached_property
    @restore_reader_state
    def _log_derivative_analysis_plot_data(self):
        line = self.reader.goto_line_equals(f"log derivativve data for plotting, l={0:>2d}\n")
        self.reader.goto_line_number(line.number - 1)
        if self.relativistic:
            return self._parse_rel_log_derivative_analysis_plot_data()
        return self._parse_nonrel_log_derivative_analysis_plot_data()

    @cached_property
    @restore_reader_state
    def plot_data(self):
        """Get the data printed to the output file by the `run_plot` subroutine of ONCVPSP.
        This includes:
            - Semilocal potentials (including the local potential)
            - Charge densities (all-electron core, pseudo-valence, model core)
            - Kinetic energy densities (for METAPSP only)
            - AE and PS bound and scattering wavefunctions and VKB projectors
            - Kinetic energy residual convergence profiles
            - Phase-shift-like log derivatives

        Returns:
            dict[str, Any]: A dictionary containing the plot data.
        """
        return {
            "semilocal_potentials": self._semilocal_potentials_plot_data,
            "charge_density": self._charge_density_plot_data,
            "kinetic_energy_density": self._kinetic_energy_density_plot_data,
            "wavefunctions_and_vkb_projectors": self._wavefunction_and_vkb_projector_plot_data,
            "convergence_profiles": self._convergence_profiles_plot_data,
            "log_derivative_analysis": self._log_derivative_analysis_plot_data,
        }

    @cached_property
    @restore_reader_state
    def all_data(self):
        """Get all the data parsed from the ONCVPSP text output file.

        Returns:
            dict[str, Any]: Generation data and plot data.
        """
        return {
            **self.generation_data,
            **self.plot_data,
        }

    def parse(self):
        """Parse the file and initialize all cached properties.

        Returns:
            "OncvpspTextOutput": This object.
        """
        # This will populate the cached propertoes
        self.all_data  # pylint: disable=pointless-statement
        try:
            self.reader.close()
        except Exception as e:  # pylint: disable=broad-exception-caught
            LOGGER.warning("Failed to close reader: %s", e)
        return self

    @property
    def input(self) -> OncvpspInput:
        """OncvpspInput: The input parameters used to generate this pseudopotential."""
        at_def: dict[str, str | int | float] = self._atom_definition_input
        ref_conf: list[dict[str, int | float | dict[int, float]]] = self._reference_configuration
        psp_opt: dict[str, int | list[dict[str, int | float]]] = self._pseudopotential_and_optimization_input
        vkb_projs: list[dict[str, int | float]] = self._vkb_projectors_input
        loc_pot: dict[str, int | float] = self._local_potential_input
        mod_core: dict[str, int | float] = self._model_core_charge_input
        out_grid: dict[str, float] = self._output_grid_input
        log_der: dict[str, float | None] = self._log_derivative_analysis_input
        test_confs: list[list[dict[str, int | float]]] = self._test_configurations_input
        input_data = {
            "oncvpsp": {
                "atsym": at_def["atsym"],
                "z": at_def["z"],
                "nc": at_def["nc"],
                "nv": at_def["nv"],
                "iexc": at_def["iexc"],
            },
            "reference_configuration": {
                "n": [s["n"] for s in ref_conf],
                "l": [s["l"] for s in ref_conf],
                "f": [s["f"] for s in ref_conf],
            },
            "pseudopotentials": {
                "lmax": psp_opt["lmax"],
                "l": [psp["l"] for psp in psp_opt["pseudopotentials"]],  # type: ignore
                "rc": [psp["rc"] for psp in psp_opt["pseudopotentials"]],  # type: ignore
                "ep": [psp["ep"] for psp in psp_opt["pseudopotentials"]],  # type: ignore
                "ncon": [psp["ncon"] for psp in psp_opt["pseudopotentials"]],  # type: ignore
                "nbas": [psp["nbas"] for psp in psp_opt["pseudopotentials"]],  # type: ignore
                "qcut": [psp["qcut"] for psp in psp_opt["pseudopotentials"]],  # type: ignore
            },
            "local_potential": {
                "lloc": loc_pot["lloc"],
                "lpopt": loc_pot["lpopt"],
                "rcloc": loc_pot["rc(5)"],
                "dvloc0": loc_pot["dvloc0"],
            },
            "vkb_projectors": {
                "l": [proj["l"] for proj in vkb_projs],
                "nproj": [proj["nproj"] for proj in vkb_projs],
                "debl": [proj["debl"] for proj in vkb_projs],
            },
            "model_core_charge": {
                "icmod": mod_core["icmod"],
                "fcfact": mod_core["fcfact"],
                "rcfact": mod_core["rcfact"],
            },
            "log_derivative_analysis": {
                "epsh1": log_der["epsh1"],
                "epsh2": log_der["epsh2"],
                "depsh": log_der["depsh"],
                "rxpsh": log_der.get("rxpsh"),
            },
            "test_configurations": [
                {
                    "n": [s["n"] for s in config],
                    "l": [s["l"] for s in config],
                    "f": [s["f"] for s in config],
                }
                for config in test_confs
            ],
            "linear_mesh": {
                "drl": out_grid["drl"],
                "rlmax": out_grid["rlmax"],
            },
            "pp_output": {
                "psfile": at_def["psfile"],
            },
        }
        return OncvpspInput(**input_data)
