#!/usr/bin/env python3
# %%
import argparse
import logging
import pathlib as pl
import sys
import typing as ty
from dataclasses import dataclass

import toml

logging.basicConfig(level=logging.ERROR)

ATOL: float = 1e-12


def fort_float(s: str) -> float:
    """Convert a Fortran-style float string to a Python float.

    Args:
        s: String to convert.

    Returns:
        Converted float.
    """
    return float(s.lower().replace("d", "e"))


@dataclass
class StateConfiguration:
    """Configuration of an atomic state."""

    "Principal quantum number"
    n: int
    "Angular quantum number"
    l: int
    "Occupation number"
    f: float

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, StateConfiguration):
            raise NotImplementedError(
                f"Cannot compare StateConfiguration with {type(other)}"
            )
        if self.n != other.n:
            logging.debug("n differs: %d != %d", self.n, other.n)
            return False
        if self.l != other.l:
            logging.debug("l differs: %d != %d", self.l, other.l)
            return False
        if abs(self.f - other.f) >= ATOL:
            logging.debug("f differs: %f != %f", self.f, other.f)
            return False
        return True

    @classmethod
    def from_str(cls, s: str) -> "StateConfiguration":
        """Create a StateConfiguration from a string.

        Args:
            s: String in the format "n l f".

        Returns:
            StateConfiguration instance.
        """
        parts = s.split()
        if len(parts) < 3:
            raise ValueError(f"Invalid state configuration: {s}")
        if len(parts) > 3:
            logging.warning("Ignoring extra parts in state configuration: %s", s)
        n, l, f = int(parts[0]), int(parts[1]), fort_float(parts[2])
        return cls(n=n, l=l, f=f)

    @classmethod
    def to_str(cls, state: "StateConfiguration", with_comment: bool = False) -> str:
        """Convert a StateConfiguration to a string.

        Args:
            state: StateConfiguration instance.

        Returns:
            String in the format "n l f".
        """
        s = f"{state.n:2d} {state.l:2d} {state.f:10.5f}"
        if with_comment:
            s = f"# {'n':>2} {'l':>2} {'f':>10}\n" + s
        return s

    def to_dict(self) -> dict[str, ty.Any]:
        """Convert a StateConfiguration to a dictionary.

        Returns:
            Dictionary with keys "n", "l", and "f".
        """
        return {"n": self.n, "l": self.l, "f": self.f}


@dataclass
class PseudopotentialConfiguration:
    """Configuration of a pseudopotential."""

    "Angular momentum"
    l: int
    "Core radii"
    rc: float
    "Energy"
    ep: float
    "Number of matching constraints"
    ncon: int
    "Number of basis functions"
    nbas: int
    "Cutoff wavevectors"
    qcut: float

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, PseudopotentialConfiguration):
            raise NotImplementedError(
                f"Cannot compare PseudopotentialConfiguration with {type(other)}"
            )
        if self.l != other.l:
            logging.debug("l differs: %d != %d", self.l, other.l)
            return False
        if abs(self.rc - other.rc) >= ATOL:
            logging.debug("rc differs: %f != %f", self.rc, other.rc)
            return False
        if abs(self.ep - other.ep) >= ATOL:
            logging.debug("ep differs: %f != %f", self.ep, other.ep)
            return False
        if self.ncon != other.ncon:
            logging.debug("ncon differs: %d != %d", self.ncon, other.ncon)
            return False
        if self.nbas != other.nbas:
            logging.debug("nbas differs: %d != %d", self.nbas, other.nbas)
            return False
        if abs(self.qcut - other.qcut) >= ATOL:
            logging.debug("qcut differs: %f != %f", self.qcut, other.qcut)
            return False
        return True

    @classmethod
    def from_str(cls, s: str) -> "PseudopotentialConfiguration":
        """Create a PseudopotentialConfiguration from a string.

        Args:
            s: String in the format "l rc ep ncon nbas qcut"

        Returns:
            PseudopotentialConfiguration instance.
        """
        parts = s.split()
        if len(parts) != 6:
            raise ValueError(f"Invalid pseudopotential configuration: {s}")
        l, rc, ep, ncon, nbas, qcut = (
            int(parts[0]),
            fort_float(parts[1]),
            fort_float(parts[2]),
            int(parts[3]),
            int(parts[4]),
            fort_float(parts[5]),
        )
        return cls(l=l, rc=rc, ep=ep, ncon=ncon, nbas=nbas, qcut=qcut)

    @classmethod
    def to_str(
        cls, pp: "PseudopotentialConfiguration", with_comment: bool = False
    ) -> str:
        """Convert a PseudopotentialConfiguration to a string.

        Args:
            pp: PseudopotentialConfiguration instance.

        Returns:
            String in the format "l rc ep ncon nbas qcut".
        """
        s = f"{pp.l:2d} {pp.rc:10.5f} {pp.ep:10.5f} {pp.ncon:2d} {pp.nbas:2d} {pp.qcut:10.5f}"
        if with_comment:
            s = (
                f"# {'l':>2} {'rc':>10} {'ep':>10} {'ncon':>2} {'nbas':>2} {'qcut':>10}\n"
                + s
            )
        return s

    def to_dict(self) -> dict[str, ty.Any]:
        """Convert a PseudopotentialConfiguration to a dictionary.

        Returns:
            Dictionary with keys "l", "rc", "ep", "ncon", "nbas", and "qcut".
        """
        return {
            "l": self.l,
            "rc": self.rc,
            "ep": self.ep,
            "ncon": self.ncon,
            "nbas": self.nbas,
            "qcut": self.qcut,
        }


@dataclass
class LocalPotentialConfiguration:
    """Configuration of the local potential."""

    "Local potential angular momentum"
    lloc: int
    "Local potential optimization level"
    lpopt: int
    "Local potential cutoff radius"
    rcloc: float
    "Local potential delta value"
    dvloc0: float

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, LocalPotentialConfiguration):
            raise NotImplementedError(
                f"Cannot compare LocalPotentialConfiguration with {type(other)}"
            )
        if self.lloc != other.lloc:
            logging.debug("lloc differs: %d != %d", self.lloc, other.lloc)
            return False
        if self.lpopt != other.lpopt:
            logging.debug("lpopt differs: %d != %d", self.lpopt, other.lpopt)
            return False
        if abs(self.rcloc - other.rcloc) >= ATOL:
            logging.debug("rcloc differs: %f != %f", self.rcloc, other.rcloc)
            return False
        if abs(self.dvloc0 - other.dvloc0) >= ATOL:
            logging.debug("dvloc0 differs: %f != %f", self.dvloc0, other.dvloc0)
            return False
        return True

    @classmethod
    def from_str(cls, s: str) -> "LocalPotentialConfiguration":
        """Create a LocalPotentialConfiguration from a string.

        Args:
            s: String in the format "lloc lpopt rcloc dvloc0"

        Returns:
            LocalPotentialConfiguration instance.
        """
        parts = s.split()
        if len(parts) != 4:
            raise ValueError(f"Invalid local potential configuration: {s}")
        lloc, lpopt, rcloc, dvloc0 = (
            int(parts[0]),
            int(parts[1]),
            fort_float(parts[2]),
            fort_float(parts[3]),
        )
        return cls(lloc=lloc, lpopt=lpopt, rcloc=rcloc, dvloc0=dvloc0)


@dataclass
class KBProjectorConfiguration:
    """Configuration of the Kleinman-Bylander projectors."""

    "Angular momentum"
    l: int
    "Number of projectors"
    nproj: int
    "Energy shift for unbound states"
    debl: float
    'n of the first well state with well_tail = "confined" (0 for the default)'
    well_n: int = 0

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, KBProjectorConfiguration):
            raise NotImplementedError(
                f"Cannot compare KBProjectorConfiguration with {type(other)}"
            )
        if self.l != other.l:
            logging.debug("l differs: %d != %d", self.l, other.l)
            return False
        if self.nproj != other.nproj:
            logging.debug("nproj differs: %d != %d", self.nproj, other.nproj)
            return False
        if abs(self.debl - other.debl) >= ATOL:
            logging.debug("debl differs: %f != %f", self.debl, other.debl)
            return False
        if self.well_n != other.well_n:
            logging.debug("well_n differs: %d != %d", self.well_n, other.well_n)
            return False
        return True

    @classmethod
    def from_str(cls, s: str) -> "KBProjectorConfiguration":
        """Create a KBProjectorConfiguration from a string.

        Args:
            s: String in the format "l nproj debl [well_n]"

        Returns:
            KBProjectorConfiguration instance.
        """
        parts = s.split()
        if len(parts) not in (3, 4):
            raise ValueError(f"Invalid KB projector configuration: {s}")
        l, nproj, debl = int(parts[0]), int(parts[1]), fort_float(parts[2])
        well_n = int(parts[3]) if len(parts) == 4 else 0
        return cls(l=l, nproj=nproj, debl=debl, well_n=well_n)


@dataclass
class ModelCoreChargeConfiguration:
    """Configuration of the model core charge."""

    "Model core charge type"
    icmod: int
    "Scaling factor for the model core charge"
    fcfact: float | None = None
    "Minimum scaling factor for the model core charge (for optimization)"
    fcfact_min: float | None = None
    "Maximum scaling factor for the model core charge (for optimization)"
    fcfact_max: float | None = None
    "Step size for the scaling factor for the model core charge (for optimization)"
    fcfact_step: float | None = None
    "Cutoff radius scaling factor"
    rcfact: float | None = None
    "Minimum cutoff radius scaling factor (for optimization)"
    rcfact_min: float | None = None
    "Maximum cutoff radius scaling factor (for optimization)"
    rcfact_max: float | None = None
    "Step size for the cutoff radius scaling factor (for optimization)"
    rcfact_step: float | None = None
    "Share of the cutoff the model core must fit (icmod 6, 7)"
    slack: float | None = None
    'Cutoff of the pseudo-waves, "max" or "weighted" (icmod 6, 7)'
    budget_cutoff: str | None = None

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, ModelCoreChargeConfiguration):
            raise NotImplementedError(
                f"Cannot compare ModelCoreChargeConfiguration with {type(other)}"
            )
        if self.icmod != other.icmod:
            logging.debug("icmod differs: %d != %d", self.icmod, other.icmod)
            return False
        if (self.fcfact is None) != (other.fcfact is None):
            logging.debug(
                "fcfact presence differs: %s != %s", self.fcfact, other.fcfact
            )
            return False
        elif self.fcfact is not None and other.fcfact is not None:
            if abs(self.fcfact - other.fcfact) >= ATOL:
                logging.debug("fcfact differs: %f != %f", self.fcfact, other.fcfact)
                return False
        if (self.fcfact_min is None) != (other.fcfact_min is None):
            logging.debug(
                "fcfact_min presence differs: %s != %s",
                self.fcfact_min,
                other.fcfact_min,
            )
            return False
        elif self.fcfact_min is not None and other.fcfact_min is not None:
            if abs(self.fcfact_min - other.fcfact_min) >= ATOL:
                logging.debug(
                    "fcfact_min differs: %f != %f", self.fcfact_min, other.fcfact_min
                )
                return False
        if (self.fcfact_max is None) != (other.fcfact_max is None):
            logging.debug(
                "fcfact_max presence differs: %s != %s",
                self.fcfact_max,
                other.fcfact_max,
            )
            return False
        elif self.fcfact_max is not None and other.fcfact_max is not None:
            if abs(self.fcfact_max - other.fcfact_max) >= ATOL:
                logging.debug(
                    "fcfact_max differs: %f != %f", self.fcfact_max, other.fcfact_max
                )
                return False
        if (self.fcfact_step is None) != (other.fcfact_step is None):
            logging.debug(
                "fcfact_step presence differs: %s != %s",
                self.fcfact_step,
                other.fcfact_step,
            )
            return False
        elif self.fcfact_step is not None and other.fcfact_step is not None:
            if abs(self.fcfact_step - other.fcfact_step) >= ATOL:
                logging.debug(
                    "fcfact_step differs: %f != %f", self.fcfact_step, other.fcfact_step
                )
                return False
        if (self.rcfact is None) != (other.rcfact is None):
            logging.debug(
                "rcfact presence differs: %s != %s", self.rcfact, other.rcfact
            )
            return False
        elif self.rcfact is not None and other.rcfact is not None:
            if abs(self.rcfact - other.rcfact) >= ATOL:
                logging.debug("rcfact differs: %f != %f", self.rcfact, other.rcfact)
                return False
        if (self.rcfact_min is None) != (other.rcfact_min is None):
            logging.debug(
                "rcfact_min presence differs: %s != %s",
                self.rcfact_min,
                other.rcfact_min,
            )
            return False
        elif self.rcfact_min is not None and other.rcfact_min is not None:
            if abs(self.rcfact_min - other.rcfact_min) >= ATOL:
                logging.debug(
                    "rcfact_min differs: %f != %f", self.rcfact_min, other.rcfact_min
                )
                return False
        if (self.rcfact_max is None) != (other.rcfact_max is None):
            logging.debug(
                "rcfact_max presence differs: %s != %s",
                self.rcfact_max,
                other.rcfact_max,
            )
            return False
        elif self.rcfact_max is not None and other.rcfact_max is not None:
            if abs(self.rcfact_max - other.rcfact_max) >= ATOL:
                logging.debug(
                    "rcfact_max differs: %f != %f", self.rcfact_max, other.rcfact_max
                )
                return False
        if (self.rcfact_step is None) != (other.rcfact_step is None):
            logging.debug(
                "rcfact_step presence differs: %s != %s",
                self.rcfact_step,
                other.rcfact_step,
            )
            return False
        elif self.rcfact_step is not None and other.rcfact_step is not None:
            if abs(self.rcfact_step - other.rcfact_step) >= ATOL:
                logging.debug(
                    "rcfact_step differs: %f != %f", self.rcfact_step, other.rcfact_step
                )
                return False
        if (self.slack is None) != (other.slack is None):
            logging.debug("slack presence differs: %s != %s", self.slack, other.slack)
            return False
        elif self.slack is not None and other.slack is not None:
            if abs(self.slack - other.slack) >= ATOL:
                logging.debug("slack differs: %f != %f", self.slack, other.slack)
                return False
        if (self.budget_cutoff or "max") != (other.budget_cutoff or "max"):
            logging.debug(
                "budget_cutoff differs: %s != %s",
                self.budget_cutoff,
                other.budget_cutoff,
            )
            return False
        return True

    @classmethod
    def from_str(cls, s: str) -> "ModelCoreChargeConfiguration":
        """Create a ModelCoreChargeConfiguration from a string.

        Args:
            s: String in the format
               "icmod [fcfact [rcfact [fcfact_min fcfact_max fcfact_step
               rcfact_min rcfact_max rcfact_step]]]" (the grid with icmod 4, 5)

        Returns:
            ModelCoreChargeConfiguration instance.
        """
        parts = s.split()
        icmod = int(parts[0])
        if icmod == 0:
            return cls(icmod=icmod)
        elif icmod == 1:
            return cls(icmod=icmod, fcfact=fort_float(parts[1]))
        elif icmod == 2:
            return cls(
                icmod=icmod,
                fcfact=fort_float(parts[1]),
            )
        elif icmod == 3:
            return cls(
                icmod=icmod, fcfact=fort_float(parts[1]), rcfact=fort_float(parts[2])
            )
        elif icmod in (4, 5):
            if len(parts) < 9:
                return cls(icmod=icmod)
            elif len(parts) == 9:
                return cls(
                    icmod=icmod,
                    fcfact_min=fort_float(parts[3]),
                    fcfact_max=fort_float(parts[4]),
                    fcfact_step=fort_float(parts[5]),
                    rcfact_min=fort_float(parts[6]),
                    rcfact_max=fort_float(parts[7]),
                    rcfact_step=fort_float(parts[8]),
                )
            else:
                raise ValueError(f"Invalid model core charge configuration: {s}")
        elif icmod in (6, 7):
            # icmod, fcfact (cutoff budget, Ha), rcfact (xc energy error, Ha; default
            # 1e-5), slack (share of the budget the core must fit; default 1),
            # budget_cutoff (0 "max", 1 "weighted"; default 0)
            return cls(
                icmod=icmod,
                fcfact=fort_float(parts[1]),
                rcfact=fort_float(parts[2]) if len(parts) > 2 else None,
                slack=fort_float(parts[3]) if len(parts) > 3 else None,
                budget_cutoff=(
                    ("weighted" if int(parts[4]) == 1 else "max")
                    if len(parts) > 4
                    else None
                ),
            )
        else:
            raise ValueError(f"Invalid model core charge configuration: {s}")


@dataclass
class LogDerivativeConfiguration:
    """Configuration of the log-derivative grid."""

    "Inner log-derivative boundary"
    epsh1: float
    "Outer log-derivative boundary"
    epsh2: float
    "Log-derivative step size"
    depsh: float
    "Exchange-correlation radius (if any)"
    rxpsh: float | None = None

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, LogDerivativeConfiguration):
            raise NotImplementedError(
                f"Cannot compare LogDerivativeConfiguration with {type(other)}"
            )
        if (self.rxpsh is None) != (other.rxpsh is None):
            logging.debug("rxpsh presence differs: %s != %s", self.rxpsh, other.rxpsh)
            return False
        if self.rxpsh is not None and other.rxpsh is not None:
            if abs(self.rxpsh - other.rxpsh) >= ATOL:
                logging.debug("rxpsh differs: %f != %f", self.rxpsh, other.rxpsh)
                return False
        if abs(self.epsh1 - other.epsh1) >= ATOL:
            logging.debug("epsh1 differs: %f != %f", self.epsh1, other.epsh1)
            return False
        if abs(self.epsh2 - other.epsh2) >= ATOL:
            logging.debug("epsh2 differs: %f != %f", self.epsh2, other.epsh2)
            return False
        if abs(self.depsh - other.depsh) >= ATOL:
            logging.debug("depsh differs: %f != %f", self.depsh, other.depsh)
            return False
        return True

    @classmethod
    def from_str(cls, s: str) -> "LogDerivativeConfiguration":
        """Create a LogDerivativeConfiguration from a string.

        Args:
            s: String in the format "epsh1 epsh2 depsh [rxpsh]"

        Returns:
            LogDerivativeConfiguration instance.
        """
        parts = s.split()
        if len(parts) not in (3, 4):
            raise ValueError(f"Invalid log-derivative configuration: {s}")
        epsh1, epsh2, depsh = (
            fort_float(parts[0]),
            fort_float(parts[1]),
            fort_float(parts[2]),
        )
        rxpsh = fort_float(parts[3]) if len(parts) == 4 else None
        return cls(epsh1=epsh1, epsh2=epsh2, depsh=depsh, rxpsh=rxpsh)


@dataclass
class LinearMeshConfiguration:
    """Configuration of the linear mesh for output."""

    "Linear radial mesh step size"
    drl: float
    "Maximum radius for the linear mesh"
    rlmax: float
    "Most points of the log mesh of the UPF file (0 for the linear mesh), [pp_output] in TOML"
    upf_mesh: int = 0

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, LinearMeshConfiguration):
            raise NotImplementedError(
                f"Cannot compare LinearMeshConfiguration with {type(other)}"
            )
        if abs(self.drl - other.drl) >= ATOL:
            logging.debug("drl differs: %f != %f", self.drl, other.drl)
            return False
        if abs(self.rlmax - other.rlmax) >= ATOL:
            logging.debug("rlmax differs: %f != %f", self.rlmax, other.rlmax)
            return False
        if self.upf_mesh != other.upf_mesh:
            logging.debug("upf_mesh differs: %d != %d", self.upf_mesh, other.upf_mesh)
            return False
        return True

    @classmethod
    def from_str(cls, s: str) -> "LinearMeshConfiguration":
        """Create a LinearMeshConfiguration from a string.

        Args:
            s: String in the format "rlmax drl [upf_mesh]"

        Returns:
            LinearMeshConfiguration instance.
        """
        parts = s.split()
        if len(parts) not in (2, 3):
            raise ValueError(f"Invalid linear mesh configuration: {s}")
        rlmax, drl = fort_float(parts[0]), fort_float(parts[1])
        upf_mesh = int(parts[2]) if len(parts) == 3 else 0
        return cls(drl=drl, rlmax=rlmax, upf_mesh=upf_mesh)


@dataclass
class LogMeshConfiguration:
    """Configuration of the log radial mesh, the optional last line of the text input. A value of
    0 or less selects ONCVPSP's built-in value."""

    "Ratio of successive mesh points"
    amesh: float
    "First mesh point"
    rr1: float
    "Outer radius of the mesh"
    rmax: float

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, LogMeshConfiguration):
            raise NotImplementedError(
                f"Cannot compare LogMeshConfiguration with {type(other)}"
            )
        for name in ("amesh", "rr1", "rmax"):
            if abs(getattr(self, name) - getattr(other, name)) >= ATOL:
                logging.debug(
                    "%s differs: %f != %f",
                    name,
                    getattr(self, name),
                    getattr(other, name),
                )
                return False
        return True

    @classmethod
    def from_str(cls, s: str) -> "LogMeshConfiguration":
        """Create a LogMeshConfiguration from a string.

        Args:
            s: String in the format "amesh rr1 rmax"

        Returns:
            LogMeshConfiguration instance.
        """
        parts = s.split()
        if len(parts) != 3:
            raise ValueError(f"Invalid log mesh configuration: {s}")
        return cls(*(fort_float(p) for p in parts))


@dataclass
class OncvpspInput:
    """Input parameters for ONCVPSP."""

    atsym: str
    zz: float
    nc: int
    nv: int
    iexc: int
    psfile: str
    reference_configuration: list[StateConfiguration]
    lmax: int
    pseudopotentials: list[PseudopotentialConfiguration]
    local_potential: LocalPotentialConfiguration
    vkb_projectors: list[KBProjectorConfiguration]
    model_core_charge: ModelCoreChargeConfiguration
    log_derivative_analysis: LogDerivativeConfiguration
    linear_mesh: LinearMeshConfiguration
    test_configurations: list[list[StateConfiguration]]
    log_mesh: LogMeshConfiguration | None = None
    c_light: float | None = None
    "The optional values of the lmax line, at their defaults when absent"
    ep_as_given: bool = False
    well_onset: float | None = None
    well_tail: str = "well"
    wave_regularization: float = 0.0
    'None, or "space" or "rc": the normalization of the well states'
    well_norm: str | None = None

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, OncvpspInput):
            raise NotImplementedError(f"Cannot compare OncvpspInput with {type(other)}")
        if self.atsym != other.atsym:
            logging.debug("atsym differs: %s != %s", self.atsym, other.atsym)
            return False
        if abs(self.zz - other.zz) >= ATOL:
            logging.debug("zz differs: %f != %f", self.zz, other.zz)
            return False
        if self.nc != other.nc:
            logging.debug("nc differs: %d != %d", self.nc, other.nc)
            return False
        if self.nv != other.nv:
            logging.debug("nv differs: %d != %d", self.nv, other.nv)
            return False
        if self.iexc != other.iexc:
            logging.debug("iexc differs: %d != %d", self.iexc, other.iexc)
            return False
        if self.psfile != other.psfile:
            logging.debug("psfile differs: %s != %s", self.psfile, other.psfile)
            return False
        if self.c_light != other.c_light:
            logging.debug("c_light differs: %s != %s", self.c_light, other.c_light)
            return False
        if len(self.reference_configuration) != len(other.reference_configuration):
            logging.debug(
                "reference_configuration length differs: %d != %d",
                len(self.reference_configuration),
                len(other.reference_configuration),
            )
            return False
        for i, (s1, s2) in enumerate(
            zip(self.reference_configuration, other.reference_configuration)
        ):
            if s1 != s2:
                logging.debug(
                    "reference_configuration entry at index %d differs: %s != %s",
                    i,
                    s1,
                    s2,
                )
                return False
        if self.lmax != other.lmax:
            logging.debug("lmax differs: %d != %d", self.lmax, other.lmax)
            return False
        for name in (
            "ep_as_given",
            "well_onset",
            "well_tail",
            "wave_regularization",
            "well_norm",
        ):
            if getattr(self, name) != getattr(other, name):
                logging.debug(
                    "%s differs: %s != %s",
                    name,
                    getattr(self, name),
                    getattr(other, name),
                )
                return False
        if len(self.pseudopotentials) != len(other.pseudopotentials):
            logging.debug(
                "pseudopotentials length differs: %d != %d",
                len(self.pseudopotentials),
                len(other.pseudopotentials),
            )
            return False
        for pp1, pp2 in zip(self.pseudopotentials, other.pseudopotentials):
            if pp1 != pp2:
                logging.debug("pseudopotentials entry differs: %s != %s", pp1, pp2)
                return False
        if self.local_potential != other.local_potential:
            logging.debug(
                "local_potential differs: %s != %s",
                self.local_potential,
                other.local_potential,
            )
            return False
        if len(self.vkb_projectors) != len(other.vkb_projectors):
            logging.debug(
                "vkb_projectors length differs: %d != %d",
                len(self.vkb_projectors),
                len(other.vkb_projectors),
            )
            return False
        for proj1, proj2 in zip(self.vkb_projectors, other.vkb_projectors):
            if proj1 != proj2:
                logging.debug("vkb_projectors entry differs: %s != %s", proj1, proj2)
                return False
        if self.model_core_charge != other.model_core_charge:
            logging.debug(
                "model_core_charge differs: %s != %s",
                self.model_core_charge,
                other.model_core_charge,
            )
            return False
        if self.log_derivative_analysis != other.log_derivative_analysis:
            logging.debug(
                "log_derivative_analysis differs: %s != %s",
                self.log_derivative_analysis,
                other.log_derivative_analysis,
            )
            return False
        if self.linear_mesh != other.linear_mesh:
            logging.debug(
                "linear_mesh differs: %s != %s",
                self.linear_mesh,
                other.linear_mesh,
            )
            return False
        if len(self.test_configurations) != len(other.test_configurations):
            logging.debug(
                "test_configurations length differs: %d != %d",
                len(self.test_configurations),
                len(other.test_configurations),
            )
            return False
        for cfg1, cfg2 in zip(self.test_configurations, other.test_configurations):
            if len(cfg1) != len(cfg2):
                logging.debug(
                    "test_configurations entry length differs: %d != %d",
                    len(cfg1),
                    len(cfg2),
                )
                return False
            for s1, s2 in zip(cfg1, cfg2):
                if s1 != s2:
                    logging.debug("test_configurations entry differs: %s != %s", s1, s2)
                    return False
        if self.log_mesh != other.log_mesh:
            logging.debug("log_mesh differs: %s != %s", self.log_mesh, other.log_mesh)
            return False
        return True

    @classmethod
    def from_text_file(cls, file: ty.Any) -> "OncvpspInput":
        """Create an OncvpspInput from a file-like object.

        Args:
            file: File-like object to read from.

        Returns:
            OncvpspInput instance.
        """
        if isinstance(file, (str, pl.Path)):
            with open(file, "r", encoding="utf-8") as f:
                return cls.from_lines(f.readlines())
        else:
            return cls.from_lines(file.readlines())

    @classmethod
    def from_str(cls, s: str) -> "OncvpspInput":
        """Create an OncvpspInput from a string.

        Args:
            s: String to read from.

        Returns:
            OncvpspInput instance.
        """
        return cls.from_lines(s.splitlines())

    @classmethod
    def from_lines(cls, lines: list[str]) -> "OncvpspInput":
        """Create an OncvpspInput from a list of lines.

        Args:
            lines: List of lines to read from.
        Returns:
            OncvpspInput instance.
        """

        def non_comment_lines(lines: list[str]) -> ty.Generator[str, None, None]:
            for line in lines:
                line = line.strip()
                # Skip empty lines and full-line comments
                if line and not line.startswith("#"):
                    # Remove inline comments
                    line = line.split("#", 1)[0].strip()
                    yield line

        lines = non_comment_lines(lines)

        first_line = next(lines)
        words = first_line.split()
        atsym = words[0]
        zz = fort_float(words[1])
        nc = int(words[2])
        nv = int(words[3])
        iexc = int(words[4])
        psfile = words[5]
        c_light = fort_float(words[6]) if len(words) > 6 else None

        reference_configuration = [
            StateConfiguration.from_str(next(lines)) for _ in range(nc + nv)
        ]

        # lmax, and optionally ep_as_given (0 or 1), well_onset (0 or less for
        # the default), well_tail (0 for "well", 1 for "free", 2 for
        # "confined"), wave_regularization and well_norm (0 for that of
        # well_tail, 1 for "space", 2 for "rc")
        words = next(lines).split()
        lmax = int(words[0])
        ep_as_given = len(words) > 1 and int(words[1]) == 1
        well_onset = fort_float(words[2]) if len(words) > 2 else 0.0
        well_tail = (
            ("well", "free", "confined")[int(words[3])] if len(words) > 3 else "well"
        )
        wave_regularization = fort_float(words[4]) if len(words) > 4 else 0.0
        well_norm = (None, "space", "rc")[int(words[5])] if len(words) > 5 else None

        pseudopotentials = [
            PseudopotentialConfiguration.from_str(next(lines)) for _ in range(lmax + 1)
        ]

        local_potential = LocalPotentialConfiguration.from_str(next(lines))

        vkb_projectors = [
            KBProjectorConfiguration.from_str(next(lines)) for _ in range(lmax + 1)
        ]

        model_core_charge = ModelCoreChargeConfiguration.from_str(next(lines))

        log_derivative_analysis = LogDerivativeConfiguration.from_str(next(lines))

        linear_mesh = LinearMeshConfiguration.from_str(next(lines))

        ncnf = int(next(lines))
        test_configurations = []
        for _ in range(ncnf):
            nvcnf = int(next(lines))
            test_configurations.append(
                [StateConfiguration.from_str(next(lines)) for _ in range(nvcnf)]
            )

        log_mesh_line = next(lines, None)
        log_mesh = (
            LogMeshConfiguration.from_str(log_mesh_line) if log_mesh_line else None
        )

        return cls(
            atsym,
            zz,
            nc,
            nv,
            iexc,
            psfile,
            reference_configuration,
            lmax,
            pseudopotentials,
            local_potential,
            vkb_projectors,
            model_core_charge,
            log_derivative_analysis,
            linear_mesh,
            test_configurations,
            log_mesh,
            c_light,
            ep_as_given,
            well_onset if well_onset > 0 else None,
            well_tail,
            wave_regularization,
            well_norm,
        )

    @classmethod
    def from_dict(cls, d: dict[str, ty.Any]) -> "OncvpspInput":
        """Create an OncvpspInput from a dictionary.

        Args:
            d: Dictionary to read from.

        Returns:
            OncvpspInput instance.
        """
        return cls(
            atsym=d["oncvpsp"]["atsym"],
            zz=d["oncvpsp"]["z"],
            nc=d["oncvpsp"]["nc"],
            nv=d["oncvpsp"]["nv"],
            iexc=d["oncvpsp"]["iexc"],
            psfile=d["pp_output"]["psfile"],
            reference_configuration=[
                StateConfiguration(n=s[0], l=s[1], f=s[2])
                for s in zip(
                    d["reference_configuration"]["n"],
                    d["reference_configuration"]["l"],
                    d["reference_configuration"]["f"],
                )
            ],
            lmax=d["pseudopotentials"]["lmax"],
            pseudopotentials=[
                PseudopotentialConfiguration(pp[0], pp[1], pp[2], pp[3], pp[4], pp[5])
                for pp in zip(
                    d["pseudopotentials"]["l"],
                    d["pseudopotentials"]["rc"],
                    d["pseudopotentials"]["ep"],
                    d["pseudopotentials"]["ncon"],
                    d["pseudopotentials"]["nbas"],
                    d["pseudopotentials"]["qcut"],
                )
            ],
            local_potential=LocalPotentialConfiguration(
                lloc=d["local_potential"]["lloc"],
                lpopt=d["local_potential"]["lpopt"],
                rcloc=d["local_potential"]["rcloc"],
                dvloc0=d["local_potential"]["dvloc0"],
            ),
            vkb_projectors=[
                KBProjectorConfiguration(proj[0], proj[1], proj[2], proj[3])
                for proj in zip(
                    d["vkb_projectors"]["l"],
                    d["vkb_projectors"]["nproj"],
                    d["vkb_projectors"]["debl"],
                    d["vkb_projectors"].get(
                        "well_n", [0] * len(d["vkb_projectors"]["l"])
                    ),
                )
            ],
            model_core_charge=ModelCoreChargeConfiguration(
                icmod=d["model_core_charge"]["icmod"],
                fcfact=d["model_core_charge"].get("fcfact"),
                fcfact_min=d["model_core_charge"].get("fcfact_min"),
                fcfact_max=d["model_core_charge"].get("fcfact_max"),
                fcfact_step=d["model_core_charge"].get("fcfact_step"),
                rcfact=d["model_core_charge"].get("rcfact"),
                rcfact_min=d["model_core_charge"].get("rcfact_min"),
                rcfact_max=d["model_core_charge"].get("rcfact_max"),
                rcfact_step=d["model_core_charge"].get("rcfact_step"),
                slack=d["model_core_charge"].get("slack"),
                budget_cutoff=d["model_core_charge"].get("budget_cutoff"),
            ),
            log_derivative_analysis=LogDerivativeConfiguration(
                epsh1=d["log_derivative_analysis"]["epsh1"],
                epsh2=d["log_derivative_analysis"]["epsh2"],
                depsh=d["log_derivative_analysis"]["depsh"],
                rxpsh=d["log_derivative_analysis"].get("rxpsh"),
            ),
            linear_mesh=LinearMeshConfiguration(
                drl=d["linear_mesh"]["a"],
                rlmax=d["linear_mesh"]["rmax"],
                upf_mesh=d["pp_output"].get("upf_mesh", 0),
            ),
            test_configurations=[
                [
                    StateConfiguration(n=s[0], l=s[1], f=s[2])
                    for s in zip(cfg["n"], cfg["l"], cfg["f"])
                ]
                for cfg in d["test_configurations"]
            ],
            log_mesh=(
                LogMeshConfiguration(
                    amesh=d["log_mesh"].get("amesh", -1.0),
                    rr1=d["log_mesh"].get("rr1", -1.0),
                    rmax=d["log_mesh"].get("rmax", -1.0),
                )
                if "log_mesh" in d
                else None
            ),
            c_light=d["oncvpsp"].get("c_light"),
            ep_as_given=d["pseudopotentials"].get("ep_as_given", False),
            well_onset=d["pseudopotentials"].get("well_onset"),
            well_tail=d["pseudopotentials"].get("well_tail", "well"),
            wave_regularization=d["pseudopotentials"].get("wave_regularization", 0.0),
            well_norm=d["pseudopotentials"].get("well_norm"),
        )

    def to_dict(self) -> dict[str, ty.Any]:
        """Convert an OncvpspInput to a dictionary in the style of the TOML input.

        Returns:
            Dictionary representation of the OncvpspInput.
        """
        d = {
            "oncvpsp": {
                "atsym": self.atsym,
                "z": self.zz,
                "nc": self.nc,
                "nv": self.nv,
                "iexc": self.iexc,
                **({} if self.c_light is None else {"c_light": self.c_light}),
            },
            "linear_mesh": {
                "a": self.linear_mesh.drl,
                "rmax": self.linear_mesh.rlmax,
            },
            "reference_configuration": {
                "n": [s.n for s in self.reference_configuration],
                "l": [s.l for s in self.reference_configuration],
                "f": [s.f for s in self.reference_configuration],
            },
            "pseudopotentials": {
                "lmax": self.lmax,
                "l": [pp.l for pp in self.pseudopotentials],
                "rc": [pp.rc for pp in self.pseudopotentials],
                "ep": [pp.ep for pp in self.pseudopotentials],
                "ncon": [pp.ncon for pp in self.pseudopotentials],
                "nbas": [pp.nbas for pp in self.pseudopotentials],
                "qcut": [pp.qcut for pp in self.pseudopotentials],
            },
            "local_potential": {
                "lloc": self.local_potential.lloc,
                "lpopt": self.local_potential.lpopt,
                "rcloc": self.local_potential.rcloc,
                "dvloc0": self.local_potential.dvloc0,
            },
            "vkb_projectors": {
                "l": [proj.l for proj in self.vkb_projectors],
                "nproj": [proj.nproj for proj in self.vkb_projectors],
                "debl": [proj.debl for proj in self.vkb_projectors],
            },
            "model_core_charge": {
                "icmod": self.model_core_charge.icmod,
                "fcfact": self.model_core_charge.fcfact,
                "fcfact_min": self.model_core_charge.fcfact_min,
                "fcfact_max": self.model_core_charge.fcfact_max,
                "fcfact_step": self.model_core_charge.fcfact_step,
                "rcfact": self.model_core_charge.rcfact,
                "rcfact_min": self.model_core_charge.rcfact_min,
                "rcfact_max": self.model_core_charge.rcfact_max,
                "rcfact_step": self.model_core_charge.rcfact_step,
                "slack": self.model_core_charge.slack,
                "budget_cutoff": self.model_core_charge.budget_cutoff,
            },
            "log_derivative_analysis": {
                "epsh1": self.log_derivative_analysis.epsh1,
                "epsh2": self.log_derivative_analysis.epsh2,
                "depsh": self.log_derivative_analysis.depsh,
                "rxpsh": self.log_derivative_analysis.rxpsh,
            },
            "pp_output": {
                "psfile": self.psfile,
            },
            "test_configurations": [
                {
                    "n": [s.n for s in cfg],
                    "l": [s.l for s in cfg],
                    "f": [s.f for s in cfg],
                }
                for cfg in self.test_configurations
            ],
        }
        if self.ep_as_given:
            d["pseudopotentials"]["ep_as_given"] = True
        if self.well_onset is not None:
            d["pseudopotentials"]["well_onset"] = self.well_onset
        if self.well_tail != "well":
            d["pseudopotentials"]["well_tail"] = self.well_tail
        if self.wave_regularization != 0.0:
            d["pseudopotentials"]["wave_regularization"] = self.wave_regularization
        if self.well_norm is not None:
            d["pseudopotentials"]["well_norm"] = self.well_norm
        if any(proj.well_n for proj in self.vkb_projectors):
            d["vkb_projectors"]["well_n"] = [
                proj.well_n for proj in self.vkb_projectors
            ]
        if self.linear_mesh.upf_mesh != 0:
            d["pp_output"]["upf_mesh"] = self.linear_mesh.upf_mesh
        if self.log_mesh is not None:
            d["log_mesh"] = {
                name: value
                for name in ("amesh", "rr1", "rmax")
                if (value := getattr(self.log_mesh, name)) > 0
            }
        return d


# %%
def main() -> None:
    parser = argparse.ArgumentParser(
        description="Convert ONCVPSP input file to TOML format."
    )
    parser.add_argument(
        "-i",
        "--input",
        type=str,
        required=True,
        help="Input ONCVPSP text file (default: stdin)",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        default=None,
        help="Output TOML file (default: stdout)",
    )
    parser.add_argument(
        "--check",
        action="store_true",
        help="Check conversion by performing a round-trip",
    )
    args = parser.parse_args()

    input_path = pl.Path(args.input)
    if not input_path.is_file():
        print(f"Error: Input file '{input_path}' does not exist.", file=sys.stderr)
        sys.exit(1)

    oncvpsp_input = OncvpspInput.from_text_file(input_path)
    toml_dict = oncvpsp_input.to_dict()
    toml_str = toml.dumps(toml_dict)

    if args.output:
        output_path = pl.Path(args.output)
        with open(output_path, "w", encoding="utf-8") as f:
            f.write(toml_str)
    else:
        print(toml_str)

    if args.check:
        toml_dict_loaded = toml.loads(toml_str)
        oncvpsp_input_loaded = OncvpspInput.from_dict(toml_dict_loaded)
        if oncvpsp_input != oncvpsp_input_loaded:
            print("Error: Round-trip conversion check failed.", file=sys.stderr)
            sys.exit(1)


if __name__ == "__main__":
    main()
