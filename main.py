from dataclasses import dataclass
from typing import ClassVar, Literal, Self
import sys

from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import numpy as np
import matplotlib.pyplot as plt
from scipy.io import loadmat



@dataclass
class Data:

    MATLAB_DATA: ClassVar[str] = "data/eric.mat"
    RINGNMR_PATH_TMPL: ClassVar[str] = "data/ring_{}.txt"

    tau: float
    eric_nu1: np.ndarray
    eric_r1rho: np.ndarray
    ringnmr_nu1: np.ndarray
    ringnmr_r1rho: np.ndarray

    @classmethod
    def fetch(cls, temp: Literal[37, 57, 67, 77]) -> Self:
        matlab = loadmat(cls.MATLAB_DATA)
        tau = matlab["tau{}".format(temp)][0][0]
        eric_nu1 = matlab["w1"][0] / (2.0 * np.pi * 1.0e3)
        eric_r1rho = matlab["curve{}".format(temp)][0]
        ringnmr_nu1 = []
        ringnmr_r1rho = []
        with open(cls.RINGNMR_PATH_TMPL.format(temp), "r") as fh:
            txt = fh.read()
            for line in txt.split("\n"):
                if line == "":
                    continue
                nu1, r1rho = line.split(" ")
                ringnmr_nu1.append(float(nu1))
                ringnmr_r1rho.append(float(r1rho))
        ringnmr_nu1 = np.array(ringnmr_nu1)
        ringnmr_r1rho = np.pow(10.0, np.array(ringnmr_r1rho))
        return cls(tau, eric_nu1, eric_r1rho, ringnmr_nu1, ringnmr_r1rho)

    @property
    def tau_str(self) -> str:
        if 1.0e-9 <= self.tau < 1.0e-6:
            return "{:.1f} ns".format(1.0e9 * self.tau)
        elif 1.0e-6 <= self.tau < 1.0e-3:
            return "{:.1f} μs".format(1.0e6 * self.tau)
        elif 1.0e-3 <= self.tau < 1.0:
            return "{:.1f} ms".format(1.0e3 * self.tau)
        else:
            return "{:.1f} s".format(self.tau)


def get_data() -> dict[int, Data]:
    data = {}
    for temp in (37, 57, 67, 77):
        data[temp] = Data.fetch(temp)
    return data


def main():
    data = get_data()
    fig, ax = plt.subplots(
        gridspec_kw=dict(left=0.1, right=0.99, bottom=0.1, top=0.99)
    )
    legend_info: list[tuple[Line2D | Patch, str]] = [
        (Line2D([0], [0], color="k", alpha=0.5), "Keeler, McDermott"),
        (Line2D([0], [0], color="k", ls=":"), "RING NMR"),
    ]
    xmin, xmax = float("inf"), float("-inf")
    for temp, datum in data.items():
        xmin = min(min(xmin, datum.eric_nu1[0]), datum.ringnmr_nu1[0])
        xmax = max(max(xmax, datum.eric_nu1[-1]), datum.ringnmr_nu1[-1])
        color = ax.plot(
            datum.eric_nu1,
            datum.eric_r1rho,
            alpha=0.5,
        )[0].get_color()
        ax.plot(
            datum.ringnmr_nu1,
            datum.ringnmr_r1rho,
            color=color,
            ls=":",
        )
        legend_info.append(
            (
                Patch(color=color),
                "{}°C ($\\tau_{{\\text{{c}}}} =$ {})".format(temp, datum.tau_str)
            )
        )

    ax.set_xlim(xmin, xmax)
    ax.set_xlabel("$\\nu_1$ (kHz)")
    ax.set_ylabel("$R_{1\\rho}$ (s$^{-1}$)")
    ax.set_yscale("log")
    ax.legend([x[0] for x in legend_info], [x[1] for x in legend_info])
    fig.savefig("r1rho.pdf")


if __name__ == "__main__":
    main()
