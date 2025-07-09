#!/usr/bin/env python

"""
Downsamples a PARCHG file. Might work on related files??
Used for making the files smaller.
Does two things: reduces the recorded precision of the numbers,
and reduces the FFT mesh (summing adjacent points).

Makes use of code snippets adapted from
- pymatgen (https://pymatgen.org)
- fortranformat (https://pypi.org/project/fortranformat/)

Copyright 2025 Bernard Field
Copyright 2012 Materials Project (part of Parchg.to_file)
Copyright 2011 Brendan Arnold (part of format_float)

MIT License
"""

import argparse
import math
import itertools

import numpy as np

class Parchg:
    def __init__(self, header:list[str], data:np.ndarray):
        self.header = header
        self.data = data

    @property
    def shape(self) -> tuple[int, int, int]:
        return self.data.shape
    
    @classmethod
    def from_file(cls, f):
        # If f is a filename, open the file.
        if isinstance(f, str):
            with open(f, 'r') as f2:
                return cls.from_file(f2)
        # Read the first few lines which specify the unit cell
        header1 = [f.readline() for _ in range(7)]
        # This 7th line is also the count of elements,
        # which tells us the length of the next part of the header:
        # the atomic structure.
        natoms = sum(int(x) for x in header1[6].split())
        header2 = [f.readline() for _ in range(natoms+2)]
        # +2 because there's also the Direct/Cartesian line,
        # and the blank line at the end.
        # Read out the shape of the FFT grid.
        dim = tuple(int(x) for x in f.readline().split())
        ngrid_pts = dim[0] * dim[1] * dim[2]
        # Combine the headers
        header = header1 + header2
        # Now, let us read the charge density data.
        data_count = 0
        dataset = np.zeros(dim)
        for line in f:
            # Use a code snippet from Pymatgen
            for tok in line.split():
                if data_count < ngrid_pts:
                    # This complicated procedure is necessary because
                    # VASP outputs x as the fastest index, followed by y
                    # then z.
                    no_x = data_count // dim[0]
                    try:
                        dataset[data_count % dim[0], no_x % dim[1], no_x // dim[1]] = float(tok)
                    except ValueError:
                        dataset[data_count % dim[0], no_x % dim[1], no_x // dim[1]] = np.nan
                    data_count += 1
            if data_count >= ngrid_pts:
                break
        # Return the data
        return cls(header, dataset)
    
    def to_file(self, f, digits:int = 5):
        """
        Writes the PARCHG file.

        May specify the precision of the numbers, `digits`.
        (It is not meaningful to increase digits beyond the precision in the source file,
        which is 5 for PARCHG.)
        """
        # If f is a filename, open the file.
        if isinstance(f, str):
            with open(f, 'w') as f2:
                return self.to_file(f2)
        # Write the header. It comes with newlines already
        f.writelines(self.header)
        # Write the dimensions
        # From VASP source code: format is IU,'(3I5)'
        dim = self.data.shape
        f.write('{:5d}{:5d}{:5d}\n'.format(*dim))
        # Write the grid.
        # From VASP source code, default format is 1X,G11.5 with 10 entries per row.
        # {:#11.5G}, except Python does it differently to Fortran.
        entries_per_row = 10
        # This one from pymatgen.io.vasp.outputs.VolumetricData.write_file
        lines = []
        count = 0
        for k, j, i in itertools.product(list(range(dim[2])), list(range(dim[1])), list(range(dim[0]))):
            lines.append(format_float(self.data[i, j, k], 6+digits, digits))
            count += 1
            if count % entries_per_row == 0:
                f.write(" " + "".join(lines) + "\n")
                lines = []
            else:
                lines.append(" ")
        if count % entries_per_row != 0:
            f.write(" " + "".join(lines) + " \n")
    
    def downsample(self, nx:int, ny:int, nz:int):
        """
        Downsamples the data, averaging neighbouring points.
        Reduces the FFT grid by a factor of nx, ny, nz in those respective directions.
        nx, ny, nz must be factors of the FFT grid size.

        N.B. The downsampling is currently a cheap downsample, in that
        the center of mass changes slightly, by a grid point. But, like, it'll do.
        """
        # Confirm factorisation.
        dim = self.data.shape
        if dim[0] / nx != dim[0] // nx:
            raise ValueError("x-grid ", dim[0], " is not divisible by ", nx)
        if dim[1] / ny != dim[1] // ny:
            raise ValueError("y-grid ", dim[1], " is not divisible by ", ny)
        if dim[2] / nz != dim[2] // nz:
            raise ValueError("z-grid ", dim[2], " is not divisible by ", nz)
        # Create new grid
        new_grid = np.zeros((dim[0]//nx, dim[1]//ny, dim[2]//nz))
        # Downsample
        # Probably take slice steps, then average all the different offsets of those slices.
        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    new_grid += self.data[i::nx, j::ny, k::nz] / nx / ny / nz
        # Record the data
        self.data = new_grid


def format_float(val:float, width:int, digits:int) -> str:
    """
    The exponential formatting part of this function is adapted from Pymatgen,
    but with extra logic to handle arbitrary widths and digits.
    The logic to check whether to do E or F mode is from https://pypi.org/project/fortranformat/

    A note: my VASP data adds an extra trailing zero to 10.000 and 100.00
    (to be 10.0000 and 100.000), at least on Gadi, but that is not the behaviour
    here or on NERSC's Fortran. This is a compiler-dependent thing.
    """
    # Catch the error case, where the original data was bad.
    if val is np.nan:
        return '*' * width
    # How much leading padding have we got? width - digits - 2 (ones place and decimal) - 4 (E+00)
    pad = width - digits - 2 - 4
    if pad < -1: # Not enough room!
        return '*' * width
    aval = abs(val)
    exp_d = 10 ** digits
    # Use exponential processing because we are sufficiently far from 1.
    # Is it less than 0.1 (to within rounding) or greater than 10**d (to within rounding)?
    if (0.0 <= aval < (0.1 - 0.05 / exp_d)) or (aval >= (exp_d - 0.5)):
        # Use Python's internal float processing. It, however, puts the most significant digit before
        # the decimal place, rather than after like in Fortran.
        fmt = f'{{:{width-1}.{digits-1}E}}'
        flt_str = fmt.format(val)
        # To make if Fortran-like, we just need to move the decimal place over, increment the exponent,
        # and add the leading zero.
        # Positive number
        if val >= 0:
            if pad >= 0:
                return pad * " " + f"0.{flt_str[pad]}{flt_str[pad+2:pad+digits+1]}E{int(flt_str[pad+digits+2:]) + 1:+03}"
            # pad == -1
            # Fortran will squeeze out an extra character by dropping the leading zero if it must.
            return f".{flt_str[0]}{flt_str[2:digits+1]}E{int(flt_str[digits+2:]) + 1:+03}"
        # Negative number.
        # Note that the minus sign takes up some of this padding.
        if pad == 0:
            # If there is no room, the minus sign takes the ones place.
            return f"-.{flt_str[1]}{flt_str[3:digits+2]}E{int(flt_str[digits+3:]) + 1:+03}"
        else:
            pad -= 1
            return pad * " " + f"-0.{flt_str[pad+1]}{flt_str[pad+3:pad+digits+2]}E{int(flt_str[pad+digits+3:]) + 1:+03}"
    else:
        # We are close to one, so do float processing.
        magnitude = math.floor(math.log10(aval) + math.log10(1+0.5/exp_d)) + 1
        # If close to the next magnitude (within rounding), make it the next magnitude.
        # (I note that not every Fortran compilation accounts for this rounding. But I shall.)
        # The +1 is to make it Fortran-like.
        # Fortran places the float such that its leading digit aligns with the E-format decimal.
        pad = width - digits - 1 - 4
        if val < 0:
            pad -= 1
        if magnitude <= 0:
            # The leading zero sits behind
            pad -= 1
        if digits - magnitude == 0:
            # Fortran adds the decimal place even if there is nothing after it
            fmt = ' ' * pad + '{:.0f}.' + ' ' * 4
        else:
            fmt = ' ' * pad + f'{{:.{digits - magnitude}f}}' + ' ' * 4
        flt_str = fmt.format(val)
        if pad == -1:
            # Too many digits. Strip the leading zero.
            if val >= 0:
                return flt_str[1:]
            else:
                return '-' + flt_str[2:]
        return flt_str


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('parchg', help="PARCHG file to read")
    parser.add_argument('output', help="Destination for downsampled parchg file.")
    parser.add_argument('-x', type=int, default=2, help="Downsample x-grid by this amount.")
    parser.add_argument('-y', type=int, default=2, help="Downsample y-grid by this amount.")
    parser.add_argument('-z', type=int, default=4, help="Downsample z-grid by this amount.")
    parser.add_argument('-d', '--digits', type=int, default=2,
                        help="Print this many significant digits. (PARCHG uses 5 normally.)")
    args = parser.parse_args()

    print("Reading PARCHG")
    chg = Parchg.from_file(args.parchg)
    chg.downsample(args.x, args.y, args.z)
    print("Writing new PARCHG")
    chg.to_file(args.output, digits=args.digits)