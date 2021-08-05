import numpy as np
import pathlib
import loadSP as lsp

"""
Find near-exact and close matches between pairs of RAMSES HOP
position (grpxxxx.pos) files. Return the list of indices for
matching halos based on location and DM particle counts.
"""


class FindSameHalos:
    def __init__(self, num1, prefix1, num2, prefix2, v=False):
        """
        Initialize the instance
        """
        self.file1 = pathlib.Path(prefix1 + "grp{:05d}.pos".format(num1))
        self.file2 = pathlib.Path(prefix2 + "grp{:05d}.pos".format(num2))
        self.verbose = v
        if self.verbose:
            print("Loading {} and {}".format(self.file1, self.file2))
        if not self.file1.is_file() or not self.file2.is_file():
            print("{} or {} not found.".format(self.file1, self.file2))
        return

    def load(self):
        """"""
        self.halos1 = np.loadtxt(self.file1, skiprows=1)
        self.halos1_pos = self.halos1[:, 2:5]  # position info
        self.halos1_part = self.halos1[:, 1]
        self.halos2 = np.loadtxt(self.file2, skiprows=1)
        self.halos2_pos = self.halos2[:, 2:5]  # position info
        self.halos2_part = self.halos2[:, 1]
        if self.verbose:
            print("Print first 5 entries...")
            for i in range(5):
                print(
                    "h1 {}= {} {:.3e} {:.3e} {:.3e}".format(
                        i,
                        self.halos1_part[i],
                        self.halos1_pos[i][0],
                        self.halos1_pos[i][1],
                        self.halos1_pos[i][2],
                    )
                )
                print(
                    "h2 {}= {} {:.3e} {:.3e} {:.3e}".format(
                        i,
                        self.halos2_part[i],
                        self.halos2_pos[i][0],
                        self.halos2_pos[i][1],
                        self.halos2_pos[i][2],
                    )
                )

        return

    def find1(self, loc):
        """
        Attempts to find the location nearest loc in the file file2
        The intent of this routine is to find the 'same' halo (@ loc)
        in another file... typically a different redshift
        """
        for i in range(self.halos2_part.size):
            dxyz = loc - self.halos2_pos[i]
            if self.verbose:
                print(
                    "h0   = {:.3e} {:.3e} {:.3e}".format(
                        loc[0],
                        loc[1],
                        loc[2],
                    )
                )
                print(
                    "h2 {}= {:.3e} {:.3e} {:.3e}".format(
                        i,
                        self.halos2_pos[i][0],
                        self.halos2_pos[i][1],
                        self.halos2_pos[i][2],
                    )
                )
                print("dr= ({:.3e} {:.3e} {:.3e})".format(dxyz[0], dxyz[1], dxyz[2]))
            dr = np.linalg.norm(dxyz)
            if self.verbose:
                print("dr for h0 and h{} is {:.5e}".format(i, dr))
                print(self.halos2_part[i] - np.array(loc))
            if dr <= 0.005:
                # We found an near-exact match. Call it good and move on.
                if self.verbose:
                    print("Match! 0 {} -- dr={}".format(i, dr))
                    print("dp = {}".format(np.abs(self.halos2_part[i] - loc)))
                return i
        # Didn't find it!
        return -1

    def find(self):
        """
        Tries to find a very close match between halo locations in two files
        Returns a list of indices of matching halos
        """
        #
        # Approach:
        # Loop through halo locations in halos1.
        # Look at index +-10 for matching halos in halos2.
        # A match is found when norm(x1,y1,z1 - x2,y2,z2) < 10%
        #
        pairings = []
        close = False
        for o in range(self.halos1_part.size):
            # Look for a match (in the halo2 list) in the interval of
            # size 100 around the current halo1 halo.
            for i in range(max(0, o - 50), min(self.halos2_part.size, o + 50)):
                if len(pairings) and np.isin(i, np.array(pairings)[:, 1]):
                    if self.verbose:
                        print("Halo indx already paired {}".format(i))
                    continue
                dxyz = self.halos1_pos[o] - self.halos2_pos[i]
                if self.verbose:
                    print(
                        "h1 {}= {:.3e} {:.3e} {:.3e}".format(
                            o,
                            self.halos1_pos[o][0],
                            self.halos1_pos[o][1],
                            self.halos1_pos[o][2],
                        )
                    )
                    print(
                        "h2 {}= {:.3e} {:.3e} {:.3e}".format(
                            i,
                            self.halos2_pos[i][0],
                            self.halos2_pos[i][1],
                            self.halos2_pos[i][2],
                        )
                    )
                    print(
                        "dr= ({:.3e} {:.3e} {:.3e})".format(dxyz[0], dxyz[1], dxyz[2])
                    )
                dr = np.linalg.norm(dxyz)
                if self.verbose:
                    print("dr for h{} and h{} is {:.5e}".format(o, i, dr))
                if (
                    dr <= 0.005
                    and np.abs(self.halos2_part[i] - self.halos1_part[o])
                    < self.halos2_part[i] * 0.1
                ):
                    # We found an near-exact match. Call it good and move on.
                    if self.verbose:
                        print("Match! {} {} -- dr={}".format(o, i, dr))
                        print(
                            "dp = {}".format(
                                np.abs(self.halos2_part[i] - self.halos1_part[o])
                            )
                        )
                    pairings.append([o, i])
                    close = False  # Reset any close match found
                    break  # go on to next halo in halos1
                elif (
                    dr <= 0.01
                    and np.abs(self.halos2_part[i] - self.halos1_part[o])
                    < self.halos2_part[i] * 0.1
                ):
                    # We found a close match, but maybe we'll find something better
                    if not close:
                        # First close match found...
                        closeMatch = [o, i]
                    elif np.abs(i - o) < np.abs(closeMatch[0] - closeMatch[1]):
                        # Favor pairs closer together in the list
                        closeMatch = [o, i]
                        close = True
            # If we get here we haven't found an exact match. Maybe a close one?
            if close:
                if self.verbose:
                    print("Near Match! {} {}".format(closeMatch[0], closeMatch[1]))
                    print(
                        "dp = {}".format(
                            np.abs(
                                self.halos2_part[closeMatch[0]]
                                - self.halos1_part[closeMatch[1]]
                            )
                        )
                    )
                pairings.append(closeMatch)
                close = False
            if len(pairings) > 50:
                break
        return np.array(pairings)

    def getHalos1(self):
        return self.halos1_pos, self.halos1_part

    def getHalos2(self):
        return self.halos2_pos, self.halos2_part
