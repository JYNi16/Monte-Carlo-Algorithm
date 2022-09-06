import numpy as np

def cos(angle):
    """transforms angle from [0,1] to cos(2pi[0,1])"""
    return np.cos(2 * np.pi * angle)


def sin(angle):
    """transforms angle from [0,1] to sin(2pi[0,1])"""
    return np.sin(2 * np.pi * angle)


class XYMetropolis:

    def __init__(self,
                 lattice_shape,
                 beta=1,
                 J=5,
                 random_state=5,
                 initial_state='hot'):
        self.beta = beta
        self.J = J
        self.rs = np.random.RandomState(seed=random_state)

        # matrix of lattice angles
        if initial_state == 'hot':
            self.A = self.rs.rand(*lattice_shape)
        elif initial_state == 'cold':
            self.A = np.zeros(lattice_shape)
        else:
            raise ValueError('initial_state must be cold or hot')
        self.time = 0

        # Matrix of winding numbers
        self.V = np.zeros((lattice_shape[0] - 1, lattice_shape[1] - 1))

        # Correlations (since we have torus topology, we can start from the left top)
        self.corr_range = int(self.A.shape[0] / 2)

        # Correlation array of length corr_range
        self.C = []

        # Magnetization
        self.M = 0

        # Squared magnetization
        self.M2 = 0

        # Vortex density
        self.Vdensity = 0

    def step(self):
        """Perform one step of metropolis algorithm"""
        pos = tuple(self.rs.randint(_) for _ in self.A.shape)
        value = self.rs.rand()
        delta_H = self.dH(pos, value)
        if (delta_H < 0) or (self.rs.rand() < np.exp(-self.beta * delta_H)):
            self.A[pos] = value

    def dH(self, pos, val):
        """Calculate delta energy"""
        delta = 0
        old_val = self.A[pos]
        pos_list = list(pos)
        incr_delta = lambda pos: cos(self.A[pos] - val) - cos(self.A[pos] - old_val)
        for i in range(len(self.A.shape)):
            pos_list[i] += 1
            pos_list[i] %= self.A.shape[i]
            delta += incr_delta(tuple(pos_list))
            pos_list[i] -= 2
            pos_list[i] %= self.A.shape[i]
            delta += incr_delta(tuple(pos_list))
            pos_list[i] += 1
            pos_list[i] %= self.A.shape[i]
        return -delta * self.J

    def get_V(self):
        """Update matrix of winding numbers (Vortex matrix)"""
        for i in range(self.V.shape[0]):
            for j in range(self.V.shape[1]):

                # create list of angles from the below-down square
                a = [self.A[i, j],
                     self.A[i, j + 1],
                     self.A[i + 1, j + 1],
                     self.A[i + 1, j]]

                # run clockwise and calculate sum of angles
                a_sum = 0
                for k in range(len(a)):
                    d = a[k] - a[(k + 1) % len(a)]
                    if abs(d) > 0.5:
                        d -= np.sign(d)
                    a_sum += d

                self.V[i, j] = a_sum

        return self.V

    def get_Vdensity(self):
        self.Vdensity = np.sum(abs(self.V)) / 2 / np.prod(self.A.shape)
        return self.Vdensity

    def get_C(self):
        """Update correlations"""
        corrs_d = []  # correlations for each dim
        self.C = []

        # compute correlations over each shift (here means distance)
        for r in range(int(self.corr_range)):
            # and each axis
            for d in range(len(self.A.shape)):
                # calculate mean over all spins
                corr = np.mean(cos(self.A - np.roll(self.A, r, axis=d)))
                corrs_d.append(corr)
            self.C.append(np.mean(corrs_d))
        return self.C

    def get_M2(self):
        """Get squared magnetization"""
        self.M2 = ((np.sum(cos(self.A))) ** 2 + (np.sum(sin(self.A))) ** 2) / (np.prod(self.A.shape)) ** 2
        return self.M2

    def simulate(self, steps):
        for _ in range(steps):
            self.step()
        self.get_V()
        self.get_C()
        self.get_M2()
        self.get_Vdensity()


if __name__=="__main__":
    sim = XYMetropolis((50,50),
                 beta=1,
                 J=5,
                 random_state=5,
                 initial_state='hot')