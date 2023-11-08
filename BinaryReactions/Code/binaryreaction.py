"""
Contains the BinaryReaction class, which is used to store information about a
binary chemical reaction. It is used to simulate the reaction and plot the
results.
"""
import numpy as np
from scipy.integrate import solve_ivp
from typing import Any, Iterable

class BinaryReaction:
    def __init__(self,
                 init_conc:Iterable[float],
                 tau: Iterable[float],
                 k: Iterable[float],
                 s: Iterable[float],
                 ) -> None:
        """
        Initializes a BinaryReaction object with the given parameters.

        A + A <--> A + A* (equilibrium) with rate constants p and q
        A* --> B + C (non-equilibrium) with rate constant r

        Parameters
        ----------
        init_conc : Iterable[float]
            The initial concentrations of the species.
        tau : Iterable[float]
            Dimensionless time at which to evaluate the solution.
        k : Iterable[float]
            Dimensionless equilibrium constants.
        s : Iterable[float]
            Dimensionless rate constants.
        """
        self.init_conc = init_conc
        self.tau = tau
        # If k and s are not iterable, make them iterable.
        if not hasattr(k, '__iter__'):
            k = [k]
        if not hasattr(s, '__iter__'):
            s = [s]
        self.k = k
        self.s = s

        self.solve_for = (0, 0)

    def rate_eq(self, t:Iterable[float], y:Iterable[float]):
        """
        Define the rate equations for the BinaryReaction object.
        These are in a dimensionless form.
        """
        # Unpack the state vector
        a, a_star, b, c = y
        
        # Get rate constants
        k = self.k[self.solve_for[0]]
        s = self.s[self.solve_for[1]]

        # Pack a vector to give to the rate equations
        pack = np.array([a, a_star, k, s])

        # Define the rate equations
        da = self._da(*pack)
        da_star = self._da_star(*pack)
        db = self._db(*pack)
        dc = self._dc(*pack)

        # Return the rate equations
        return np.array([da, da_star, db, dc])
    
    def solve(self) -> np.ndarray:
        """
        Solves the reaction system and returns the result.

        Returns
        -------
        np.ndarray
            The solution to the reaction system.
        """
        sol = solve_ivp(self.rate_eq, (self.tau[0], self.tau[-1]), self.init_conc, t_eval=self.tau)
        return sol
    
    def _da(self, a, a_star, k, s):
        return 1/2*k*a*a_star - a**2 + 1/2*a**2
    
    def _da_star(self, a, a_star, k, s):
        return 1/2*a**2 - 1/2*k*a*a_star - s*k*a_star
    
    def _db(self, a, a_star, k, s):
        return 1/2*k*s*a_star
    
    def _dc(self, a, a_star, k, s):
        return 1/2*k*s*a_star

    def __repr__(self) -> str:
        """
        Returns a string representation of the BinaryReaction object.

        Returns
        -------
        str
            A string representation of the BinaryReaction object.
        """
        return f"BinaryReaction({self.p}, {self.q}, {self.r})"
    
    def __str__(self) -> str:
        """
        Returns a string representation of the BinaryReaction object.

        Returns
        -------
        str
            A string representation of the BinaryReaction object.
        """
        return f"BinaryReaction({self.p}, {self.q}, {self.r})"
    
    def __eq__(self, other) -> bool:
        """
        Returns True if the BinaryReaction object is equal to the other object.

        Returns
        -------
        bool
            True if the BinaryReaction object is equal to the other object.
        """
        if not isinstance(other, BinaryReaction):
            return False
        return (self.p == other.p and
                self.q == other.q and
                self.r == other.r)
    
    def __ne__(self, other) -> bool:
        """
        Returns True if the BinaryReaction object is not equal to the other object.

        Returns
        -------
        bool
            True if the BinaryReaction object is not equal to the other object.
        """
        return not self.__eq__(other)
    
    def __call__(self, t:float, y:Iterable[float]) -> Any:
        """
        Returns the rate equations for the BinaryReaction object.

        Parameters
        ----------
        t : float
            The time.
        y : Iterable[float]
            The state vector.

        Returns
        -------
        Any
            The rate equations for the BinaryReaction object.
        """
        return self.rate_eq(t, y)
    

class BinaryEquilibrium(BinaryReaction):
    def __init__(self, init_conc: Iterable[float], tau: Iterable[float], k: Iterable[float], s: Iterable[float]) -> None:
        super().__init__(init_conc, tau, k, s)

    def _da(self, a, a_star, k, s):
        return a**3/(a - s) - a**2

    def _da_star(self, a, a_star, k, s):
        return NotImplemented

    def _db(self, a, a_star, k, s):
        return s*a**2/(a - s)

    def _dc(self, a, a_star, k, s):
        return s*a**2/(a - s)

    def a_star(self, a_solution: Iterable[float]):
        """
        Returns the equilibrium value of a_star.
        NOT THE DERIVATIVE OF a_star.

        Returns
        -------
        float
            The equilibrium value of a_star.
        """
        a = a_solution

        # Get rate constants
        k = self.k[self.solve_for[0]]
        s = self.s[self.solve_for[1]]

        return a**2/(k*(a-s))
    
    def rate_eq(self, t:Iterable[float], y:Iterable[float]):
        """
        Define the rate equations for the BinaryReaction object.
        These are in a dimensionless form. In the case of the
        equilibrium approximation, the rate equations are different
        and we do not need to solve for a_star.
        """
        # Unpack the state vector
        a, b, c = y
        
        # Get rate constants
        k = self.k[self.solve_for[0]]
        s = self.s[self.solve_for[1]]

        # Pack a vector to give to the rate equations
        pack = np.array([a, None, k, s])

        # Define the rate equations
        da = self._da(*pack)
        db = self._db(*pack)
        dc = self._dc(*pack)

        # Return the rate equations
        return np.array([da, db, dc])

    def solve(self) -> np.ndarray:
        """
        Solves the reaction system and returns the result.

        Returns
        -------
        np.ndarray
            The solution to the reaction system.
        """
        x0 = [self.init_conc[0], self.init_conc[2], self.init_conc[3]]
        sol = solve_ivp(self.rate_eq, (self.tau[0], self.tau[-1]), x0, t_eval=self.tau)

        sol = [sol.y[0], self.a_star(sol.y[0]), sol.y[1], sol.y[2]]

        return sol

class BinarySingular(BinaryReaction):
    def __init__(self, init_conc: Iterable[float], tau: Iterable[float], k: Iterable[float], s: Iterable[float]) -> None:
        super().__init__(init_conc, tau, k, s)

    def sol_b(self, sol_a):

        a = sol_a
        # Get rate constants
        k = self.k[self.solve_for[0]]
        s = self.s[self.solve_for[1]]

        return s**3 * np.log((a-s)/(1-s)) + \
               s**2 * (a - 1) + \
               s/2 * (a**2 - 1)
    
    def sol_c(self, sol_a):
            
            a = sol_a
            # Get rate constants
            k = self.k[self.solve_for[0]]
            s = self.s[self.solve_for[1]]
    
            return s**3 * np.log((a-s)/(1-s)) + \
                s**2 * (a - 1) + \
                s/2 * (a**2 - 1)
    
    def sol_a_star(self, a_solution: Iterable[float]):
        """
        Returns the equilibrium value of a_star.
        NOT THE DERIVATIVE OF a_star.

        Returns
        -------
        float
            The equilibrium value of a_star.
        """
        a = a_solution

        # Get rate constants
        k = self.k[self.solve_for[0]]
        s = self.s[self.solve_for[1]]

        return a**2/(k*(a-s))
    
    def _da(self, t, a):
        # Get rate constants
        k = self.k[self.solve_for[0]]
        s = self.s[self.solve_for[1]]

        return a**3/(a - s) - a**2
    
    def sol_a(self):
        """
        Solves the reaction system and returns the result.

        Returns
        -------
        np.ndarray
            The solution to the reaction system.
        """
        y0 = [self.init_conc[0]]
        sol = solve_ivp(self._da, (self.tau[0], self.tau[-1]), y0, t_eval=self.tau)
        return sol
    
    def solve(self):
        sol_a = self.sol_a()
        sol = [sol_a.y[0], self.sol_a_star(sol_a.y[0]), self.sol_b(sol_a.y[0]), self.sol_c(sol_a.y[0])]
        return sol
    