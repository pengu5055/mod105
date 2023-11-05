"""
Contains the BinaryReaction class, which is used to store information about a
binary chemical reaction. This class is used by the ReactionSystem class.
"""
import numpy as np
from typing import Any, Iterable

class BinaryReaction:
    def __init__(self,
                 p: float,
                 q: float,
                 r: float,
                 ) -> None:
        """
        Initializes a BinaryReaction object with the given parameters.

        A + A <--> A + A* (equilibrium) with rate constants p and q
        A* --> B + C (non-equilibrium) with rate constant r

        Parameters
        ----------
        p : float
            The rate constant for the forward direction
            of the equilibrium reaction.
        q : float
            The rate constant for the backward direction
            of the equilibrium reaction.
        r : float
            The rate constant for the forward direction
            of the non-equilibrium reaction.
        """
        self.p = p
        self.q = q
        self.r = r

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
        
    
    def rate_eq(self, t:Iterable[float], y:Iterable[float]):
        """
        Define the rate equations for the BinaryReaction object.
        """
        # Unpack the state vector
        A, A_star, B, C = y
        
        # Get rate constants
        p = self.p
        q = self.q
        r = self.r

        # Define the rate equations
        dAdt = self._dAdt(A, A_star, p, q)
        dA_stardt = self._dA_stardt(A, A_star, p, q, r)
        dBdt = self._dBdt(A_star, r)
        dCdt = self._dCdt(A_star, r)

        # Return the rate equations
        return np.array([dAdt, dA_stardt, dBdt, dCdt])
    
    def _dAdt(self, A, A_star, p, q):
        return -p*A**2 + q*A*A_star
    
    def _dA_stardt(self, A, A_star, p, q, r):
        return p*A**2 - q*A*A_star - r*A_star
    
    def _dBdt(self, A_star, r):
        return r*A_star
    
    def _dCdt(self, A_star, r):
        return r*A_star



