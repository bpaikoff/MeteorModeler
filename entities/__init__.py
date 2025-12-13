# entities/__init__.py
"""Data structures and classes for meteor simulation."""

from entities.meteor import Meteor
from entities.fragment import FragmentModel
from entities.particles import Particle
from entities.compositions import Compositions

__all__ = [
    "Meteor",
    "FragmentModel",
    "Particle",
    "Compositions",
]
