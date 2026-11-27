"""Errors and warnings that carry the MATLAB message identifiers.

The MATLAB toolbox raises errors such as
``error('GrowthPredict:fit_model:rBounds', ...)`` and its tests check the
identifier. The Python port keeps the identifier on the exception
(``exc.identifier``) and on the warning category message, so code and tests
can check it the same way.
"""

from __future__ import annotations

import warnings


class GrowthPredictError(Exception):
    """Error with a MATLAB-style identifier, e.g. ``GrowthPredict:get_Rt:type``."""

    def __init__(self, identifier, message=""):
        self.identifier = identifier
        self.message = message
        super().__init__("%s: %s" % (identifier, message) if message else identifier)


class GrowthPredictWarning(UserWarning):
    """Warning with a MATLAB-style identifier (``w.identifier``)."""

    def __init__(self, identifier, message=""):
        self.identifier = identifier
        self.message = message
        super().__init__("%s: %s" % (identifier, message) if message else identifier)


def warn(identifier, message, stacklevel=3):
    """Issue a :class:`GrowthPredictWarning` with the given identifier."""
    warnings.warn(GrowthPredictWarning(identifier, message), stacklevel=stacklevel)
