"""Custom exceptions for your_project."""

class ShapleyError(Exception):
    """Base exception for all your_project errors.

    Users can catch this to handle any library-specific error.
    """
    pass


class InvalidDataError(ShapleyError):
    """Raised when input data is invalid or cannot be parsed."""
    pass

class ComputationError(ShapleyError):
    """Raised when a calculation fails or produces invalid results."""
    pass

class ConfigError(ShapleyError):
    """Raised when configuration is missing or incorrect."""
    pass

class InforError(ShapleyError):
    """Raised when there is an error related to information processing."""
    pass

class AmbiguousFormulaError(ShapleyError):
    """Raised when a boolean formula is ambiguous or cannot be interpreted."""
    pass

class NetworkError(ShapleyError):
    """Raised when there is an error related to network processing."""
    pass

class UnknownError(ShapleyError):
    """Raised for unexpected errors."""
    pass 