"""
This module defines some common decorators used in the comp aero UI
"""


from typing import Callable


def error_message_decorator(func) -> Callable:
    """
    Decorator that wraps a function in  try / catch statement and shows that error in a pop up
    message window in the UI

    Args:
        func (callable): The function to be wrapped

    Returns: The wrapped function
    """

    def _wrapped_func(self):
        result = None
        try:
            result = func(self)
        except ValueError as e:
            self.show_error(f"Error Occured! \n {e}")
        return result

    return _wrapped_func

