from typing import Callable

def error_message_decorator(func) -> Callable:
    
    def _wrapped_func(self):
        result = None
        try:
            result = func(self)
        except ValueError as e:
            self.show_error(f"Error Occured! \n {e}")
        return result
    return _wrapped_func