import cProfile
import io
import pstats


def profile_function(func, *args, **kwargs):
    """
    Profile a function using cProfile and return stats as string.
    """
    pr = cProfile.Profile()
    pr.enable()
    result = func(*args, **kwargs)
    pr.disable()
    s = io.StringIO()
    ps = pstats.Stats(pr, stream=s).sort_stats("cumulative")
    ps.print_stats(20)
    return s.getvalue(), result
