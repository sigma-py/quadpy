from .line_segment import integrate_adaptive


def quad(f, interval):
    return integrate_adaptive(f, interval, eps=1.0e-10)
