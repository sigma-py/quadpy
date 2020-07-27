import warnings

import numpy


class QuadratureScheme:
    def __init__(self, name, weights, points, degree, source, tol=1.0e-14):
        if tol > 1.0e-12:
            warnings.warn(f"The {name} scheme has low precision ({tol:.3e}).")

        self.test_tolerance = tol
        self.name = name
        self.degree = degree
        self.source = source

        # assert weights.shape[0] == points.shape[1], (
        #     f"Shape mismatch for {name}: "
        #     f"weights.shape = {weights.shape}, points.shape = {points.shape}"
        # )

        if weights.dtype == numpy.float64:
            self.weights = weights
        else:
            assert weights.dtype in [numpy.dtype("O"), numpy.int_]
            self.weights = weights.astype(numpy.float64)
            self.weights_symbolic = weights

        if points.dtype == numpy.float64:
            self.points = points
        else:
            assert points.dtype in [numpy.dtype("O"), numpy.int_]
            self.points = points.astype(numpy.float64)
            self.points_symbolic = points

    def savefig(self, filename, *args, **kwargs):
        from matplotlib import pyplot as plt

        self.plot(*args, **kwargs)
        # mpl keeps a hidden background patch that renders bbox_inches ineffective.
        # keep an eye out for https://stackoverflow.com/q/61712551/353337
        plt.savefig(filename, transparent=True, bbox_inches="tight", pad_inches=0)

    def show(self, *args, **kwargs):
        from matplotlib import pyplot as plt

        self.plot(*args, **kwargs)
        plt.show()

    def __str__(self):
        message = f"<quadrature scheme for {self.domain}>"
        message += f"\n  name:                 {self.name}"
        try:
            source = self.source
        except AttributeError:
            source = None

        if source:
            message += f"\n  source:               {source.title}"
            if source.authors:
                authors = ", ".join(source.authors)
                message += f"\n                        {authors}"

            string = []
            try:
                if source.journal:
                    string.append(source.journal)
            except AttributeError:
                pass
            try:
                if source.volume:
                    string.append(f"vol. {source.volume}")
            except AttributeError:
                pass
            try:
                if source.number:
                    string.append(f"no. {source.number}")
            except AttributeError:
                pass
            try:
                if source.pages:
                    string.append(f"pp. {source.pages}")
            except AttributeError:
                pass
            try:
                if source.year:
                    string.append(source.year)
            except AttributeError:
                pass

            if string:
                s = ", ".join(string)
                message += f"\n                        {s}"
            if source.url:
                message += f"\n                        {source.url}"

        message += f"\n  degree:               {self.degree}"
        message += f"\n  num points/weights:   {len(self.weights)}"

        try:
            pi = self.points_inside()
        except AttributeError:
            pi = None
        if pi:
            if pi.all():
                point_position = "strictly inside"
            elif self.points_inside_or_boundary().all():
                point_position = "inside + boundary"
            else:
                point_position = "outside"
            message += f"\n  point position:       {point_position}"

        weights_positive = all(self.weights > 0)
        message += f"\n  all weights positive: {weights_positive}"
        return message
