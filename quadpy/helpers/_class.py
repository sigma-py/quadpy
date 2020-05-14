class QuadratureScheme:
    def __init__(self):
        return

    def __str__(self):
        message = f"<quadrature scheme for {self.domain}>"
        message += f"\n  name:                 {self.name}"
        if self.source:
            message += f"\n  source:               {self.source.title}"
            if self.source.authors:
                authors = ", ".join(self.source.authors)
                message += f"\n                        {authors}"

            string = []
            try:
                if self.source.journal:
                    string.append(self.source.journal)
            except AttributeError:
                pass
            try:
                if self.source.volume:
                    string.append(f"vol. {self.source.volume}")
            except AttributeError:
                pass
            try:
                if self.source.number:
                    string.append(f"no. {self.source.number}")
            except AttributeError:
                pass
            try:
                if self.source.pages:
                    string.append(f"pp. {self.source.pages}")
            except AttributeError:
                pass
            try:
                if self.source.year:
                    string.append(self.source.year)
            except AttributeError:
                pass
            if string:
                s = ", ".join(string)
                message += f"\n                        {s}"
            if self.source.url:
                message += f"\n                        {self.source.url}"

        message += f"\n  degree:               {self.degree}"
        message += f"\n  num points/weights:   {len(self.weights)}"
        if self.points_inside():
            point_position = "strictly inside"
        elif self.points_inside_or_boundary():
            point_position = "inside + boundary"
        else:
            point_position = "outside"
        message += f"\n  point position:       {point_position}"
        weights_positive = all(self.weights > 0)
        message += f"\n  all weights positive: {weights_positive}"
        return message
