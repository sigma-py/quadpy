from ..helpers import article

# ERR All schemes only have degree 1.
source = article(
    authors=["K. Sham Sunder", "R.A. Cookson"],
    title="Integration points for triangles and tetrahedrons obtained from the Gaussian quadrature points for a line",
    journal="Computers & Structures",
    year="1985",
    volume="21",
    number="5",
    pages="881-885",
    url="https://doi.org/10.1016/0045-7949%2885%2990198-1",
)


# def sunder_cookson_01():
#     weights, points = s3([1])
#     return T2Scheme("Sunder-Cookson 1", weights, points, 1, source)
#
#
# def sunder_cookson_02():
#     # ERR
#     warnings.warn("Sunder-Cookson claim degree 2, but the scheme is only degree 1.")
#     alpha = (1 - 1 / sqrt(3)) / 2
#     l2 = alpha / (alpha + 1)
#     weights, points = s2([1 / 3, l2])
#     return T2Scheme("Sunder-Cookson 2", weights, points, 1, source)
#
#
# def sunder_cookson_03():
#     # ERR
#     warnings.warn("Sunder-Cookson claim degree 3, but the scheme is only degree 1.")
#     weights, points = concat(s3(4 / 9), s2([5 / 27, 0.101286507323456]))
#     return T2Scheme("Sunder-Cookson 3", weights, points, 1, source)
#
#
# def sunder_cookson_04():
#     # ERR
#     warnings.warn("Sunder-Cookson claim degree 4, but the scheme is only degree 1.")
#     weights, points = concat(
#         s2(0.115951615045818, 0.064924047829079),
#         s2(0.217381718287515, 0.248125658963212),
#     )
#     return T2Scheme("Sunder-Cookson 4", weights, points, 1, source)
#
#
# def sunder_cookson_05():
#     # ERR
#     warnings.warn("Sunder-Cookson claim degree 5, but the scheme is only degree 1.")
#     weights, points = concat(
#         s3(64 / 225),
#         s2(0.078975628352063, 0.044808124460621),
#         s2(0.159542890166455, 0.187497434742173),
#     )
#     return T2Scheme("Sunder-Cookson 5", weights, points, 1, source)
