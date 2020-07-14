import json
import math
import pathlib
import warnings

import numpy

from ...helpers import article, fsd, untangle
from .._helpers import U3Scheme, cartesian_to_spherical, untangle2

source = article(
    authors=["Sangwoo Heo", "Yuan Xu"],
    title="Constructing Fully Symmetric Cubature Formulae for the Sphere",
    journal="Mathematics of Computation",
    volume="70",
    number="233",
    month="jan",
    year="2001",
    pages="269-279",
    url="https://doi.org/10.1090/S0025-5718-00-01198-4",
)


def _read(filename):
    this_dir = pathlib.Path(__file__).resolve().parent

    with open(this_dir / filename, "r") as f:
        data = json.load(f)

    name = data.pop("name")
    degree = data.pop("degree")
    tol = data.pop("test-tolerance")

    if tol > 1.0e-12:
        warnings.warn(f"The {name} scheme has low precision ({tol:.3e}).")

    points, weights = untangle2(data)
    theta_phi = cartesian_to_spherical(points)
    return U3Scheme(name, weights, points, theta_phi, degree, source, tol=tol)


def heo_xu_13():
    return _read("heo_xu_13.json")


def heo_xu_15():
    return _read("heo_xu_15.json")


def heo_xu_17():
    return _read("heo_xu_17.json")


def heo_xu_19a():
    return _read("heo_xu_19a.json")


def heo_xu_19b():
    return _read("heo_xu_19b.json")


def heo_xu_21a():
    return _read("heo_xu_21a.json")


def heo_xu_21b():
    return _read("heo_xu_21b.json")


def heo_xu_21c():
    return _read("heo_xu_21c.json")


def heo_xu_21d():
    return _read("heo_xu_21d.json")


def heo_xu_21e():
    return _read("heo_xu_21e.json")


def heo_xu_21f():
    return _read("heo_xu_21f.json")


def heo_xu_23_1():
    warnings.warn("The Heo-Xu schemes are only single-precision.")

    name = "Heo-Xu 23-1"
    degree = 23
    data = [
        (0.005026500922, _f((1.0, 1))),
        (0.005279416073, _f2(0.176588660459)),
        (0.003732271633, _f2(0.339207318490)),
        (0.006051284349, _f2(0.498904016243)),
        (0.005561610887, _f2(0.679838773734)),
        (0.005177363547, _f1(0.615520670749)),
        (0.005381929440, _f1(0.364554848325)),
        (0.004613082753, _f11(0.491903042583, 0.842732170863)),
    ]
    points, weights = untangle(data)
    theta_phi = cartesian_to_spherical(points)
    return U3Scheme(name, weights, points, theta_phi, degree, source, 1.299e-10)


def heo_xu_23_2():
    warnings.warn("The Heo-Xu schemes are only single-precision.")

    name = "Heo-Xu 23-2"
    degree = 23
    data = [
        (0.005651017861, _f((math.sqrt(1.0 / 3.0), 3))),
        (0.004259841569, _f2(0.115535209070)),
        (0.005294395887, _f2(0.282433358777)),
        (0.005588219406, _f2(0.441560530469)),
        (0.005591297404, _f2(0.670525624125)),
        (0.002936895883, _f2(0.706832372661)),
        (0.005051846065, _f1(0.345770219761)),
        (0.005530248916, _f11(0.525118572444, 0.836036015482)),
    ]
    points, weights = untangle(data)
    theta_phi = cartesian_to_spherical(points)
    return U3Scheme(name, weights, points, theta_phi, degree, source, 8.587e-11)


def heo_xu_23_3():
    warnings.warn("The Heo-Xu schemes are only single-precision.")

    name = "Heo-Xu 23-3"

    degree = 23
    data = [
        (-0.013079151392, _f((math.sqrt(1.0 / 3.0), 3))),
        (+0.002613177651, _f2(0.083820743273)),
        (+0.004696071564, _f2(0.208890425565)),
        (+0.010074474289, _f2(0.527146296056)),
        (+0.005747985286, _f2(0.684194032927)),
        (+0.005781262714, _f1(0.599358474983)),
        (+0.005394066441, _f1(0.364297979633)),
        (+0.005859672926, _f11(0.461647695180, 0.847134675079)),
    ]
    points, weights = untangle(data)
    theta_phi = cartesian_to_spherical(points)
    return U3Scheme(name, weights, points, theta_phi, degree, source, 1.222e-10)


def heo_xu_25_1():
    warnings.warn("The Heo-Xu schemes are only single-precision.")

    name = "Heo-Xu 25-1"
    degree = 25
    data = [
        (0.004313243133, _f((math.sqrt(1.0 / 3.0), 3))),
        (0.003986365505, _f2(0.111691690919)),
        (0.003663031548, _f2(0.315067166823)),
        (0.004204049922, _f2(0.459462014542)),
        (0.004269004376, _f2(0.660753497156)),
        (0.004203472415, _f2(0.702154945166)),
        (0.004142483118, _f1(0.532020255731)),
        (0.004090305599, _f11(0.519695051509, 0.822359911686)),
        (0.003789950437, _f11(0.329337385202, 0.938200966027)),
    ]
    points, weights = untangle(data)
    theta_phi = cartesian_to_spherical(points)
    return U3Scheme(name, weights, points, theta_phi, degree, source, 1.755e-10)


def heo_xu_25_2():
    warnings.warn("The Heo-Xu schemes are only single-precision.")

    name = "Heo-Xu 25-2"
    degree = 25
    data = [
        (0.003694297843, _f2(0.107086858755)),
        (0.003835709610, _f2(0.333879222938)),
        (0.004019086734, _f2(0.515654412063)),
        (0.003295936329, _f2(0.668811941305)),
        (0.004023268501, _f2(0.702834079289)),
        (0.003428775431, _f11(0.520513375926, 0.792782385892)),
        (0.003917073182, _f11(0.523649799991, 0.845078822417)),
        (0.004053335212, _f11(0.320409052387, 0.940317942977)),
    ]
    points, weights = untangle(data)
    theta_phi = cartesian_to_spherical(points)
    return U3Scheme(name, weights, points, theta_phi, degree, source, 1.266e-10)


def heo_xu_27_1():
    warnings.warn("The Heo-Xu schemes are only single-precision.")

    name = "Heo-Xu 27-1"
    degree = 27
    data = [
        (+0.004205508418, _f((math.sqrt(1.0 / 3.0), 3))),
        (+0.003927799571, _f2(0.110768319347)),
        (-0.000407112852, _f2(0.222696255452)),
        (+0.003694205329, _f2(0.322320222672)),
        (+0.004136341725, _f2(0.462107300704)),
        (+0.004202512176, _f2(0.660712667463)),
        (+0.004176738239, _f2(0.702450665599)),
        (+0.004229582701, _f1(0.525731112119)),
        (+0.004071467594, _f11(0.524493924092, 0.819343388819)),
        (+0.004080914226, _f11(0.323348454269, 0.939227929750)),
    ]
    points, weights = untangle(data)
    theta_phi = cartesian_to_spherical(points)
    return U3Scheme(name, weights, points, theta_phi, degree, source, 1.418e-10)


def heo_xu_27_2():
    warnings.warn("The Heo-Xu schemes are only single-precision.")

    name = "Heo-Xu 27-2"
    degree = 27
    data = [
        (+0.004145413998, _f((math.sqrt(1.0 / 3.0), 3))),
        (-0.001001399850, _f((1.0, 1))),
        (+0.004007770760, _f2(0.103271889407)),
        (+0.003609101265, _f2(0.326015048812)),
        (+0.004065377803, _f2(0.464476869175)),
        (+0.004167019227, _f2(0.661042838405)),
        (+0.004176827652, _f2(0.702427075066)),
        (+0.003963497360, _f11(0.523633133581, 0.818640861296)),
        (+0.002254418391, _f11(0.526385749330, 0.849564329470)),
        (+0.004036641877, _f11(0.324014265315, 0.939146832366)),
    ]
    points, weights = untangle(data)
    theta_phi = cartesian_to_spherical(points)
    return U3Scheme(name, weights, points, theta_phi, degree, source, 1.642e-10)


def heo_xu_27_3():
    warnings.warn("The Heo-Xu schemes are only single-precision.")

    name = "Heo-Xu 27-3"
    degree = 27
    data = [
        (0.003893829077, _f2(0.110332978624)),
        (0.003606286203, _f2(0.319075401552)),
        (0.003808504359, _f2(0.453117779552)),
        (0.002421634085, _f2(0.614431551407)),
        (0.004077606558, _f2(0.702545464357)),
        (0.004062279727, _f1(0.530118512908)),
        (0.002263516691, _f11(0.622283431805, 0.708196634681)),
        (0.003782934834, _f11(0.517277385525, 0.826495160269)),
        (0.003851811803, _f11(0.326096877477, 0.939054487386)),
    ]
    points, weights = untangle(data)
    theta_phi = cartesian_to_spherical(points)
    return U3Scheme(name, weights, points, theta_phi, degree, source, 1.419e-10)


def heo_xu_29():
    warnings.warn("The Heo-Xu schemes are only single-precision.")

    name = "Heo-Xu 29"
    degree = 29

    data = [
        (0.001894697146, _f2(0.072505573442)),
        (0.003783811492, _f2(0.310481161354)),
        (0.003012182159, _f2(0.440992773670)),
        (0.002724403361, _f2(0.527582710733)),
        (0.003559107906, _f2(0.659924119492)),
        (0.003227454716, _f2(0.699729990451)),
        (0.002024884516, _f11(0.594169250176, 0.802748431448)),
        (0.003467968868, _f11(0.531107498266, 0.806604859340)),
        (0.003229411196, _f11(0.426665999347, 0.898209890946)),
        (0.003010240364, _f11(0.247454899976, 0.963715742772)),
    ]
    points, weights = untangle(data)
    theta_phi = cartesian_to_spherical(points)
    return U3Scheme(name, weights, points, theta_phi, degree, source, 2.362e-10)


def heo_xu_31():
    warnings.warn("The Heo-Xu schemes are only single-precision.")

    name = "Heo-Xu 31"
    degree = 31
    data = [
        (0.000578329494, _f((math.sqrt(1.0 / 3.0), 3))),
        (0.003061522104, _f2(0.097855721318)),
        (0.002631322890, _f2(0.336734048041)),
        (0.002821112765, _f2(0.521197545604)),
        (0.002963046287, _f2(0.658338802723)),
        (0.002854170188, _f1(0.635957331845)),
        (0.002771653333, _f1(0.475363041502)),
        (0.002626546226, _f1(0.291089031268)),
        (0.002893976172, _f11(0.622288265573, 0.760613900822)),
        (0.002837749321, _f11(0.505117362658, 0.786669092004)),
        (0.002715804588, _f11(0.461391815178, 0.869489953690)),
        (0.002424728107, _f11(0.289362045905, 0.943297248354)),
    ]
    points, weights = untangle(data)
    theta_phi = cartesian_to_spherical(points)
    return U3Scheme(name, weights, points, theta_phi, degree, source, 2.324e-10)


def heo_xu_33():
    warnings.warn("The Heo-Xu schemes are only single-precision.")

    name = "Heo-Xu 33"
    degree = 33
    data = [
        (0.002848140682, _f((math.sqrt(1.0 / 3.0), 3))),
        (0.002449327062, _f2(0.087642362514)),
        (0.002179377451, _f2(0.244826259453)),
        (0.002653179507, _f2(0.373613468997)),
        (0.002817792197, _f2(0.483960471286)),
        (0.002832643891, _f2(0.648358940509)),
        (0.002776115721, _f2(0.692365152833)),
        (0.002634851612, _f1(0.424680986592)),
        (0.001044959201, _f11(0.663530717023, 0.747173970345)),
        (0.002552267148, _f11(0.560342562719, 0.821357965819)),
        (0.002771951327, _f11(0.544103582859, 0.784667458480)),
        (0.002540900745, _f11(0.408251323723, 0.892327016644)),
        (0.002276921078, _f11(0.259902436546, 0.962160309302)),
    ]
    points, weights = untangle(data)
    theta_phi = cartesian_to_spherical(points)
    return U3Scheme(name, weights, points, theta_phi, degree, source, 1.714e-10)


def heo_xu_35():
    warnings.warn("The Heo-Xu schemes are only single-precision.")

    name = "Heo-Xu 35"
    degree = 35
    data = [
        (0.002515482567, _f((math.sqrt(1.0 / 3.0), 3))),
        (0.001527515529, _f2(0.069156813118)),
        (0.002054028840, _f2(0.175148458557)),
        (0.002318417781, _f2(0.285287793163)),
        (0.002451618442, _f2(0.392405552644)),
        (0.002504293398, _f2(0.491306203394)),
        (0.002513606412, _f2(0.645641884498)),
        (0.002529886683, _f2(0.690921150829)),
        (0.001275574306, _f2(0.707101476221)),
        (0.002417442376, _f1(0.471598691154)),
        (0.001910951282, _f1(0.210272522872)),
        (0.002512236855, _f11(0.590515704894, 0.799927854385)),
        (0.002496644054, _f11(0.555015236112, 0.771746262687)),
        (0.002416930044, _f11(0.450233038264, 0.868946032283)),
        (0.002236607760, _f11(0.334436314543, 0.937180985852)),
    ]
    points, weights = untangle(data)
    theta_phi = cartesian_to_spherical(points)
    return U3Scheme(name, weights, points, theta_phi, degree, source, 2.104e-09)


def heo_xu_37():
    warnings.warn("The Heo-Xu schemes are only single-precision.")

    name = "Heo-Xu 37"
    degree = 37
    data = [
        (0.001436589472, _f((math.sqrt(1.0 / 3.0), 3))),
        (0.002233871811, _f2(0.181665204347)),
        (0.002119180525, _f2(0.303427242195)),
        (0.002281458727, _f2(0.483529149430)),
        (0.001864035223, _f2(0.625463680619)),
        (0.001858409063, _f2(0.705074619796)),
        (0.002336486555, _f1(0.444572576129)),
        (0.001818751796, _f11(0.630713490341, 0.746274388124)),
        (0.001961713367, _f11(0.596705492981, 0.724560936181)),
        (0.001611967438, _f11(0.587046661123, 0.807169311670)),
        (0.001942087580, _f11(0.505464222397, 0.844656624412)),
        (0.002245940979, _f11(0.458746299673, 0.824994538226)),
        (0.001967307858, _f11(0.369539839929, 0.915094022626)),
        (0.001561575961, _f11(0.277639225331, 0.958974697835)),
        (0.001137835823, _f11(0.116956662074, 0.992865475735)),
    ]
    points, weights = untangle(data)
    theta_phi = cartesian_to_spherical(points)
    return U3Scheme(name, weights, points, theta_phi, degree, source, 9.131e-10)


def heo_xu_39_1():
    warnings.warn("The Heo-Xu schemes are only single-precision.")

    name = "Heo-Xu 39-1"
    degree = 39
    data = [
        (0.001461069347, _f2(0.067102856429)),
        (0.002081064425, _f2(0.338906615896)),
        (0.002003883131, _f2(0.448234256905)),
        (0.001868396116, _f2(0.701808440747)),
        (0.001974704892, _f1(0.633404727273)),
        (0.001120631568, _f11(0.563818465348, 0.651676642607)),
        (0.001626513803, _f11(0.628884552131, 0.731765201351)),
        (0.001822889111, _f11(0.568134999733, 0.808896715135)),
        (0.001943309150, _f11(0.490729623093, 0.818248053726)),
        (0.001766607242, _f11(0.589988659798, 0.705342105582)),
        (0.001361185877, _f11(0.493911282058, 0.868009534637)),
        (0.001852624903, _f11(0.415275416293, 0.891987868825)),
        (0.001383169283, _f11(0.350457251849, 0.935213276436)),
        (0.001483291414, _f11(0.278497898816, 0.940831443705)),
        (0.001778552026, _f11(0.205495716318, 0.975518961862)),
    ]
    points, weights = untangle(data)
    theta_phi = cartesian_to_spherical(points)
    return U3Scheme(name, weights, points, theta_phi, degree, source, 2.291e-07)


def heo_xu_39b():
    return _read("heo_xu_39b.json")


def _f(*items):
    return fsd(3, *items)


def _f1(val):
    return fsd(3, (val, 1), (numpy.sqrt(1.0 - val ** 2), 1))


def _f2(val):
    return fsd(3, (val, 2), (numpy.sqrt(1.0 - 2 * val ** 2), 1))


def _f11(val1, val2):
    return fsd(3, (val1, 1), (val2, 1), (numpy.sqrt(1.0 - val1 ** 2 - val2 ** 2), 1))
