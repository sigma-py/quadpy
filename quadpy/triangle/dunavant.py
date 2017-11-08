# -*- coding: utf-8 -*-
#
from sympy import Rational

from .helpers import _s3, _s21, untangle3, weights_from_points
from ..helpers import untangle


class Dunavant(object):
    '''
    D.A. Dunavant,
    High Degree Efficient Symmetrical Gaussian Quadrature Rules for the
    Triangle,
    Article in International Journal for Numerical Methods in Engineering,
    21(6):1129-1148, June 1985,
    <https://doi.org/10.1002/nme.1620210612>.
    '''
    def __init__(self, index):
        self.name = 'Dunavant(%d)' % index
        if index == 1:
            self.degree = 1
            data = [(1, _s3())]
            self.bary, self.weights = untangle(data)
        elif index == 2:
            self.degree = 2
            data = [(Rational(1, 3), _s21(Rational(1, 6)))]
            self.bary, self.weights = untangle(data)
        elif index == 3:
            self.degree = 3
            data = [
                (-Rational(9, 16), _s3()),
                (Rational(25, 48), _s21(Rational(1, 5))),
                ]
            self.bary, self.weights = untangle(data)
        elif index == 4:
            self.degree = 4
            point_data = {
                's2': [
                    [0.445948490915965],
                    [0.091576213509771],
                    ],
                }
            weight_data = weights_from_points(point_data, self.degree)
            self.bary, self.weights = untangle3(point_data, weight_data)
        elif index == 5:
            self.degree = 5
            point_data = {
                's3': [[]],
                's2': [
                    [0.4701420641051],
                    [0.101286507323456],
                    ],
                }
            weight_data = weights_from_points(point_data, self.degree)
            self.bary, self.weights = untangle3(point_data, weight_data)
        elif index == 6:
            self.degree = 6
            point_data = {
                's2': [
                    [0.249286745170910],
                    [0.063089014491502],
                    ],
                's1': [
                    [0.053145049844817, 0.310352451033784],
                    ],
                }
            weight_data = weights_from_points(point_data, self.degree)
            self.bary, self.weights = untangle3(point_data, weight_data)
        elif index == 7:
            self.degree = 7
            point_data = {
                's3': [[]],
                's2': [
                    [0.260345966079040],
                    [0.065130102902216],
                    ],
                's1': [
                    [0.048690315425316, 0.312865496004874],
                    ],
                }
            weight_data = weights_from_points(point_data, self.degree)
            self.bary, self.weights = untangle3(point_data, weight_data)
        elif index == 8:
            self.degree = 8
            point_data = {
                's3': [[]],
                's2': [
                    [0.459292588292723],
                    [0.170569307751760],
                    [0.050547228317031],
                    ],
                's1': [
                    [0.008394777409958, 0.263112829634638],
                    ]
                }
            weight_data = weights_from_points(point_data, self.degree)
            self.bary, self.weights = untangle3(point_data, weight_data)
        elif index == 9:
            self.degree = 9
            point_data = {
                's3': [[]],
                's2': [
                    [0.489682519198738],
                    [0.437089591492937],
                    [0.188203535619033],
                    [0.044729513394453],
                    ],
                's1': [
                    [0.036838412054736, 0.221962989160766],
                    ],
                }
            weight_data = weights_from_points(point_data, self.degree)
            self.bary, self.weights = untangle3(point_data, weight_data)
        elif index == 10:
            self.degree = 10
            point_data = {
                's3': [[]],
                's2': [
                    [0.485577633383657],
                    [0.109481575485037],
                    ],
                's1': [
                    [0.141707219414880, 0.307939838764121],
                    [0.025003534762686, 0.246672560639903],
                    [0.009540815400299, 0.066803251012200],
                    ],
                }
            weight_data = weights_from_points(point_data, self.degree)
            self.bary, self.weights = untangle3(point_data, weight_data)
        elif index == 11:
            self.degree = 11
            point_data = {
                's2': [
                    [0.534611048270758],
                    [0.398969302965855],
                    [0.203309900431282],
                    [0.119350912282581],
                    [0.032364948111276],
                    ],
                's1': [
                    [0.050178138310495, 0.356620648261293],
                    [0.021022016536166, 0.171488980304042],
                    ],
                }
            weight_data = weights_from_points(point_data, self.degree)
            self.bary, self.weights = untangle3(point_data, weight_data)
        elif index == 12:
            self.degree = 12
            point_data = {
                's2': [
                    [0.488217389773805],
                    [0.439724392294460],
                    [0.271210385012116],
                    [0.127576145541586],
                    [0.021317350453210],
                    ],
                's1': [
                    [0.115343494534698, 0.275713269685514],
                    [0.022838332222257, 0.281325580989940],
                    [0.025734050548330, 0.116251915907597],
                    ],
                }
            weight_data = weights_from_points(point_data, self.degree)
            self.bary, self.weights = untangle3(point_data, weight_data)
        elif index == 13:
            self.degree = 13
            point_data = {
                's3': [[]],
                's2': [
                    [0.495048184939705],
                    [0.468716635109574],
                    [0.414521336801277],
                    [0.229399572042831],
                    [0.114424495196330],
                    [0.024811391363459],
                    ],
                's1': [
                    [0.094853828379579, 0.268794997058761],
                    [0.018100773278807, 0.291730066734288],
                    [0.022233076674090, 0.126357385491669],
                    ],
                }
            weight_data = weights_from_points(point_data, self.degree)
            self.bary, self.weights = untangle3(point_data, weight_data)
        elif index == 14:
            self.degree = 14
            point_data = {
                's2': [
                    [0.488963910362179],
                    [0.417644719340454],
                    [0.273477528308839],
                    [0.177205532412543],
                    [0.061799883090873],
                    [0.019390961248701],
                    ],
                's1': [
                    [0.057124757403648, 0.172266687821356],
                    [0.092916249356972, 0.336861459796345],
                    [0.014646950055654, 0.298372882136258],
                    [0.001268330932872, 0.118974497696957],
                    ],
                }
            weight_data = weights_from_points(point_data, self.degree)
            self.bary, self.weights = untangle3(point_data, weight_data)
        elif index == 15:
            self.degree = 15
            point_data = {
                's2': [
                    [0.506972916858243],
                    [0.431406354283023],
                    [0.277693644847144],
                    [0.126464891041254],
                    [0.070808385974686],
                    [0.018965170241073],
                    ],
                's1': [
                    [+0.133734161966621, 0.261311371140087],
                    [+0.036366677396917, 0.575586555512814],
                    [-0.010174883126571, 0.285712220049916],
                    [+0.036843869875878, 0.215599664072284],
                    [+0.012459809331199, 0.103575616576386],
                    ],
                }
            weight_data = weights_from_points(point_data, self.degree)
            self.bary, self.weights = untangle3(point_data, weight_data)
        elif index == 16:
            self.degree = 16
            point_data = {
                's3': [[]],
                's2': [
                    [0.497380541948438],
                    [0.413469438549352],
                    [0.470458599066991],
                    [0.240553749969521],
                    [0.147965794222573],
                    [0.075465187657474],
                    [0.016596402623025],
                    ],
                's1': [
                    [+0.103575692245252, 0.296555596579887],
                    [+0.020083411655416, 0.337723063403079],
                    [-0.004341002614139, 0.204748281642812],
                    [+0.041941786468010, 0.189358492130623],
                    [+0.014317320230681, 0.085283615682657],
                    ],
                }
            weight_data = weights_from_points(point_data, self.degree)
            self.bary, self.weights = untangle3(point_data, weight_data)
        elif index == 17:
            self.degree = 17
            point_data = {
                's3': [[]],
                's2': [
                    [0.497170540556774],
                    [0.482176322624625],
                    [0.450239969020782],
                    [0.400266239377397],
                    [0.252141267970953],
                    [0.162047004658461],
                    [0.075875882260746],
                    [0.015654726967822],
                    ],
                's1': [
                    [0.334319867363658, 0.655493203809423],
                    [0.292221537796944, 0.572337590532020],
                    [0.319574885423190, 0.626001190286228],
                    [0.190704224192292, 0.796427214974071],
                    [0.180483211648746, 0.752351005937729],
                    [0.080711313679564, 0.904625504095608],
                    ],
                }
            weight_data = weights_from_points(point_data, self.degree)
            self.bary, self.weights = untangle3(point_data, weight_data)
        elif index == 18:
            self.degree = 18
            point_data = {
                's3': [[]],
                's2': [
                    [0.493344808630921],
                    [0.469210594241957],
                    [0.436281395887006],
                    [0.394846170673416],
                    [0.249794568803157],
                    [0.161432193743843],
                    [0.076598227485371],
                    [0.024252439353450],
                    [0.043146367216965],
                    ],
                's1': [
                    [0.358911494940944, 0.632657968856636],
                    [0.294402476751957, 0.574410971510855],
                    [0.325017801641814, 0.624779046792512],
                    [0.184737559666046, 0.748933176523037],
                    [0.218796800013321, 0.769207005420443],
                    [0.101179597136408, 0.883962302273467],
                    [0.020874755282586, 1.014347260005363],
                    ],
                }
            weight_data = weights_from_points(point_data, self.degree)
            self.bary, self.weights = untangle3(point_data, weight_data)
        elif index == 19:
            self.degree = 19
            point_data = {
                's3': [[]],
                's2': [
                    [0.489609987073006],
                    [0.454536892697893],
                    [0.401416680649431],
                    [0.255551654403098],
                    [0.177077942152130],
                    [0.110061053227952],
                    [0.055528624251840],
                    [0.012621863777229],
                    ],
                's1': [
                    [0.395754787356943, 0.600633794794645],
                    [0.307929983880436, 0.557603261588784],
                    [0.264566948406520, 0.720987025817365],
                    [0.358539352205951, 0.594527068955871],
                    [0.157807405968595, 0.839331473680839],
                    [0.075050596975911, 0.701087978926173],
                    [0.142421601113383, 0.822931324069857],
                    [0.065494628082938, 0.924344252620784],
                    ],
                }
            weight_data = weights_from_points(point_data, self.degree)
            self.bary, self.weights = untangle3(point_data, weight_data)
        else:
            assert index == 20
            self.degree = 20
            point_data = {
                's3': [[]],
                's2': [
                    [0.500950464352200],
                    [0.488212957934729],
                    [0.455136681950283],
                    [0.401996259318289],
                    [0.255892909759421],
                    [0.176488255995106],
                    [0.104170855336758],
                    [0.053068963840930],
                    [0.041618715196029],
                    [0.011581921406822],
                    ],
                's1': [
                    [0.344855770229001, 0.606402646106160],
                    [0.377843269594854, 0.615842614456541],
                    [0.306635479062357, 0.559048000390295],
                    [0.249419362774742, 0.736606743262866],
                    [0.212775724802802, 0.711675142287434],
                    [0.146965436053239, 0.861402717154987],
                    [0.137726978828923, 0.835586957912363],
                    [0.059696109149007, 0.929756171556853],
                    ]
                }
            weight_data = weights_from_points(point_data, self.degree)
            self.bary, self.weights = untangle3(point_data, weight_data)

        # convert barycentric coordinates to reference triangle
        self.points = self.bary[:, 1:]
        return
