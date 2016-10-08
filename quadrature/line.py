# -*- coding: utf-8 -*-
#
import math
import numpy
import sympy


def integrate(f, a, b, rule):
    out = math.fsum([
        weight * f(0.5 * (point + 1) * (b-a) + a)
        for point, weight in zip(rule.points, rule.weights)
        ])
    return 0.5 * (b - a) * out


def show(a, b, scheme, render=True):
    from matplotlib import pyplot as plt
    pts = 0.5 * (scheme.points + 1) * (b-a) + a
    plt.plot([0.0, 1.0], [0.0, 0.0], '-k')
    plt.bar(
        pts, scheme.weights,
        color='b',
        alpha=0.5,
        width=(b-a) / len(scheme.weights)
        )
    # plt.axis('equal')
    if render:
        plt.show()
    return


class Midpoint(object):
    def __init__(self):
        self.weights = [2.0]
        self.points = numpy.array([
            0.0
            ])
        self.degree = 1
        return


class Trapezoidal(object):
    def __init__(self):
        self.weights = [1.0, 1.0]
        self.points = numpy.array([
            -1.0,
            1.0
            ])
        self.degree = 1
        return


class GaussLegendre(object):
    '''
    Gauß-Legendre quadrature.
    '''
    def __init__(self, order):
        self.degree = 2*order - 1
        self.points, self.weights = numpy.polynomial.legendre.leggauss(order)
        return


class GaussPatterson(object):
    '''
    Gauß-Patterson quadrature.
    <https://people.sc.fsu.edu/~jburkardt/datasets/quadrature_rules_patterson/quadrature_rules_patterson.html>
    '''
    def __init__(self, index):
        self.degree = 3*2**index - 1
        if index == 0:
            self.degree = 1  # override degree
            self.weights = [2.0]
            self.points = numpy.array([
                0.0
                ])
        elif index == 1:
            self.weights = numpy.array([
                5.0 / 9.0,
                8.0 / 9.0,
                5.0 / 9.0,
                ])
            self.points = numpy.array([
                -0.7745966692414834,
                0.0,
                0.7745966692414834,
                ])
        elif index == 2:
            self.weights = numpy.array([
                0.1046562260264673,
                0.2684880898683334,
                0.4013974147759622,
                0.4509165386584741,
                0.4013974147759622,
                0.2684880898683334,
                0.1046562260264673,
                ])
            self.points = numpy.array([
                -0.9604912687080203,
                -0.7745966692414834,
                -0.4342437493468025,
                0.0,
                0.4342437493468025,
                0.7745966692414834,
                0.9604912687080203,
                ])
        elif index == 3:
            self.weights = numpy.array([
                0.1700171962994026E-01,
                0.5160328299707974E-01,
                0.9292719531512454E-01,
                0.1344152552437842,
                0.1715119091363914,
                0.2006285293769890,
                0.2191568584015875,
                0.2255104997982067,
                0.2191568584015875,
                0.2006285293769890,
                0.1715119091363914,
                0.1344152552437842,
                0.9292719531512454E-01,
                0.5160328299707974E-01,
                0.1700171962994026E-01,
                ])
            self.points = numpy.array([
                -0.9938319632127550,
                -0.9604912687080203,
                -0.8884592328722570,
                -0.7745966692414834,
                -0.6211029467372264,
                -0.4342437493468025,
                -0.2233866864289669,
                0.0,
                0.2233866864289669,
                0.4342437493468025,
                0.6211029467372264,
                0.7745966692414834,
                0.8884592328722570,
                0.9604912687080203,
                0.9938319632127550,
                ])
        elif index == 4:
            self.weights = numpy.array([
                0.2544780791561875E-02,
                0.8434565739321106E-02,
                0.1644604985438781E-01,
                0.2580759809617665E-01,
                0.3595710330712932E-01,
                0.4646289326175799E-01,
                0.5697950949412336E-01,
                0.6720775429599070E-01,
                0.7687962049900353E-01,
                0.8575592004999034E-01,
                0.9362710998126447E-01,
                0.1003142786117956,
                0.1056698935802348,
                0.1095784210559246,
                0.1119568730209535,
                0.1127552567207687,
                0.1119568730209535,
                0.1095784210559246,
                0.1056698935802348,
                0.1003142786117956,
                0.9362710998126447E-01,
                0.8575592004999034E-01,
                0.7687962049900353E-01,
                0.6720775429599070E-01,
                0.5697950949412336E-01,
                0.4646289326175799E-01,
                0.3595710330712932E-01,
                0.2580759809617665E-01,
                0.1644604985438781E-01,
                0.8434565739321106E-02,
                0.2544780791561875E-02,
                ])
            self.points = numpy.array([
                -0.9990981249676676,
                -0.9938319632127550,
                -0.9815311495537401,
                -0.9604912687080203,
                -0.9296548574297401,
                -0.8884592328722570,
                -0.8367259381688688,
                -0.7745966692414834,
                -0.7024962064915271,
                -0.6211029467372264,
                -0.5313197436443756,
                -0.4342437493468025,
                -0.3311353932579768,
                -0.2233866864289669,
                -0.1124889431331866,
                0.0,
                0.1124889431331866,
                0.2233866864289669,
                0.3311353932579768,
                0.4342437493468025,
                0.5313197436443756,
                0.6211029467372264,
                0.7024962064915271,
                0.7745966692414834,
                0.8367259381688688,
                0.8884592328722570,
                0.9296548574297401,
                0.9604912687080203,
                0.9815311495537401,
                0.9938319632127550,
                0.9990981249676676,
                ])
        elif index == 5:
            self.weights = numpy.array([
                0.3632214818455306E-03,
                0.1265156556230068E-02,
                0.2579049794685688E-02,
                0.4217630441558855E-02,
                0.6115506822117246E-02,
                0.8223007957235930E-02,
                0.1049824690962132E-01,
                0.1290380010035127E-01,
                0.1540675046655950E-01,
                0.1797855156812827E-01,
                0.2059423391591271E-01,
                0.2323144663991027E-01,
                0.2586967932721475E-01,
                0.2848975474583355E-01,
                0.3107355111168797E-01,
                0.3360387714820773E-01,
                0.3606443278078257E-01,
                0.3843981024945553E-01,
                0.4071551011694432E-01,
                0.4287796002500773E-01,
                0.4491453165363220E-01,
                0.4681355499062801E-01,
                0.4856433040667320E-01,
                0.5015713930589954E-01,
                0.5158325395204846E-01,
                0.5283494679011652E-01,
                0.5390549933526606E-01,
                0.5478921052796287E-01,
                0.5548140435655936E-01,
                0.5597843651047632E-01,
                0.5627769983125430E-01,
                0.5637762836038471E-01,
                0.5627769983125430E-01,
                0.5597843651047632E-01,
                0.5548140435655936E-01,
                0.5478921052796287E-01,
                0.5390549933526606E-01,
                0.5283494679011652E-01,
                0.5158325395204846E-01,
                0.5015713930589954E-01,
                0.4856433040667320E-01,
                0.4681355499062801E-01,
                0.4491453165363220E-01,
                0.4287796002500773E-01,
                0.4071551011694432E-01,
                0.3843981024945553E-01,
                0.3606443278078257E-01,
                0.3360387714820773E-01,
                0.3107355111168797E-01,
                0.2848975474583355E-01,
                0.2586967932721475E-01,
                0.2323144663991027E-01,
                0.2059423391591271E-01,
                0.1797855156812827E-01,
                0.1540675046655950E-01,
                0.1290380010035127E-01,
                0.1049824690962132E-01,
                0.8223007957235930E-02,
                0.6115506822117246E-02,
                0.4217630441558855E-02,
                0.2579049794685688E-02,
                0.1265156556230068E-02,
                0.3632214818455306E-03,
                ])
            self.points = numpy.array([
                -0.9998728881203576,
                -0.9990981249676676,
                -0.9972062593722220,
                -0.9938319632127550,
                -0.9886847575474295,
                -0.9815311495537401,
                -0.9721828747485818,
                -0.9604912687080203,
                -0.9463428583734029,
                -0.9296548574297401,
                -0.9103711569570043,
                -0.8884592328722570,
                -0.8639079381936905,
                -0.8367259381688688,
                -0.8069405319502176,
                -0.7745966692414834,
                -0.7397560443526947,
                -0.7024962064915271,
                -0.6629096600247806,
                -0.6211029467372264,
                -0.5771957100520458,
                -0.5313197436443756,
                -0.4836180269458411,
                -0.4342437493468025,
                -0.3833593241987304,
                -0.3311353932579768,
                -0.2777498220218243,
                -0.2233866864289669,
                -0.1682352515522075,
                -0.1124889431331866,
                -0.5634431304659279E-01,
                0.0,
                0.5634431304659279E-01,
                0.1124889431331866,
                0.1682352515522075,
                0.2233866864289669,
                0.2777498220218243,
                0.3311353932579768,
                0.3833593241987304,
                0.4342437493468025,
                0.4836180269458411,
                0.5313197436443756,
                0.5771957100520458,
                0.6211029467372264,
                0.6629096600247806,
                0.7024962064915271,
                0.7397560443526947,
                0.7745966692414834,
                0.8069405319502176,
                0.8367259381688688,
                0.8639079381936905,
                0.8884592328722570,
                0.9103711569570043,
                0.9296548574297401,
                0.9463428583734029,
                0.9604912687080203,
                0.9721828747485818,
                0.9815311495537401,
                0.9886847575474295,
                0.9938319632127550,
                0.9972062593722220,
                0.9990981249676676,
                0.9998728881203576,
                ])
        elif index == 6:
            self.weights = numpy.array([
                0.5053609520786252E-04,
                0.1807395644453884E-03,
                0.3777466463269846E-03,
                0.6326073193626335E-03,
                0.9383698485423815E-03,
                0.1289524082610417E-02,
                0.1681142865421470E-02,
                0.2108815245726633E-02,
                0.2568764943794020E-02,
                0.3057753410175531E-02,
                0.3572892783517299E-02,
                0.4111503978654693E-02,
                0.4671050372114322E-02,
                0.5249123454808859E-02,
                0.5843449875835640E-02,
                0.6451900050175737E-02,
                0.7072489995433555E-02,
                0.7703375233279742E-02,
                0.8342838753968157E-02,
                0.8989275784064136E-02,
                0.9641177729702537E-02,
                0.1029711695795636E-01,
                0.1095573338783790E-01,
                0.1161572331995513E-01,
                0.1227583056008277E-01,
                0.1293483966360737E-01,
                0.1359157100976555E-01,
                0.1424487737291678E-01,
                0.1489364166481518E-01,
                0.1553677555584398E-01,
                0.1617321872957772E-01,
                0.1680193857410386E-01,
                0.1742193015946417E-01,
                0.1803221639039129E-01,
                0.1863184825613879E-01,
                0.1921990512472777E-01,
                0.1979549504809750E-01,
                0.2035775505847216E-01,
                0.2090585144581202E-01,
                0.2143898001250387E-01,
                0.2195636630531782E-01,
                0.2245726582681610E-01,
                0.2294096422938775E-01,
                0.2340677749531401E-01,
                0.2385405210603854E-01,
                0.2428216520333660E-01,
                0.2469052474448768E-01,
                0.2507856965294977E-01,
                0.2544576996546477E-01,
                0.2579162697602423E-01,
                0.2611567337670610E-01,
                0.2641747339505826E-01,
                0.2669662292745036E-01,
                0.2695274966763303E-01,
                0.2718551322962479E-01,
                0.2739460526398143E-01,
                0.2757974956648187E-01,
                0.2774070217827968E-01,
                0.2787725147661370E-01,
                0.2798921825523816E-01,
                0.2807645579381725E-01,
                0.2813884991562715E-01,
                0.2817631903301660E-01,
                0.2818881418019236E-01,
                0.2817631903301660E-01,
                0.2813884991562715E-01,
                0.2807645579381725E-01,
                0.2798921825523816E-01,
                0.2787725147661370E-01,
                0.2774070217827968E-01,
                0.2757974956648187E-01,
                0.2739460526398143E-01,
                0.2718551322962479E-01,
                0.2695274966763303E-01,
                0.2669662292745036E-01,
                0.2641747339505826E-01,
                0.2611567337670610E-01,
                0.2579162697602423E-01,
                0.2544576996546477E-01,
                0.2507856965294977E-01,
                0.2469052474448768E-01,
                0.2428216520333660E-01,
                0.2385405210603854E-01,
                0.2340677749531401E-01,
                0.2294096422938775E-01,
                0.2245726582681610E-01,
                0.2195636630531782E-01,
                0.2143898001250387E-01,
                0.2090585144581202E-01,
                0.2035775505847216E-01,
                0.1979549504809750E-01,
                0.1921990512472777E-01,
                0.1863184825613879E-01,
                0.1803221639039129E-01,
                0.1742193015946417E-01,
                0.1680193857410386E-01,
                0.1617321872957772E-01,
                0.1553677555584398E-01,
                0.1489364166481518E-01,
                0.1424487737291678E-01,
                0.1359157100976555E-01,
                0.1293483966360737E-01,
                0.1227583056008277E-01,
                0.1161572331995513E-01,
                0.1095573338783790E-01,
                0.1029711695795636E-01,
                0.9641177729702537E-02,
                0.8989275784064136E-02,
                0.8342838753968157E-02,
                0.7703375233279742E-02,
                0.7072489995433555E-02,
                0.6451900050175737E-02,
                0.5843449875835640E-02,
                0.5249123454808859E-02,
                0.4671050372114322E-02,
                0.4111503978654693E-02,
                0.3572892783517299E-02,
                0.3057753410175531E-02,
                0.2568764943794020E-02,
                0.2108815245726633E-02,
                0.1681142865421470E-02,
                0.1289524082610417E-02,
                0.9383698485423815E-03,
                0.6326073193626335E-03,
                0.3777466463269846E-03,
                0.1807395644453884E-03,
                0.5053609520786252E-04,
                ])
            self.points = numpy.array([
                -0.9999824303548916,
                -0.9998728881203576,
                -0.9995987996719107,
                -0.9990981249676676,
                -0.9983166353184074,
                -0.9972062593722220,
                -0.9957241046984072,
                -0.9938319632127550,
                -0.9914957211781061,
                -0.9886847575474295,
                -0.9853714995985203,
                -0.9815311495537401,
                -0.9771415146397057,
                -0.9721828747485818,
                -0.9666378515584165,
                -0.9604912687080203,
                -0.9537300064257611,
                -0.9463428583734029,
                -0.9383203977795929,
                -0.9296548574297401,
                -0.9203400254700124,
                -0.9103711569570043,
                -0.8997448997769401,
                -0.8884592328722570,
                -0.8765134144847053,
                -0.8639079381936905,
                -0.8506444947683502,
                -0.8367259381688688,
                -0.8221562543649804,
                -0.8069405319502176,
                -0.7910849337998483,
                -0.7745966692414834,
                -0.7574839663805136,
                -0.7397560443526947,
                -0.7214230853700989,
                -0.7024962064915271,
                -0.6829874310910792,
                -0.6629096600247806,
                -0.6422766425097595,
                -0.6211029467372264,
                -0.5994039302422429,
                -0.5771957100520458,
                -0.5544951326319325,
                -0.5313197436443756,
                -0.5076877575337166,
                -0.4836180269458411,
                -0.4591300119898323,
                -0.4342437493468025,
                -0.4089798212298887,
                -0.3833593241987304,
                -0.3574038378315322,
                -0.3311353932579768,
                -0.3045764415567140,
                -0.2777498220218243,
                -0.2506787303034832,
                -0.2233866864289669,
                -0.1958975027111002,
                -0.1682352515522075,
                -0.1404242331525602,
                -0.1124889431331866,
                -0.8445404008371088E-01,
                -0.5634431304659279E-01,
                -0.2818464894974569E-01,
                0.0,
                0.2818464894974569E-01,
                0.5634431304659279E-01,
                0.8445404008371088E-01,
                0.1124889431331866,
                0.1404242331525602,
                0.1682352515522075,
                0.1958975027111002,
                0.2233866864289669,
                0.2506787303034832,
                0.2777498220218243,
                0.3045764415567140,
                0.3311353932579768,
                0.3574038378315322,
                0.3833593241987304,
                0.4089798212298887,
                0.4342437493468025,
                0.4591300119898323,
                0.4836180269458411,
                0.5076877575337166,
                0.5313197436443756,
                0.5544951326319325,
                0.5771957100520458,
                0.5994039302422429,
                0.6211029467372264,
                0.6422766425097595,
                0.6629096600247806,
                0.6829874310910792,
                0.7024962064915271,
                0.7214230853700989,
                0.7397560443526947,
                0.7574839663805136,
                0.7745966692414834,
                0.7910849337998483,
                0.8069405319502176,
                0.8221562543649804,
                0.8367259381688688,
                0.8506444947683502,
                0.8639079381936905,
                0.8765134144847053,
                0.8884592328722570,
                0.8997448997769401,
                0.9103711569570043,
                0.9203400254700124,
                0.9296548574297401,
                0.9383203977795929,
                0.9463428583734029,
                0.9537300064257611,
                0.9604912687080203,
                0.9666378515584165,
                0.9721828747485818,
                0.9771415146397057,
                0.9815311495537401,
                0.9853714995985203,
                0.9886847575474295,
                0.9914957211781061,
                0.9938319632127550,
                0.9957241046984072,
                0.9972062593722220,
                0.9983166353184074,
                0.9990981249676676,
                0.9995987996719107,
                0.9998728881203576,
                0.9999824303548916,
                ])
        else:
            raise ValueError('Illegal Gauss-Patterson order')


class ClenshawCurtis(object):
    '''
    Clenshaw-Curtis quadrature.
    <https://people.sc.fsu.edu/~jburkardt/datasets/quadrature_rules_clenshaw_curtis/quadrature_rules_clenshaw_curtis.html>
    '''
    def __init__(self, order):
        if order == 1:
            self.weights = [2.0]
            self.points = numpy.array([
                0.0
                ])
            self.degree = 1
        elif order == 2:
            self.weights = numpy.array([
                1.0,
                1.0
                ])
            self.points = numpy.array([
                -1.0,
                1.0
                ])
            self.degree = 1
        elif order == 3:
            self.weights = numpy.array([
                1.0/3.0,
                4.0/3.0,
                1.0/3.0,
                ])
            self.points = numpy.array([
                -1.0,
                0.0,
                1.0
                ])
            self.degree = 3
        elif order == 4:
            self.weights = numpy.array([
                1.0/9.0,
                8.0/9.0,
                8.0/9.0,
                1.0/9.0,
                ])
            self.points = numpy.array([
                -1.0,
                -0.5,
                0.5,
                1.0,
                ])
            self.degree = 3
        elif order == 5:
            self.weights = numpy.array([
                1.0/15.0,
                8.0/15.0,
                0.8,
                8.0/15.0,
                1.0/15.0,
                ])
            self.points = numpy.array([
                -1.0,
                -numpy.sqrt(0.5),
                0.0,
                numpy.sqrt(0.5),
                1.0
                ])
            self.degree = 5
        elif order == 9:
            self.weights = numpy.array([
                0.1587301587301588E-01,
                0.1462186492160182,
                0.2793650793650794,
                0.3617178587204897,
                0.3936507936507936,
                0.3617178587204898,
                0.2793650793650794,
                0.1462186492160182,
                0.1587301587301588E-01,
                ])
            self.points = numpy.array([
                -1.0,
                -0.9238795325112867,
                -numpy.sqrt(0.5),
                -0.3826834323650897,
                0.0,
                0.3826834323650898,
                numpy.sqrt(0.5),
                0.9238795325112867,
                1.0,
                ])
            self.degree = 9
        elif order == 17:
            self.weights = numpy.array([
                0.3921568627450983E-02,
                0.3736870283720560E-01,
                0.7548233154315186E-01,
                0.1089055525818909,
                0.1389564683682331,
                0.1631726642817033,
                0.1814737842364934,
                0.1925138646129257,
                0.1964101258218905,
                0.1925138646129257,
                0.1814737842364934,
                0.1631726642817033,
                0.1389564683682331,
                0.1089055525818909,
                0.7548233154315184E-01,
                0.3736870283720570E-01,
                0.3921568627450983E-02,
                ])
            self.points = numpy.array([
                -1.0000000000000000,
                -0.9807852804032304,
                -0.9238795325112867,
                -0.8314696123025453,
                -0.7071067811865475,
                -0.5555702330196020,
                -0.3826834323650897,
                -0.1950903220161282,
                0.00000000000000,
                0.1950903220161283,
                0.3826834323650898,
                0.5555702330196023,
                0.7071067811865475,
                0.8314696123025452,
                0.9238795325112867,
                0.9807852804032304,
                1.0000000000000000,
                ])
            self.degree = 17
        elif order == 33:
            self.weights = numpy.array([
                0.9775171065493659E-03,
                0.9393197962955013E-02,
                0.1923424513268114E-01,
                0.2845791667723369E-01,
                0.3759434191404722E-01,
                0.4626276283775175E-01,
                0.5455501630398032E-01,
                0.6227210954529399E-01,
                0.6942757563043547E-01,
                0.7588380044138848E-01,
                0.8163481765493850E-01,
                0.8657753844182743E-01,
                0.9070611286772098E-01,
                0.9394324443876872E-01,
                0.9629232594548820E-01,
                0.9769818820805558E-01,
                0.9817857778176831E-01,
                0.9769818820805558E-01,
                0.9629232594548819E-01,
                0.9394324443876871E-01,
                0.9070611286772098E-01,
                0.8657753844182743E-01,
                0.8163481765493850E-01,
                0.7588380044138850E-01,
                0.6942757563043547E-01,
                0.6227210954529399E-01,
                0.5455501630398031E-01,
                0.4626276283775177E-01,
                0.3759434191404721E-01,
                0.2845791667723370E-01,
                0.1923424513268119E-01,
                0.9393197962955048E-02,
                0.9775171065493659E-03,
                ])
            self.points = numpy.array([
                -1.0000000000000000,
                -0.9951847266721968,
                -0.9807852804032304,
                -0.9569403357322088,
                -0.9238795325112867,
                -0.8819212643483549,
                -0.8314696123025453,
                -0.7730104533627370,
                -0.7071067811865475,
                -0.6343932841636454,
                -0.5555702330196020,
                -0.4713967368259977,
                -0.3826834323650897,
                -0.2902846772544622,
                -0.1950903220161282,
                -0.9801714032956065E-01,
                0.000000000000000,
                0.9801714032956077E-01,
                0.1950903220161283,
                0.2902846772544623,
                0.3826834323650898,
                0.4713967368259978,
                0.5555702330196023,
                0.6343932841636455,
                0.7071067811865475,
                0.7730104533627370,
                0.8314696123025452,
                0.8819212643483550,
                0.9238795325112867,
                0.9569403357322088,
                0.9807852804032304,
                0.9951847266721969,
                1.0000000000000000,
                ])
            self.degree = 33
        elif order == 65:
            self.weights = numpy.array([
                0.2442002442002449E-03,
                0.2351490675311702E-02,
                0.4831465448790911E-02,
                0.7192693161736115E-02,
                0.9582338795283791E-02,
                0.1192339471421277E-01,
                0.1425206043235199E-01,
                0.1653498765728959E-01,
                0.1878652974179578E-01,
                0.2098627442973744E-01,
                0.2314069493435819E-01,
                0.2523506498175476E-01,
                0.2727225714146840E-01,
                0.2924065319746835E-01,
                0.3114129710406762E-01,
                0.3296454656997634E-01,
                0.3471049818092511E-01,
                0.3637092028663919E-01,
                0.3794545992128482E-01,
                0.3942698871295609E-01,
                0.4081501340035782E-01,
                0.4210333111141810E-01,
                0.4329151496169082E-01,
                0.4437417923925731E-01,
                0.4535110955166067E-01,
                0.4621766751092559E-01,
                0.4697395904661415E-01,
                0.4761604458525018E-01,
                0.4814443257251221E-01,
                0.4855584485714105E-01,
                0.4885125664306610E-01,
                0.4902801843102554E-01,
                0.4908762351494248E-01,
                0.4902801843102556E-01,
                0.4885125664306610E-01,
                0.4855584485714105E-01,
                0.4814443257251221E-01,
                0.4761604458525019E-01,
                0.4697395904661414E-01,
                0.4621766751092559E-01,
                0.4535110955166067E-01,
                0.4437417923925733E-01,
                0.4329151496169082E-01,
                0.4210333111141811E-01,
                0.4081501340035782E-01,
                0.3942698871295609E-01,
                0.3794545992128483E-01,
                0.3637092028663919E-01,
                0.3471049818092511E-01,
                0.3296454656997635E-01,
                0.3114129710406762E-01,
                0.2924065319746836E-01,
                0.2727225714146839E-01,
                0.2523506498175477E-01,
                0.2314069493435821E-01,
                0.2098627442973743E-01,
                0.1878652974179578E-01,
                0.1653498765728961E-01,
                0.1425206043235200E-01,
                0.1192339471421278E-01,
                0.9582338795283809E-02,
                0.7192693161736120E-02,
                0.4831465448790926E-02,
                0.2351490675311698E-02,
                0.2442002442002449E-03,
                ])
            self.points = numpy.array([
                -1.0000000000000000,
                -0.9987954562051724,
                -0.9951847266721968,
                -0.9891765099647810,
                -0.9807852804032304,
                -0.9700312531945440,
                -0.9569403357322088,
                -0.9415440651830207,
                -0.9238795325112867,
                -0.9039892931234433,
                -0.8819212643483549,
                -0.8577286100002720,
                -0.8314696123025453,
                -0.8032075314806448,
                -0.7730104533627370,
                -0.7409511253549589,
                -0.7071067811865475,
                -0.6715589548470184,
                -0.6343932841636454,
                -0.5956993044924334,
                -0.5555702330196020,
                -0.5141027441932217,
                -0.4713967368259977,
                -0.4275550934302819,
                -0.3826834323650897,
                -0.3368898533922199,
                -0.2902846772544622,
                -0.2429801799032639,
                -0.1950903220161282,
                -0.1467304744553616,
                -0.9801714032956065E-01,
                -0.4906767432741801E-01,
                0.0,
                0.4906767432741813E-01,
                0.9801714032956077E-01,
                0.1467304744553617,
                0.1950903220161283,
                0.2429801799032640,
                0.2902846772544623,
                0.3368898533922201,
                0.3826834323650898,
                0.4275550934302822,
                0.4713967368259978,
                0.5141027441932218,
                0.5555702330196023,
                0.5956993044924335,
                0.6343932841636455,
                0.6715589548470184,
                0.7071067811865475,
                0.7409511253549592,
                0.7730104533627370,
                0.8032075314806448,
                0.8314696123025452,
                0.8577286100002721,
                0.8819212643483550,
                0.9039892931234433,
                0.9238795325112867,
                0.9415440651830208,
                0.9569403357322088,
                0.9700312531945440,
                0.9807852804032304,
                0.9891765099647810,
                0.9951847266721969,
                0.9987954562051724,
                1.0000000000000000,
                ])
            self.degree = 65
        else:
            raise ValueError('Illegal Clenshaw-Curtis order')


class NewtonCotesClosed(object):
    '''
    Closed Newton-Cotes formulae.
    <https://en.wikipedia.org/wiki/Newton%E2%80%93Cotes_formulas#Closed_Newton.E2.80.93Cotes_formulae>,
    <http://mathworld.wolfram.com/Newton-CotesFormulas.html>.
    '''
    def __init__(self, index):
        self.points = numpy.linspace(-1.0, 1.0, index+1)
        self.degree = index + 1 if index % 2 == 0 else index

        # Formula (26) from
        # <http://mathworld.wolfram.com/Newton-CotesFormulas.html>.
        # Note that Sympy carries out all operations in rationals, i.e.,
        # _exactly_. Only at the end, the rational is converted into a float.
        n = index
        self.weights = numpy.empty(n+1)
        for r in range(n+1):
            t = sympy.Symbol('t')
            f = 1
            for i in range(n+1):
                if i != r:
                    f *= (t - i)
            alpha = 2 * \
                (-1)**(n-r) * sympy.integrate(f, (t, 0, n)) \
                / (math.factorial(r) * math.factorial(n-r)) \
                / index
            self.weights[r] = alpha

        return


class NewtonCotesOpen(object):
    '''
    Open Newton-Cotes formulae.
    <http://math.stackexchange.com/a/1959071/36678>
    '''
    def __init__(self, index):
        self.points = numpy.linspace(-1.0, 1.0, index+1)[1:-1]
        self.degree = index - 1 if index % 2 == 0 else index - 2
        #
        n = index
        self.weights = numpy.empty(n-1)
        for r in range(1, n):
            t = sympy.Symbol('t')
            f = 1
            for i in range(1, n):
                if i != r:
                    f *= (t - i)
            alpha = 2 * \
                (-1)**(n-r+1) * sympy.integrate(f, (t, 0, n)) \
                / (math.factorial(r-1) * math.factorial(n-1-r)) \
                / n
            self.weights[r-1] = alpha
        return
