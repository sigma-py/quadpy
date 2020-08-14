def get_good_scheme(degree):
    from ._helpers import schemes

    if degree <= 22:
        return {
            0: schemes["dunavant_0"],
            1: schemes["dunavant_0"],
            2: schemes["stroud_c2_3_1"],
            3: schemes["stroud_c2_3_1"],
            4: schemes["burnside"],
            5: schemes["burnside"],
            6: schemes["witherden_vincent_07"],
            7: schemes["witherden_vincent_07"],
            8: schemes["rabinowitz_richter_1"],
            9: schemes["rabinowitz_richter_1"],
            10: schemes["witherden_vincent_11"],
            11: schemes["witherden_vincent_11"],
            12: schemes["witherden_vincent_13"],
            13: schemes["witherden_vincent_13"],
            14: schemes["rabinowitz_richter_6"],
            15: schemes["rabinowitz_richter_6"],
            16: schemes["witherden_vincent_17"],
            17: schemes["witherden_vincent_17"],
            18: schemes["witherden_vincent_19"],
            19: schemes["witherden_vincent_19"],
            20: schemes["witherden_vincent_21"],
            21: schemes["witherden_vincent_21"],
        }[degree]()

    # degree > 50
    assert False
    return
