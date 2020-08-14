def get_good_scheme(degree):
    from ._helpers import all_schemes

    if degree <= 22:
        return {
            0: all_schemes["dunavant_0"],
            1: all_schemes["dunavant_0"],
            2: all_schemes["stroud_c2_3_1"],
            3: all_schemes["stroud_c2_3_1"],
            4: all_schemes["burnside"],
            5: all_schemes["burnside"],
            6: all_schemes["witherden_vincent_07"],
            7: all_schemes["witherden_vincent_07"],
            8: all_schemes["rabinowitz_richter_1"],
            9: all_schemes["rabinowitz_richter_1"],
            10: all_schemes["witherden_vincent_11"],
            11: all_schemes["witherden_vincent_11"],
            12: all_schemes["witherden_vincent_13"],
            13: all_schemes["witherden_vincent_13"],
            14: all_schemes["rabinowitz_richter_6"],
            15: all_schemes["rabinowitz_richter_6"],
            16: all_schemes["witherden_vincent_17"],
            17: all_schemes["witherden_vincent_17"],
            18: all_schemes["witherden_vincent_19"],
            19: all_schemes["witherden_vincent_19"],
            20: all_schemes["witherden_vincent_21"],
            21: all_schemes["witherden_vincent_21"],
        }[degree]()

    # degree > 50
    assert False
    return
