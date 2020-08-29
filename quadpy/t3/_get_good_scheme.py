def get_good_scheme(degree):
    from ._helpers import schemes

    if degree <= 14:
        return {
            0: schemes["keast_0"],
            1: schemes["keast_0"],
            2: schemes["hammer_marlowe_stroud_1"],
            3: schemes["witherden_vincent_03"],
            4: schemes["liu_vinokur_12"],
            5: schemes["liu_vinokur_12"],
            6: schemes["keast_7"],
            7: schemes["witherden_vincent_07"],
            8: schemes["witherden_vincent_08"],
            9: schemes["witherden_vincent_09"],
            10: schemes["witherden_vincent_10"],
            11: schemes["vioreanu_rokhlin_8"],
            12: schemes["vioreanu_rokhlin_9"],
            13: schemes["vioreanu_rokhlin_9"],
            14: schemes["zhang_cui_liu_2"],
        }[degree]()

    return None
