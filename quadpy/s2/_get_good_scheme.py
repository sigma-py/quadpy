def get_good_scheme(degree):
    from ._helpers import schemes

    if degree <= 19:
        return {
            0: schemes["midpoint"],
            1: schemes["midpoint"],
            2: schemes["albrecht_collatz"],
            3: schemes["albrecht_collatz"],
            4: schemes["mysovskih_1"],
            5: schemes["radon"],
            6: schemes["kim_song_6"],
            7: schemes["kim_song_6"],
            8: schemes["luo_meng_2"],
            9: schemes["luo_meng_2"],
            10: schemes["mysovskih_2"],
            11: schemes["mysovskih_2"],
            12: schemes["cools_haegemans_13_1"],
            13: schemes["cools_haegemans_13_1"],
            14: schemes["kim_song_11"],
            15: schemes["kim_song_11"],
            16: schemes["cools_kim_1"],
            17: schemes["cools_kim_1"],
            18: schemes["kim_song_15"],
            19: schemes["kim_song_15"],
        }[degree]()

    return None
