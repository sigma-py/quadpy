import helpers


def test():
    exp = helpers.get_all_exponents(dim=3, max_degree=3)
    assert exp[0] == [[0, 0, 0]]
    assert exp[1] == [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    assert exp[2] == [[2, 0, 0], [1, 1, 0], [1, 0, 1], [0, 2, 0], [0, 1, 1], [0, 0, 2]]
    assert exp[3] == [
        [3, 0, 0],
        [2, 1, 0],
        [2, 0, 1],
        [1, 2, 0],
        [1, 1, 1],
        [1, 0, 2],
        [0, 3, 0],
        [0, 2, 1],
        [0, 1, 2],
        [0, 0, 3],
    ]
    return
