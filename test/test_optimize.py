import quadpy


def test_optimize_t2():
    d = {
        "name": "test",
        "domain": "T2",
        "degree": 2,
        "data": {"s2": [[0.33336839947], [0.1669753349]]},
    }
    out, _, _ = quadpy.optimize(d)

    assert abs(out["s2"][0][0] - 1 / 3) < 1.0e-12
    assert abs(out["s2"][1][0] - 1 / 6) < 1.0e-12
