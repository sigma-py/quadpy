import warnings

from ..helpers import article
from ._helpers import C2Scheme, concat, pm, pm2

source = article(
    authors=["C.R. Morrow", "T.N.L. Patterson"],
    title="The Construction of Algebraic Cubature Formulae by the Distribution of Nodes Along Selected Lines",
    journal="SIAM J. Numer. Anal.",
    volume="22",
    number="6",
    pages="1178–1190",
    url="https://doi.org/10.1137/0722071",
)


def morrow_patterson_1():
    warnings.warn("The Morrow-Patterson schemes are only single-precision.")
    weights, points = concat(
        pm2(
            [0.1627661292, 0.2386191861, 0.8611363116],
            [0.3051478053, 0.2386191861, 0.3399810436],
            [0.0140154409, 0.6612093865, 1.112553370],
            [0.2089530426, 0.6612093865, 0.7017715829],
            [-0.1118414644e-03, 0.9324695142, 1.532186887],
            [0.0541846008, 0.9324695142, 0.8829273433],
            [0.1172517331, 0.9324695142, 0.3592258245],
        ),
        pm([0.2755861791, 0.6612093865, 0.0]),
    )
    weights /= 4
    return C2Scheme("Morrow-Patterson 1", weights, points, 11, source)


def morrow_patterson_2():
    warnings.warn("The Morrow-Patterson schemes are only single-precision.")
    weights, points = concat(
        pm2(
            [+0.153974824894e-01, 0.950125098376e-1, 0.968160239507],
            [+0.342239043212e-01, 0.950125098376e-1, 0.836031107326],
            [+0.493728555246e-01, 0.950125098376e-1, 0.613371432700],
            [+0.591743444190e-01, 0.950125098376e-1, 0.324253423403],
            [+0.219687383042e-01, 0.281603550779, 0.950894801068],
            [+0.422241228318e-01, 0.281603550779, 0.771492557752],
            [+0.557059403767e-01, 0.281603550779, 0.500482824502],
            [+0.626761415564e-01, 0.281603550779, 0.173233147858],
            [+0.284719756769e-04, 0.281603550779, 0.114815304093e1],
            [+0.107864593517e-01, 0.458016777657, 0.976231599013],
            [+0.288090801678e-01, 0.458016777657, 0.860064479607],
            [+0.449057826672e-01, 0.458016777657, 0.639857912633],
            [+0.552567682967e-01, 0.458016777657, 0.340595575643],
            [-0.238570094309e-08, 0.458016777657, 0.166045668443e1],
            [+0.152203943218e-02, 0.617876244403, 0.101356615914e1],
            [+0.186940270679e-01, 0.617876244403, 0.936005992420],
            [+0.334515581891e-01, 0.617876244403, 0.760621522239],
            [+0.448929333192e-01, 0.617876244403, 0.496236604823],
            [+0.510354307986e-01, 0.617876244403, 0.172327489281],
            [+0.941251999907e-11, 0.617876244403, 0.217081790217e1],
            [+0.839123482817e-02, 0.755404408355, 0.974824439241],
            [+0.215480412155e-01, 0.755404408355, 0.854811906217],
            [+0.328880560591e-01, 0.755404408355, 0.634285445805],
            [+0.403387733018e-01, 0.755404408355, 0.337501230340],
            [-0.192849256780e-05, 0.755404408355, 0.116365304999e1],
            [-0.211309155229e-12, 0.755404408355, 0.262105563727e1],
            [+0.373264085082e-02, 0.865631202388, 0.984091441674],
            [+0.108083624686e-01, 0.865631202388, 0.912693883002],
            [+0.206089282780e-01, 0.865631202388, 0.746762989660],
            [+0.280439642185e-01, 0.865631202388, 0.488326550840],
            [+0.319642570619e-01, 0.865631202388, 0.169703337599],
            [+0.132307491853e-13, 0.865631202388, 0.298557375978e1],
            [+0.358804444921e-06, 0.865631202388, 0.123267231547e1],
            [+0.377203097373e-02, 0.944575023073, 0.976647356726],
            [+0.694302275949e-02, 0.944575023073, 0.882226287758],
            [+0.651880370677e-02, 0.944575023073, 0.786493198166],
            [+0.151222484324e-01, 0.944575023073, 0.610059075287],
            [+0.194753613784e-01, 0.944575023073, 0.327509702777],
            [-0.155019627860e-14, 0.944575023073, 0.324790845348e1],
            [-0.339816478433e-07, 0.944575023073, 0.131269184417e1],
            [+0.102363309226e-02, 0.989400934992, 0.985756803095],
            [+0.296184898206e-02, 0.989400934992, 0.914705094008],
            [+0.463894901111e-02, 0.989400934992, 0.769393763416],
            [+0.350689331377e-02, 0.989400934992, 0.616728021409],
            [+0.655751993234e-02, 0.989400934992, 0.442872401546],
            [+0.846360999901e-02, 0.989400934992, 0.158216672941],
            [+0.203393765743e-15, 0.989400934992, 0.339722181604e1],
            [+0.508117261709e-08, 0.989400934992, 0.135173206143e1],
        ),
        pm(
            [+0.625640474012e-01, 0.950125098376e-1, 0.0],
            [+0.587968625942e-01, 0.458016777657, 0.0],
            [+0.429295886871e-01, 0.755404408355, 0.0],
            [+0.208441813388e-01, 0.944575023073, 0.0],
        ),
    )
    weights /= 4
    # TODO The article claims degree 31. Check for errors.
    return C2Scheme("Morrow-Patterson 2", weights, points, 25, source)
