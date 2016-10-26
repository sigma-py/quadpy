'''
Generate code from cartesian coordinates.
'''
import mpmath
# set precision
mpmath.mp.dps = 32

# reference triangle:
t0 = [-1, -1/mpmath.sqrt(3)]
t1 = [+1, -1/mpmath.sqrt(3)]
t2 = [0, 2/mpmath.sqrt(3)]

x = '''
     0.00000000000000000000000000000000D+00, &
      0.00000000000000000000000000000000D+00, &
      0.00000000000000000000000000000000D+00, &
     -0.39685461846296817096661395897363D+00, &
      0.00000000000000000000000000000000D+00, &
     -0.36877346231328110106712038942974D+00, &
      0.00000000000000000000000000000000D+00, &
      0.00000000000000000000000000000000D+00, &
     -0.72443410422579609569939260421083D+00, &
      0.00000000000000000000000000000000D+00/
    '''
y = '''
    -0.56396461592102123502624022053391D+00, &
     -0.47207168193213434618577598142010D+00, &
      0.35411102701218903723318259570554D+00, &
     -0.54446208086261457427222068671521D+00, &
      0.00000000000000000000000000000000D+00, &
     -0.40806649765315498179968814806897D+00, &
     -0.28109188226279360944191830212271D+00, &
      0.76131746182322280252086625652347D+00, &
     -0.53930344496737094791136532996013D+00, &
      0.10684585018901699613124809207464D+01/
    '''
w = '''
   0.65418593445945714253573279005246D-02, &
      0.21571270093488444532084979604155D-01, &
      0.30310770119514528295026716361971D-01, &
      0.23854497740070562467907639645715D-01, &
      0.11323203959116968208057898102191D-01, &
      0.48973694128817658616454407740346D-01, &
      0.30892926213314122881388117756484D-01, &
      0.20335382082811117457514058128897D-01, &
      0.20258422938614600267531787728451D-01, &
      0.52836422050728359852135506640993D-02/
    '''


def convert_to_mpf_list(multiline_str):
    ll = []
    for line in multiline_str.splitlines():
        # properly format
        l = line \
            .replace(',', '') \
            .replace('&', '') \
            .replace('/', '') \
            .replace('D', 'e') \
            .strip()
        if len(l) == 0:
            continue
        ll.append(mpmath.mp.mpf(l))
    return ll

X = convert_to_mpf_list(x)
Y = convert_to_mpf_list(y)
W = convert_to_mpf_list(w)

# generate barycentric coordinate code
XY = [[xx, yy] for xx, yy in zip(X, Y)]
T = mpmath.matrix([
    [t1[0] - t0[0], t2[0] - t0[0]],
    [t1[1] - t0[1], t2[1] - t0[1]],
    ])
tol = 1.0e-10
multiplicities = []
for k, xy in enumerate(XY):
    b = [xy[0] - t0[0], xy[1] - t0[1]]
    sol = mpmath.lu_solve(T, b)
    lmbda = [sol[0], sol[1], 1.0-sol[0]-sol[1]]
    assert abs(sum(lmbda) - 1.0) < tol
    # print('%.15e %.15e %.15e' % (lmbda[0], lmbda[1], lmbda[2]))
    diffs = [
        abs(lmbda[0] - lmbda[1]),
        abs(lmbda[1] - lmbda[2]),
        abs(lmbda[2] - lmbda[0]),
        ]
    if diffs[0] < tol and diffs[1] < tol and diffs[2] < tol:
        print('_s3(),')
        multiplicities.append(1)
    elif diffs[0] < tol:
        print('_s21(%s),' % lmbda[0])
        multiplicities.append(3)
    elif diffs[1] < tol:
        print('_s21(%s),' % lmbda[1])
        multiplicities.append(3)
    elif diffs[2] < tol:
        print('_s21(%s),' % lmbda[2])
        multiplicities.append(3)
    else:
        print('_s111(\n%s,\n%s\n),' % (lmbda[0], lmbda[1]))
        multiplicities.append(6)

# generate weight code
print
for weight, m in zip(W, multiplicities):
    # print('%.15e' % (w / 0.21934566882541541013653648363283 / m))
    print(
        '%s * numpy.ones(%d),'
        % (weight / mpmath.mp.mpf('0.21934566882541541013653648363283') / m, m)
        )
