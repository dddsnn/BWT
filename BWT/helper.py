import re

text = r'''\(\infty\) & \(\infty\) & 2136995 & none (this is standard BW compression) \\
0 & 0 & 2800868 & all \\
0 & 2 & 2270347 & 0x00, \textbackslash n, 0x1a, space, !, \&, ', ), *, +, ``,'',
-, ., 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, :, ;, =, >, ?, A, D, E, K, Q, R, U, V, X, e,
i, s, y \\
0 & 3 & 2160684 & 0x00, \textbackslash n, 0x1a, space, !, \&, ), *, +,
``,'', ., 0, 1, 5, 6, 7, 8, 9, :, ;, =, >, ?, E, K, R, U, V, X, e \\
0 & 4 & 2100693 & 0x00, \textbackslash n, 0x1a, space, !, \&, ), *, +, ``,'', .,
0, 5, 7, :, ;, =, >, ?, E, U, V, X \\
100 & 4 & 2100464 & \textbackslash n, space, !, +, ``,'', ., :, ;, >, ?, E, U \\
0 & 4.5 & 2100162 & 0x00, \textbackslash n, 0x1a, space, !, \&, ), *, +, ``,'',
., 0, 5, :, ;, =, >, ?, X \\
100 & 4.5 & 2100036 & \textbackslash n, space, !, +, ``,'', ., :, ;, >, ? \\
0 & 5 & 2100768 & 0x00, \textbackslash n, 0x1a, space, !, \&, ), *, +,
``,'', ., 0, :, ;, =, >, X \\
100 & 5 & 2100607 & \textbackslash n, space, !, +, ``,'', ., :, ;, > \\
100 & 6 & 2120885 & \textbackslash n, ``,'', . \\'''

def improvement(match):
    natural = 2136995
    size = 6150168
    ratio = natural / size
    res = match.group(1)
    res += ' ({0:.2f})'.format((ratio - (int(match.group(1)) / size)) * 100)
    return res

print(re.sub(r'([0-9]{7})', improvement, text))
