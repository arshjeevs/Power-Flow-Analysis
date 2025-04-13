      % nd   V     Ang.  Pg      Qg     PL       QL               Type
    bus_data = [
        1   1.060   0   114.17   -16.9   0       0       0       10     1;
        2   1.045   0   40.00    0      21.7     12.7    -42.0   50.0   2;
        3   1.010   0   0        0      94.2     19.1    23.4    40.0   2;
        4   1.000   0   0        0      47.8     -3.9    0       0      3;
        5   1.000   0   0        0      7.6      1.6     0       0      3;
        6   1.000   0   0        0      11.2     7.5     0       0      3;
        7   1.000   0   0        0      0        0       0       0      3;
        8   1.000   0   0        0      0        0       0       0      3;
        9   1.000   0   0        0      29.5     16.6    0       0      3;
        10  1.000   0   0        0      9.0      5.8     0       0      3;
        11  1.000   0   0        0      3.5      1.8     0       0      3;
        12  1.000   0   0        0      6.1      1.6     0       0      3;
        13  1.000   0   0        0      13.8     5.8     0       0      3;
        14  1.000   0   0        0      14.9     5.0     0       0      3;
    ];
    
      % No  F   T       R        Xl          Xc
    line_data = [
        1   1   2   0.01938   0.05917   0.02640;
        2   1   5   0.05403   0.22304   0.02190;
        3   2   3   0.04699   0.19797   0.01870;
        4   2   4   0.05811   0.17632   0.02460;
        5   2   5   0.05695   0.17388   0.01700;
        6   3   4   0.06701   0.17103   0.01730;
        7   4   5   0.01335   0.04211   0.00640;
        8   4   7   0         0.20912   0;
        9   4   9   0         0.55618   0;
        10  5   6   0         0.25202   0;
        11  6   11  0.09498   0.19890   0;
        12  6   12  0.12291   0.25581   0;
        13  6   13  0.06615   0.13027   0;
        14  7   8   0         0.17615   0;
        15  7   9   0         0.11001   0;
        16  9   10  0.03181   0.08450   0;
        17  9   14  0.12711   0.27038   0;
        18  10  11  0.08205   0.19207   0;
        19  12  13  0.22092   0.19988   0;
        20  13  14  0.17093   0.34802   0;
    ];
        
    tap_data = [
        4   7   0.978;
        4   9   0.969;
        5   6   0.932;
    ];
    
    % Shunt capacitor data from Table A.5
    shunt_data = [
        9   0.19;
    ];
