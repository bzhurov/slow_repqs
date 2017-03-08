config = {
            'name': 'repqs19reqs',
            'cdir': '/data/projects/slow_repqs/auto/',
            'resdir': '/data/projects/',
            'branch': [
                {'name': 'ss_rm', 'sc': '.time', 'c': '.ss', 'icp': ['rm'], 'UZSTOP' : {1: 0.0, 1: 0.01}, 'n_start_repeat' : 200, 'nmx': int(1e8)},
                {'name': 'ss_se', 'sb': 'ss_rm', 'sp': 'UZ', 'c': '.ss', 'UZSTOP': {2: 2000}, 'icp':['se'], 'nmx': int(1e8)},
                {'name': 'lc_rm', 'sb': 'ss_rm', 'sp': 'HB', 'c': '.lc', 'UZSTOP': {1:0.0, 1: 0.01}, 'icp': ['rm', 11], 'nmx': int(1e9)}
                {'name': 'lc_se', 'sb': 'lc_rm', 'sp': 'UZ', 'c': '.lc', 'icp':['se', 11], 'STOP': {'HB', 'LP'} 'nmx': int(1e9)},
                ] 
        }

model = {
        'ode': [
            'rm * G1 - dm * M1',
            'rm * G2 - dm * M2',
            'rm * G3  + rm * omega * ( gt - G3I ) - dm * M3',
            'rm * rmI * GI - dm * MI',
            'rp * M1 - dp * P1 - 2 * kd * P1**2 + 2 * kud * D1',
            'rp * M2 - dp * P2 - 2 * kd * P2**2 + 2 * kud * D2',
            'rp * M3 - dp * P3 - 2 * kd * P3**2 + 2 * kud * D3',
            'rp * MI - dp * PI',
            'rs * PI - ds * S - eta * (S - se) - ka * (pR-R) * S + kua * R',
            'ka * (pR-R) * S - kua * R - 2 * kRd * R**2 + 2 * kRud * DR',
            '-delta * dp * D1 + kd * P1**2 - kud * D1 - kr * D1 * G2 + kur * (gt - G2) - kr * D1 * GI + kur * (gt - GI)',
            '-delta * dp * D2 + kd * P2**2 - kud * D2 - kr * D2 * G3 + kur * (gt - G3)',
            '-delta * dp * D3 + kd * P3**2 - kud * D3 - kr * D3 * G1 + kur * (gt - G1)',
            '-delta * dp * DR + kRd * R**2 - kRud * DR - kRr * DR * G3I + kRur * (gt - G3I)',
            '-kr * D3 * G1 + kur * (gt - G1)',
            '-kr * D1 * G2 + kur * (gt - G2)',
            '-kr * D1 * GI + kur * (gt - GI)',
            '-kr * D2 * G3 + kur * (gt - G3)',
            '-kRr * DR * G3I + kRur * (gt - G3I)'
        ],
        'var': 'M1 M2 M3 MI P1 P2 P3 PI S R D1 D2 D3 DR G1 G2 GI G3 G3I',
        'par': 'rm rp dm dp delta kd kud kr kur gt omega rmI rs ds eta ka kua pR kRd kRud kRr kRur se',
        'G1': 20.0, 'G2': 20.0, 'GI': 20.0, 'G3I': 20.0, 'G3': 20.0,

#         'rm': 0.0, 'rp': 6.0, 'dm': 0.198, 'dp': 0.198,
#         'delta': 0.0, 'kr': 0.720, 'kur': 54, 'kd': 1.5, 'kud': 30.0, 'gt': 20,
#         'se': 0.0, 'omega': 0.0, 'rmI': 2.5, 'dmI': 0.35, 'rpI': 6.94, 'dpI': 0.01,
#         'rs': 0.04, 'ds': 0.001, 
#         'eta': 10, 'kua': 10, 'ka': 0.1, 'kRd': 0.05, 'kRud': 1, 'kRr': 0.05, 'kRur': 10, 'pR': 100

        'rm': 0.0, 'rp': 0.1, 'dm': 0.0033, 'dp': 0.0033,
        'delta': 0.0, 'kr': 0.012, 'kur': 0.9, 'kd': 0.025, 'kud': 0.5, 'gt': 20,
        'se': 0.0, 'omega': 0.5, 'rmI': 0.1,
        # 'rs': 0.0067, 'ds': 1.67e-5, 
        'rs': 0.0067, 'ds': 1.67e-5, 
        'eta': 0.0067, 'kua': 0.166667, 'ka': 0.001667, 'kRd': 0.000833, 'kRud': 0.016667, 'kRr': 0.000833, 'kRur': 0.16667, 'pR': 50

    }
