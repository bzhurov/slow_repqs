kr, kur = 0.012, 0.9 
kd, kud = 0.025, 0.5
ka, kua = 0.05, 0.1
kRr, kRur = 0.000833, 1.6667
kRd, kRud = 0.034, 0.034
#--------------
Cd = kd / kud
Cr = kr / kur
Ca = ka / kua
CRr = kRr / kRur
CRd = kRd / kRud

config = {
            'name': 'repqs9preqs',
            'cdir': '/data/projects/slow_repqs/auto/',
            'resdir': '/data/projects/',
            'branch': [
                {'name': 'ss_rm0', 'sc': '.time', 'c': '.ss', 'icp': ['rm'], 'UZSTOP' : {1: 1.0}, 'n_start_repeat' : 200},
                {'name': 'ss_omega', 'sp': 'UZ', 'sb': 'ss_rm0', 'c': '.ss', 'icp': ['omega'], 'UZSTOP' : {2: 0.0}},

                ] 
        }

model = {
        'ode': [
            'rm * gt / (Cr*Cd*P3**2+1) - dm * M1',
            'rm * gt / (Cr*Cd*P1**2+1) - dm * M2',
            'rm * gt / (Cr*Cd*P2**2+1)  + rm * omega * ( gt - gt / (CRr*CRd*Ca**2*pR**2*S**2+1) ) - dm * M3',
            'rm * rmI * gt / (Cr*Cd*P1**2+1) - dm * MI',
            '(rp * M1 - dp * P1) / (1 + 4*Cd*P1 + 8*gt*Cr*Cd*P1 / (Cr*Cd*P1**2 + 1)**2)',
            '(rp * M2 - dp * P2) / (1 + 4*Cd*P2 + 4*gt*Cr*Cd*P2 / (Cr*Cd*P2**2 + 1)**2)',
            '(rp * M3 - dp * P3) / (1 + 4*Cd*P3 + 4*gt*Cr*Cd*P3 / (Cr*Cd*P3**2 + 1)**2)',
            'rp * MI - dp * PI',
            '(rs * PI - ds * S - eta * (S - se)) / (1 + Ca*pR + 4*gt*CRr*CRd*Ca**2*pR**2*S / (CRr*CRd*Ca**2*pR**2*S**2 + 1)**2)'
        ],
        'var': 'M1 M2 M3 MI P1 P2 P3 PI S',
        'par': 'rm rp dm dp gt omega rmI rs ds eta pR se Cd Cr Ca CRr CRd',

        'rm': 0.0, 'rp': 0.1, 'dm': 0.0033, 'dp': 0.0033,
        'delta': 0.0, 'kr': 0.012, 'kur': 0.9, 'kd': 0.025, 'kud': 0.5, 'gt': 20,
        'se': 0.0, 'omega': 0.0, 'rmI': 0.1,
        'rs': 0.0067, 'ds': 1.67e-5, 
        'eta': 0.0067, 'pR': 5,
        'Cd': Cd, 'Cr': Cr,
        'Ca': Ca, 'CRd': CRd, 'CRr': CRr

    }
